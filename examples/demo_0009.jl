using Revise
using GeoPUPM
using Ferrite
using Geogram
using FileIO
using GLMakie
using GeometryBasics
using Comodo
using PUPM
using Rotations  # For RotX, RotY, RotZ femur geometry
par = DynamicParams()
# --- Load STL file ---
fileName_mesh = joinpath(GeoUPM_dir(), "assets", "stl", "femur_iso.stl")
M = load(fileName_mesh)

# Extract faces and vertices
F = [TriangleFace{Int64}(f) for f in faces(M)]
V = [Point{3, Float64}(v) for v in coordinates(M)]

# --- STEP 1: Scale and Rotate ---
scaleFactorSize = 1.0  # Optional user-defined global scale
# V = [Point{3, Float64}(v .* 1000.0 .* scaleFactorSize) for v in V]  # convert to mm
V = [Point{3, Float64}(v .* 1.0 .* scaleFactorSize) for v in V]  # convert to mm

# Apply rotations using RotX, RotY, RotZ (equivalent to Euler sequence)
Rz1 = RotZ(0.065π)
Rx  = RotX(-0.5π)
Rz2 = RotZ(0.36π)
RflipZ = RotX(π)  # Flip Z-up to Z-down

# Final rotation matrix
R = RflipZ * Rz2 * Rx * Rz1

# Apply rotation to all points
V = [Point{3, Float64}(R * v) for v in V]

# --- Remesh the surface ---
n1 = 1000  # target number of surface points
F1, V1 = ggremesh(F, V; nb_pts = n1)

# --- Generate tetrahedral volume mesh ---
@time E_tet, V_tet, CE, Fb, CFb_type = tetgenmesh(F1, V1)

# --- Convert to Ferrite grid ---
grid = to_grid(E_tet, V_tet)

# --- Bounding box info for boundary tagging ---
min_x = minimum(v[1] for v in V_tet)
max_x = maximum(v[1] for v in V_tet)
min_y = minimum(v[2] for v in V_tet)
max_y = maximum(v[2] for v in V_tet)
min_z = minimum(v[3] for v in V_tet)
max_z = maximum(v[3] for v in V_tet)

# tol = 0.13 * minimum([max_x - min_x, max_y - min_y, max_z - min_z])
# tol_2 = 0.2 * minimum([max_x - min_x, max_y - min_y, max_z - min_z])

# # --- Define boundary facet sets ---
# addfacetset!(grid, "top",    x -> abs(x[3] - max_z) < tol)
# addfacetset!(grid, "bottom", x -> abs(x[3] - min_z) < tol)
# addfacetset!(grid, "front",  x -> abs(x[2] - min_y) < tol)
# addfacetset!(grid, "back",   x -> abs(x[2] - max_y) < tol)
# addfacetset!(grid, "left",   x -> abs(x[1] - min_x) < tol)
# addfacetset!(grid, "right",  x -> abs(x[1] - max_x) < tol)
tol = 0.03 * minimum([max_x - min_x, max_y - min_y, max_z - min_z])
tol_b= 0.10 * minimum([max_x - min_x, max_y - min_y, max_z - min_z])

tol_2 = 0.18 * minimum([max_x - min_x, max_y - min_y, max_z - min_z])

# --- Define boundary facet sets ---
addfacetset!(grid, "top",    x -> abs(x[3] - max_z) < tol)
addfacetset!(grid, "bottom", x -> abs(x[3] - min_z) < tol_b)
addfacetset!(grid, "front",  x -> abs(x[2] - min_y) < tol)
addfacetset!(grid, "back",   x -> abs(x[2] - max_y) < tol)
addfacetset!(grid, "left",   x -> abs(x[1] - min_x) < tol)
addfacetset!(grid, "right",  x -> abs(x[1] - max_x) < tol)

addfacetset!(grid, "top_2", x ->
    abs(x[3] - max_z) < tol_2 &&     # at the back (Y direction)
    abs(x[1] - min_x) < tol_2        # strictly near left edge in X
)

# --- Define cell and facet values ---
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron, order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# --- Create DOF handler ---
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end

# --- Create Dirichlet BCs ---
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end

# --- Run setup ---
par.grid = grid
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh, grid)
par.cell_values, par.facet_values = create_values()


par.load_applications = [
    LoadApplication_3d(
        LoadCondition_3d_geo("traction_load", [0.0, 0.0,-1]),
        getfacetset(grid, "top")
    ),
    LoadApplication_3d(
        LoadCondition_3d_geo("traction_load", [0.0, 0.0, 1]),
        getfacetset(grid, "top_2")
    )
]

# par.load_applications = [
#     LoadApplication_3d(
#         LoadCondition_3d_geo("traction_load", [0.0, 0.0, 1]),
#         getfacetset(grid, "top")
#     )
# ]

# --- Material and optimization parameters ---
par.E0 = 1.0
par.E = fill(par.E0, Ferrite.getncells(grid))
par.ν = 0.35
par.tnele = length(grid.cells)
par.Emin = 1e-4
par.Emax = 1.0
par.ρ0 = 1.0
par.tol = 1e-2
par.γ = 1
par.η = π / (4.0)
par.k = 10
par.vf = 0.5

# --- Output and run ---
file_name = "linear_elasticity_3d"
dir = "/Users/aminalibakhshi/Desktop/vtu_geo/"
remove_vtk_files(dir)

par.max_itr = 300

# Run topology optimization
# @time"top time" top_upm_3d(par, file_name, dir)
@time"top time" top_upm_3d_geo(par, file_name, dir)

