using Revise
using GeoPUPM
using Ferrite
using Geogram
using FileIO
using GLMakie
using GeometryBasics
using Comodo
using PUPM

# --- Setup parameters ---
par = DynamicParams()
sampleSize = 2.0

# Load STL and generate volume mesh
fileName_mesh = joinpath(GeoUPM_dir(), "assets", "stl", "cube_hole.stl")
M = load(fileName_mesh)
F = [TriangleFace{Int64}(f) for f in faces(M)]
V = [Point{3,Float64}(v) for v in coordinates(M)]

n1 = 30000  # remesh point count
F1, V1 = ggremesh(F, V; nb_pts=n1)

# Generate tetrahedral mesh
@time E_tet, V_tet, CE, Fb, CFb_type = tetgenmesh(F1, V1)

# Convert to Ferrite grid
grid = to_grid(E_tet, V_tet)

# --- Calculate bounding box ---
min_x, max_x = extrema([v[1] for v in V_tet])
min_y, max_y = extrema([v[2] for v in V_tet])
min_z, max_z = extrema([v[3] for v in V_tet])

Lx = max_x - min_x
Ly = max_y - min_y
Lz = max_z - min_z

cx_top = (min_x + max_x) / 2
cy_top = (min_y + max_y) / 2
cz_top = max_z
cz_bottom = min_z

# --- Geometric filter for circular ring ---
function in_circular_ring(x, r_inner, r_outer, cx, cy, cz; tol=1e-3)
    z_match = abs(x[3] - cz) < tol
    r2 = (x[1] - cx)^2 + (x[2] - cy)^2
    return z_match && (r2 ≤ r_outer^2) && (r2 ≥ r_inner^2)
end

# --- Boundary creation function ---
function create_boundary(grid, min_x, max_x, min_y, max_y, min_z, max_z)
    cx_top, cy_top = (min_x + max_x) / 2, (min_y + max_y) / 2

    # Top ring parameters
    r_inner = 0.15 * (max_x - min_x)
    r_outer = 0.35 * (max_x - min_x)

    # Define the ring region at the top surface
    addfacetset!(grid, "top_ring", x -> in_circular_ring(x, r_inner, r_outer, cx_top, cy_top, max_z))

    # Bottom corner circles for support
    r_bottom = 0.2 * (max_x - min_x)
    corner_centers = [
        (min_x, min_y),
        (max_x, min_y),
        (min_x, max_y),
        (max_x, max_y)
    ]
    for (i, (cx, cy)) in enumerate(corner_centers)
        addnodeset!(grid, "bottom_corner_circle_$i", x -> (abs(x[3] - min_z) < 1e-3) &&
                                                  ((x[1] - cx)^2 + (x[2] - cy)^2 ≤ r_bottom^2))
    end
end

# --- Cell/facet values setup ---
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron, order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    return CellValues(qr, ip), FacetValues(qr_face, ip)
end

# --- DOF handler setup ---
function create_dofhandler(grid)
    dh = DofHandler(grid)
    add!(dh, :u, Lagrange{RefTetrahedron, 1}()^3)
    close!(dh)
    return dh
end

# --- Boundary condition setup ---
function create_bc(dh, grid)
    ch = ConstraintHandler(dh)
    for i in 1:4
        nodeset = getnodeset(grid, "bottom_corner_circle_$i")
        if !isempty(nodeset)
            dbc = Dirichlet(:u, nodeset, (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
            add!(ch, dbc)
        else
            println("Warning: Node set bottom_corner_circle_$i is empty")
        end
    end
    close!(ch)
    return ch
end

# --- Run setup ---
create_boundary(grid, min_x, max_x, min_y, max_y, min_z, max_z)
par.grid = grid
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh, grid)
par.cell_values, par.facet_values = create_values()

# --- Apply traction load on top ring ---
par.loads = [LoadCondition_3d("traction_load", [0.0, 0.0, 1.0])]
par.Neumann_bc = getfacetset(grid, "top_ring")

# --- Material and optimization parameters ---
par.E0 = 1.0
par.E = fill(par.E0, Ferrite.getncells(grid))
par.ν = 0.3
par.tnele = length(grid.cells)
par.Emin = 1e-4
par.Emax = 1.0
par.ρ0 = 1.0
par.tol = 1e-2
par.γ = 1
par.η = π / (3.5)
par.k = 4
par.vf = 0.5

# --- Output and run ---
file_name = "linear_elasticity_3d"
dir = "/Users/aminalibakhshi/Desktop/vtu_geo/"
remove_vtk_files(dir)

par.max_itr = 300

# Run topology optimization
@time"top time" top_upm_3d(par, file_name, dir)
