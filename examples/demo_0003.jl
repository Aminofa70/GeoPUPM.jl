using Revise
using GeoPUPM
using Ferrite
using Geogram
using FileIO
using GLMakie
using GeometryBasics
using Comodo
using PUPM
# topology optimization to generate a chair from a cube stl file (simple example)
par = DynamicParams()  # Create a dynamic parameter object
sampleSize = 2.0 # Size of the sample form -1 to 1
fileName_mesh = joinpath(GeoUPM_dir(),"assets","stl","cube.stl")
M = load(fileName_mesh)
F = [TriangleFace{Int64}(f) for f in faces(M)]
V = [Point{3,Float64}(v) for v in coordinates(M)]
n = length(V) # Original number of points

# Remeshing the surface 
n1 = 2000
F1, V1 = ggremesh(F, V; nb_pts=n1)

# Generate tetrahedral mesh
# @time E_tet, V_tet, CE, Fb, CFb_type = tetgenmesh(F1, V1)
opts = "-q1.2 -a1e-4"
@time E_tet, V_tet, CE, Fb, CFb_type = tetgenmesh(F1, V1; stringOpt = opts)
## Convert to grid in Ferrite
grid = to_grid(E_tet, V_tet, Ferrite.Tetrahedron)

# Calculate domain bounds dynamically
min_x = minimum([v[1] for v in V_tet])
max_x = maximum([v[1] for v in V_tet])
min_y = minimum([v[2] for v in V_tet])
max_y = maximum([v[2] for v in V_tet])
min_z = minimum([v[3] for v in V_tet])
max_z = maximum([v[3] for v in V_tet])

# Calculate dimensions and centers
Lx = max_x - min_x
Ly = max_y - min_y
Lz = max_z - min_z

cx_top = (min_x + max_x) / 2
cy_top = (min_y + max_y) / 2
cz_top = max_z
cz_bottom = min_z

# Function to check if a point is inside a circle on a given plane
function in_circle(x, r, cx, cy, cz; tol=1e-3)
    return (abs(x[3] - cz) < tol) && ((x[1] - cx)^2 + (x[2] - cy)^2 <= r^2)
end

# Function to define boundary conditions and nodesets
function create_boundary(grid, min_x, max_x, min_y, max_y, min_z, max_z)
    # Parameters for the circle on the top surface
    r_top = 0.3 * (max_x - min_x)  # Adjusted radius for the top circle
    
    cx_top, cy_top = (min_x + max_x) / 2, (min_y + max_y) / 2

    # Parameters for the circles on the bottom surface
    r_bottom = 0.2 * (max_x - min_x)  # Smaller radius for bottom circles
    
    # Define the bottom corner centers dynamically
    corner_centers = [
        (min_x, min_y),
        (max_x, min_y),
        (min_x, max_y),
        (max_x, max_y)
    ]

    # Add nodeset for the top surface inside the circle
    # addnodeset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, max_z))
    addfacetset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, max_z))

    # Add nodesets for each bottom corner circle
    for (i, (cx, cy)) in enumerate(corner_centers)
        addnodeset!(grid, "bottom_corner_circle_$i", x -> in_circle(x, r_bottom, cx, cy, min_z))
    end
end

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron,order}()^dim
    qr = QuadratureRule{RefTetrahedron}(1)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTetrahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    for i in 1:4
        nodeset = getnodeset(grid, "bottom_corner_circle_$i")
        if !isempty(nodeset)
            dbc = Dirichlet(:u, nodeset, (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
            add!(ch, dbc)
        else
            println("Warning: Node set bottom_corner_circle_$i is empty")
        end
    end
    Ferrite.close!(ch)
    return ch
end

# Create boundaries dynamically based on the mesh
create_boundary(grid, min_x, max_x, min_y, max_y, min_z, max_z)
par.grid = grid
par.dh = create_dofhandler(grid)
par.ch = create_bc(par.dh, grid)

par.cell_values, par.facet_values = create_values()
# par.loads = [LoadCondition_3d("nodal_load", [0.0, 0.0, 1.0])]
par.loads = [LoadCondition_3d("traction_load", [0.0, 0.0, 0.01])]

# Material properties
par.E0 = 1.0
par.E = fill(par.E0,Ferrite.getncells(grid))
par.ν = 0.3
par.tnele = length(grid.cells)
# Optimization parameters
par.Emin = 1e-4
par.Emax = 1.0
par.ρ0 = 1.0
par.tol = 1e-2
par.γ = 1
par.η = π / (3.5)
par.k = 4
par.vf = 0.5

# Neumann BC
# par.Neumann_bc = par.Neumann_bc = Ferrite.getnodeset(grid, "top_circle")  # Nodes on the edge
par.Neumann_bc = Ferrite.getfacetset(grid, "top_circle")  # Nodes on the edge
file_name = "linear_elasticity_3d"
dir = "/Users/aminalibakhshi/Desktop/vtu_geo/"
par.max_itr = 300
remove_vtk_files(dir) # optional
#fem example
# Run topology optimization
@time"top time" top_upm_3d(par, file_name, dir)