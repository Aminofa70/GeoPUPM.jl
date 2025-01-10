using Revise
using GeoPUPM
using Ferrite
using Geogram
using FileIO
using GLMakie
using GeometryBasics
using Comodo
#=
this demo test the boundary condistions of the ggremesh function

=#

fileName_mesh = joinpath(GeoUPM_dir(),"assets","stl","cube.stl")
M = load(fileName_mesh)
F = [TriangleFace{Int64}(f) for f in faces(M)]
V = [Point{3,Float64}(v) for v in coordinates(M)]
n = length(V) # Original number of points
# Remeshing the surface 
n1 = 10000
F1,V1 = ggremesh(F,V; nb_pts=n1)

# Generate tetrahedral mesh
E_tet, V_tet, CE, Fb, CFb_type = tetgenmesh(F1, V1)
## convert to grid in Ferrite
grid = to_grid(E_tet, V_tet)
min_x = minimum([v[1] for v in V_tet])
max_x = maximum([v[1] for v in V_tet])
min_y = minimum([v[2] for v in V_tet])
max_y = maximum([v[2] for v in V_tet])
min_z = minimum([v[3] for v in V_tet])
max_z = maximum([v[3] for v in V_tet])

### define facet sets
addfacetset!(grid, "top", x -> abs(x[3] - (max_z)) < 1e-2)
addfacetset!(grid, "bottom", x -> abs(x[3] - (min_z)) < 1e-2)
addfacetset!(grid, "front", x -> abs(x[2] - (min_y)) < 1e-2)
addfacetset!(grid, "back", x -> abs(x[2] - (max_y)) < 1e-2)
addfacetset!(grid, "left", x -> abs(x[1] - (min_x)) < 1e-2)
addfacetset!(grid, "right", x -> abs(x[1] - (max_x)) < 1e-2)

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron,order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
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
    # dbc = Dirichlet(:u, getfacetset(grid, "my_top"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3])
    dbc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3])
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
grid = grid
dh = create_dofhandler(grid)
ch = create_bc(dh, grid)
VTKGridFile("boundary-conditions", dh) do vtk
    Ferrite.write_constraints(vtk,ch)
end