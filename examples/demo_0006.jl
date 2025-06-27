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
par           = DynamicParams()  # Create a dynamic parameter object
fileName_mesh = joinpath(GeoUPM_dir(),"assets","stl","cube.stl")
mesh          = load(fileName_mesh)

grid_vox      = voxelize(mesh, 0.005)
V, E          = voxelgrid_to_hex8(grid_vox)
grid          = to_grid(E, V, Ferrite.Hexahedron)
