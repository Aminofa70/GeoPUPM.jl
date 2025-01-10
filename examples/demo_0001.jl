using Revise
using GeoPUPM
using Ferrite
using Geogram
using FileIO
using GLMakie
using GeometryBasics
using Comodo
#=
this demo shows how to remesh a surface using the ggremesh function
=#

fileName_mesh = joinpath(GeoUPM_dir(),"assets","stl","cube.stl")
M = load(fileName_mesh)
F = [TriangleFace{Int64}(f) for f in faces(M)]
V = [Point{3,Float64}(v) for v in coordinates(M)]
n = length(V) # Original number of points
# Remeshing the surface 
n1 = 8000 
n2 = 4000
n3 = 1000
F1,V1 = ggremesh(F,V; nb_pts=n1)
F2,V2 = ggremesh(F,V; nb_pts=n2, remesh_anisotropy=0.0, remesh_gradation = 1.0, pre_max_hole_area=100, pre_max_hole_edges=0, post_max_hole_area=100, post_max_hole_edges=0, quiet=0, suppress = true)
F3,V3 = ggremesh(F,V; nb_pts=n3)

## VISUALISATION

strokeWidth1 = 0.5

fig = Figure(size=(1200,1200))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original, $n points")
hp1=poly!(ax1,GeometryBasics.Mesh(V,F),strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Remeshed, $n1 points")
hp2 = poly!(ax2,GeometryBasics.Mesh(V1,F1),strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Remeshed, $n2 points, anisotropic")
hp3 = poly!(ax3,GeometryBasics.Mesh(V2,F2),strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Remeshed, $n3 points")
hp4 = poly!(ax4,GeometryBasics.Mesh(V3,F3),strokewidth=strokeWidth1,color=:white, shading = FastShading)

fig