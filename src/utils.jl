#=
package utils
=#

function GeoUPM_dir()
    joinpath(@__DIR__, "..")
end

function to_grid(E, V)
    # Convert vertices to Ferrite nodes
    nodes = [Ferrite.Node((v[1], v[2], v[3])) for v in V]
    # Create Ferrite cells based on the specified element type
    cells =  [Ferrite.Tetrahedron((e[1], e[2], e[3], e[4])) for e in E]
    return Ferrite.Grid(cells, nodes)
end

