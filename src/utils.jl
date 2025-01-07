#=
package utils

=#

function to_grid(E, V)
    # Create tetrahedral elements
    cells = [Ferrite.Tetrahedron((e[1], e[2], e[3], e[4])) for e in E]

    # Convert vertices (V) to Ferrite Nodes
    nodes = [Ferrite.Node((v[1], v[2], v[3])) for v in V]
    # Return the grid
    return Grid(cells, nodes)

end