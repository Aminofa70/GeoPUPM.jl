#=
package utils

=#

function GeoUPM_dir()
    joinpath(@__DIR__, "..")
end


function to_grid(E, V, type::Type)
    # Convert vertices to Ferrite nodes
    nodes = [Ferrite.Node((v[1], v[2], v[3])) for v in V]

    # Create Ferrite cells based on the specified element type
    cells = if type == Ferrite.Tetrahedron
        [Ferrite.Tetrahedron((e[1], e[2], e[3], e[4])) for e in E]
    elseif type == Ferrite.Hexahedron
        [Ferrite.Hexahedron((e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8])) for e in E]
    else
        error("Unsupported element type. Use Ferrite.Tetrahedron or Ferrite.Hexahedron.")
    end

    return Ferrite.Grid(cells, nodes)
end

struct VoxelGrid
    voxels::Array{Int8, 3}
    scale::Float64
    offset::Vector{Float64}
end

distance(p₁, p₂) = norm(p₁ - p₂)

function split_edge!(new, p₁, p₂, length::Float64, iter::Int, max_iter::Int=20)
    if iter <= max_iter && distance(p₁, p₂) > length
        new_p = (p₁ + p₂) / 2
        push!(new, new_p)
        split_edge!(new, p₁, new_p, length, iter + 1, max_iter)
        split_edge!(new, p₂, new_p, length, iter + 1, max_iter)
    end
end

# ------------------------------------------------------
function fill_voxel_interior!(voxels::Array{Int8,3})
    nx, ny, nz = size(voxels)
    visited = falses(nx, ny, nz)
    queue = Vector{NTuple{3,Int}}()

    # 1) enqueue all boundary‐zero voxels
    for i in 1:nx, j in 1:ny
        if voxels[i,j,1] == 0
            push!(queue, (i,j,1)); visited[i,j,1] = true
        end
        if voxels[i,j,nz] == 0
            push!(queue, (i,j,nz)); visited[i,j,nz] = true
        end
    end
    for i in 1:nx, k in 1:nz
        if voxels[i,1,k] == 0
            push!(queue, (i,1,k)); visited[i,1,k] = true
        end
        if voxels[i,ny,k] == 0
            push!(queue, (i,ny,k)); visited[i,ny,k] = true
        end
    end
    for j in 1:ny, k in 1:nz
        if voxels[1,j,k] == 0
            push!(queue, (1,j,k)); visited[1,j,k] = true
        end
        if voxels[nx,j,k] == 0
            push!(queue, (nx,j,k)); visited[nx,j,k] = true
        end
    end

    # 2) BFS 6‐connected flood
    nbrs = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
    while !isempty(queue)
        i,j,k = pop!(queue)
        for (di,dj,dk) in nbrs
            ii, jj, kk = i+di, j+dj, k+dk
            if 1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny && 1 ≤ kk ≤ nz &&
               !visited[ii,jj,kk] && voxels[ii,jj,kk] == 0
                visited[ii,jj,kk] = true
                push!(queue, (ii,jj,kk))
            end
        end
    end

    # 3) fill any remaining zeros (interior)
    for i in 1:nx, j in 1:ny, k in 1:nz
        if voxels[i,j,k] == 0 && !visited[i,j,k]
            voxels[i,j,k] = 1
        end
    end

    return voxels
end

# ------------------------------------------------------
function voxelize(mesh::GeometryBasics.Mesh, pitch::Float64; max_iter=10, edge_factor=2.0)
    max_edge = pitch / edge_factor
    pts = Point3[]

    # sample edges + triangle interiors
    for tri in mesh
        p₁, p₂, p₃ = tri
        split_edge!(pts, p₁, p₂, max_edge, 0, max_iter)
        split_edge!(pts, p₂, p₃, max_edge, 0, max_iter)
        split_edge!(pts, p₃, p₁, max_edge, 0, max_iter)

        # barycentric grid on face
        L = maximum((distance(p₁,p₂), distance(p₂,p₃), distance(p₃,p₁)))
        n = ceil(Int, L / max_edge)
        for i in 0:n, j in 0:(n-i)
            u = i/n; v = j/n; w = 1 - u - v
            push!(pts, u*p₂ + v*p₃ + w*p₁)
        end
    end

    # bounds → voxel indices
    vmin = Point3( typemax(Float32) )
    vmax = Point3( typemin(Float32) )
    for p in pts
        vmin = min.(p, vmin)
        vmax = max.(p, vmax)
    end
    widths = vmax .- vmin
    trans = LinearMap(UniformScaling(1/pitch)) ∘ Translation(-vmin)
    pts = map(trans, pts)

    grid_size = ceil.(Int, widths ./ pitch) .+ 1
    voxels = zeros(Int8, grid_size...)

    Threads.@threads for p in pts
        idx = floor.(Int, p .+ 1)
        if all(1 .<= idx .<= size(voxels))
            voxels[idx...] = 1
        end
    end

    fill_voxel_interior!(voxels)
    return VoxelGrid(voxels, pitch, collect(vmin))
end

# ------------------------------------------------------
function voxelgrid_to_hex8(grid::VoxelGrid)
    nx, ny, nz = size(grid.voxels)
    pitch, offset = grid.scale, grid.offset

    point_dict = Dict{Tuple{Int,Int,Int}, Int}()
    points = Point3f0[]
    cells  = NTuple{8,Int}[]

    corner_offsets = [
        (0,0,0),(1,0,0),(1,1,0),(0,1,0),
        (0,0,1),(1,0,1),(1,1,1),(0,1,1)
    ]

    for i in 1:nx-1, j in 1:ny-1, k in 1:nz-1
        if grid.voxels[i,j,k] == 1
            inds = Int[]
            for (dx,dy,dz) in corner_offsets
                key = (i+dx, j+dy, k+dz)
                if !haskey(point_dict, key)
                    # ▶ key .- 1 ensures index=1 → world = offset + 0*pitch
                    coords = offset .+ pitch .* (key .- 1)
                    x,y,z = coords
                    push!(points, Point3f0(x,y,z))
                    point_dict[key] = length(points)
                end
                push!(inds, point_dict[key])
            end
            push!(cells, Tuple(inds))
        end
    end

    return points, cells
end

