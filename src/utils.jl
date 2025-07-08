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

# Struct for load conditions in 3D
struct LoadCondition_3d_geo
    load_type::String
    load_data::Union{Nothing, Vector{Float64}, Float64}
end

# Struct for associating a load condition with its Neumann boundary set
struct LoadApplication_3d
    load_condition::LoadCondition_3d_geo
    neumann_bc::Any  # Replace `Any` with concrete type if known
end

# Struct for finite element solver results in 3D
struct FEMSolver_3d_geo
    u::Vector{Float64}
    compliance::Float64
    σ::Any
    ε::Any
    U::Vector{Float64}
    H::Vector{Float64}
    E_node::Vector{Float64}
end

function fem_solver_3d_geo(par::DynamicParams)
    # Extract parameters
    grid = par.grid
    cell_values = par.cell_values
    facet_values = par.facet_values
    dh = par.dh
    ch = par.ch
    E = par.E
    ν = par.ν

    # Assemble global stiffness matrix
    K = allocate_matrix(dh)
    assemble_global_3d!(K, dh, cell_values, E, ν)

    # Initialize external force vector
    f_ext = zeros(ndofs(dh))

    # Apply multiple loads
    if !isempty(par.load_applications)
        for app in par.load_applications
            load = app.load_condition
            Neumann_bc = app.neumann_bc

            @assert load.load_data !== nothing "Load data must not be nothing"

            if load.load_type == "traction_load"
                @assert load.load_data isa Vector{Float64}
                assemble_external_forces_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)

            elseif load.load_type == "nodal_load"
                @assert load.load_data isa Vector{Float64}
                apply_nodal_force_3d!(grid, Neumann_bc, load.load_data, f_ext, dh)

            elseif load.load_type == "pressure_load"
                @assert load.load_data isa Float64 || load.load_data isa Vector{Float64}
                assemble_external_pressure_3d!(f_ext, dh, Neumann_bc, facet_values, load.load_data)

            else
                error("Unknown load type: $(load.load_type)")
            end
        end
    else
        println("⚠️ Warning: No load applications provided. Proceeding without external forces.")
    end

    # Apply constraints and solve the system
    apply!(K, f_ext, ch)
    u = K \ f_ext
    c = 0.5 * dot(f_ext, u)

    # Post-processing
    _, σ = calculate_stresses_3d(grid, dh, cell_values, u, E, ν)
    _, ε = calculate_strains_3d(grid, dh, cell_values, u)
    U = calculate_strain_energy_3d(grid, dh, cell_values, u, E, ν)
    H = calculate_H_3d(grid, dh, cell_values, u, E, ν)
    E_node = compute_nodal_data_3D(grid, E)

    return FEMSolver_3d_geo(u, c, σ, ε, U, H, E_node)
end

function top_upm_3d_geo(par::DynamicParams, name_of_file::String, directory::String)
    grid = par.grid
    dh = par.dh
    E = par.E
    # nx = par.nx ; ny = par.ny ; nz = par.nz
    tnele = par.tnele
    E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
    k = par.k ; γ = par.γ ; volfrac = par.vf ; η = par.η ; ρ0 = par.ρ0
    max_itr = par.max_itr ; tol = par.tol

    loop = 1

    # Initial FEM solve
    fem = fem_solver_3d_geo(par)
    compliance = fem.compliance
    H = fem.H
    W_tot = sum(fem.U)

    strain_energy_vector = [W_tot, W_tot * 10]
    A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
    println("Iter $loop: C = $compliance, ΔR = $A")

    # Iterative optimization loop
    while abs(A) > tol && loop <= max_itr
        # FEM solve with current parameters
        fem = fem_solver_3d_geo(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)

        # Material update routines
        Enew = update_upm(k, E, H, Emax, Emin)
        ρ = transfer_to_density(Enew, E0, ρ0, γ)
        ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
        Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)

        # Update par to reflect the new E distribution
        par.E = Enew_frac
        E = Enew_frac

        # FEM solve after material update
        fem = fem_solver_3d_geo(par)
        compliance = fem.compliance
        H = fem.H
        W_tot = sum(fem.U)
        u = fem.u
        σ = fem.σ
        ε = fem.ε
        E_node = fem.E_node
        # Create iteration-specific filename
        loop_str = @sprintf("%03d", loop)
        filename_it = string(name_of_file, "_", loop_str)
        full_path = joinpath(directory, filename_it)

        VTKGridFile(full_path, dh) do vtk
            write_solution(vtk, dh, u)
            # Write 3D stress components
            for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
                write_cell_data(vtk, σ[j], "stress_" * key)
            end
            # Write 3D strain components
            for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
                write_cell_data(vtk, ε[j], "strain_" * key)
            end
            write_cell_data(vtk, E, "Young's modulus")
            write_node_data(vtk, E_node, "Nodal Young's modulus")
            Ferrite.write_cellset(vtk, grid)
        end

        # Update strain energy vector and A
        strain_energy_vector[1] = strain_energy_vector[2]
        strain_energy_vector[2] = W_tot
        A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

        loop += 1
        println("Iter $loop: C = $compliance, ΔR = $A")
    end

    # Handle termination
    if loop > max_itr
        compliance = -1
    end

    # Final solve and write results
    fem = fem_solver_3d_geo(par)
    compliance = fem.compliance
    u = fem.u
    σ = fem.σ
    ε = fem.ε
    E_node = fem.E_node
    full_path = joinpath(directory, name_of_file)

    VTKGridFile(full_path, dh) do vtk
        write_solution(vtk, dh, u)
        # Write 3D stress components
        for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
            write_cell_data(vtk, σ[j], "stress_" * key)
        end
        # Write 3D strain components
        for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
            write_cell_data(vtk, ε[j], "strain_" * key)
        end
        write_cell_data(vtk, E, "Young's modulus")
        write_node_data(vtk, E_node, "Nodal Young's modulus")
        Ferrite.write_cellset(vtk, grid)
    end
end






###################################################

# function calculate_element_volumes(grid, dh, cv)
#     element_volumes = zeros(getncells(grid))

#     for cell in CellIterator(dh)
#         reinit!(cv, cell)
#         cell_volume = 0.0

#         for q_point in 1:getnquadpoints(cv)
#             dΩ = getdetJdV(cv, q_point)
#             cell_volume += dΩ
#         end

#         element_volumes[cellid(cell)] = cell_volume
#     end

#     return element_volumes
# end

# function density_filter_to_vf!(density, vf, volumes, eta)
#     rhomin, rhomax = 0.01, 1.0
#     total_volume = sum(volumes)

#     function transform(rholoc, rhotr, eta, rhomin, rhomax)
#         if rholoc < rhotr
#             return rhomin
#         elseif rholoc > rhotr + 1.0 / tan(eta)
#             return rhomax
#         else
#             return tan(eta) * (rholoc - rhotr)
#         end
#     end

#     rho_vf0_bound = -1.0 / tan(eta)
#     rho_vf1_bound = 1.0
#     error = 10.0
#     rhotr = 0.0

#     while error > 0.001
#         rhotr = (rho_vf1_bound + rho_vf0_bound) / 2

#         sumdmin = sum(transform(density[i], rho_vf1_bound, eta, rhomin, rhomax) * volumes[i] for i in eachindex(density))
#         sumdmax = sum(transform(density[i], rho_vf0_bound, eta, rhomin, rhomax) * volumes[i] for i in eachindex(density))
#         sumdmid = sum(transform(density[i], rhotr, eta, rhomin, rhomax) * volumes[i] for i in eachindex(density))

#         if (sumdmin - vf * total_volume) / (sumdmid - vf * total_volume) > 0
#             rho_vf1_bound = rhotr
#         elseif (sumdmax - vf * total_volume) / (sumdmid - vf * total_volume) > 0
#             rho_vf0_bound = rhotr
#         else
#             println("Error: Bounds mismatch: ", sumdmax, sumdmin, sumdmid, vf * total_volume)
#         end

#         error = abs(vf * total_volume - sumdmid)
#     end

#     for i in eachindex(density)
#         density[i] = transform(density[i], rhotr, eta, rhomin, rhomax)
#     end

#     return density
# end

# function top_upm__geo_3d(par::DynamicParams, name_of_file::String, directory::String)
#     grid = par.grid
#     dh = par.dh
#     E = par.E
#     # nx = par.nx ; ny = par.ny ; nz = par.nz
#     tnele = par.tnele
#     E0 = par.E0 ; Emin = par.Emin ; Emax = par.Emax
#     k = par.k ; γ = par.γ ; volfrac = par.vf ; η = par.η ; ρ0 = par.ρ0
#     max_itr = par.max_itr ; tol = par.tol
#     cv = par.cell_values
#     loop = 1

#     # Initial FEM solve
#     fem = fem_solver_3d(par)
#     compliance = fem.compliance
#     H = fem.H
#     W_tot = sum(fem.U)

#     strain_energy_vector = [W_tot, W_tot * 10]
#     A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]
#     println("Iter $loop: C = $compliance, ΔR = $A")

#     # Iterative optimization loop
#     while abs(A) > tol && loop <= max_itr
#         # FEM solve with current parameters
#         fem = fem_solver_3d(par)
#         compliance = fem.compliance
#         H = fem.H
#         W_tot = sum(fem.U)
#         volumes = calculate_element_volumes(grid, dh, cv)

#         # Material update routines
#         Enew = update_upm(k, E, H, Emax, Emin)
#         ρ = transfer_to_density(Enew, E0, ρ0, γ)
#         volumes = calculate_element_volumes(grid, dh, cv)
#         ρnew = density_filter_to_vf!(ρ, volfrac, volumes, η)
#         # ρnew = filter_density_to_vf!(ρ, volfrac, tnele, η)
#         Enew_frac = transfer_to_young(ρnew, E0, ρ0, γ, Emin, Emax)

#         # Update par to reflect the new E distribution
#         par.E = Enew_frac
#         E = Enew_frac

#         # FEM solve after material update
#         fem = fem_solver_3d(par)
#         compliance = fem.compliance
#         H = fem.H
#         W_tot = sum(fem.U)
#         u = fem.u
#         σ = fem.σ
#         ε = fem.ε
#         E_node = fem.E_node
#         # Create iteration-specific filename
#         loop_str = @sprintf("%03d", loop)
#         filename_it = string(name_of_file, "_", loop_str)
#         full_path = joinpath(directory, filename_it)

#         VTKGridFile(full_path, dh) do vtk
#             write_solution(vtk, dh, u)
#             # Write 3D stress components
#             for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
#                 write_cell_data(vtk, σ[j], "stress_" * key)
#             end
#             # Write 3D strain components
#             for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
#                 write_cell_data(vtk, ε[j], "strain_" * key)
#             end
#             write_cell_data(vtk, E, "Young's modulus")
#             write_node_data(vtk, E_node, "Nodal Young's modulus")
#             Ferrite.write_cellset(vtk, grid)
#         end

#         # Update strain energy vector and A
#         strain_energy_vector[1] = strain_energy_vector[2]
#         strain_energy_vector[2] = W_tot
#         A = (strain_energy_vector[2] - strain_energy_vector[1]) / strain_energy_vector[1]

#         loop += 1
#         println("Iter $loop: C = $compliance, ΔR = $A")
#     end

#     # Handle termination
#     if loop > max_itr
#         compliance = -1
#     end

#     # Final solve and write results
#     fem = fem_solver_3d(par)
#     compliance = fem.compliance
#     u = fem.u
#     σ = fem.σ
#     ε = fem.ε
#     E_node = fem.E_node
#     full_path = joinpath(directory, name_of_file)

#     VTKGridFile(full_path, dh) do vtk
#         write_solution(vtk, dh, u)
#         # Write 3D stress components
#         for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
#             write_cell_data(vtk, σ[j], "stress_" * key)
#         end
#         # Write 3D strain components
#         for (j, key) in enumerate(("11", "22", "33", "12", "23", "13"))
#             write_cell_data(vtk, ε[j], "strain_" * key)
#         end
#         write_cell_data(vtk, E, "Young's modulus")
#         write_node_data(vtk, E_node, "Nodal Young's modulus")
#         Ferrite.write_cellset(vtk, grid)
#     end
# end
