#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Active forces
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abstract type AbstractActiveForce end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Units
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct UnitScaling
    length::Float64   # e.g., in units of R
    time::Float64     # e.g., in units of R^2 / D
    energy::Float64   # e.g., in units of k_B * T
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Particle
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractParticle
export BrownianParticle

abstract type AbstractParticle end

"""
Stores properties of a particle

    Particle(; ptype, x, y, θ, R, γₜ, γᵣ, Dₜ, Dᵣ, active_force)

# Arguments
- `ptype::Symbol = :particle`: particle type
- `x::Float64`: x position
- `y::Float64`: y position
- `z::Float64`: z position
- `θ::Float64 = 2 * pi * rand()`: in-plane orientation
- `R::Float64`: radius
- `γₜ::Float64 = -1.0`: translational friction.  If `γ_trans < 0`, use `γ_trans = 6 * pi * 8.9e-4 * R`.
- `γᵣ::Float64 = -1.0`: rotational friction. If `γ_rot < 0`, use `γ_rot = 8 * pi * 8.9e-4 * R^3`.
- `Dₜ::Float64`: translational diffusion
- `Dᵣ::Float64`: rotational diffusion
- `active_force::Union{AbstractActiveForce, Nothing} = nothing`: active force.
-  units::UnitScaling = UnitScaling(1.0, 1.0, 1.0): scaling of units for length, time, and energy.
"""
mutable struct BrownianParticle <: AbstractParticle
    ptype::Symbol # particle type

    x_d::SVector{3, Float64}
    f_d::SVector{3, Float64}

    R::Float64 # particle radius
    kᵦT::Float64 # Boltzmann constant times temperature

    γₜ::Float64 # translational friction
    γᵣ::Float64 # rotational friction
    Dₜ::Float64 # translational diffusivity
    Dᵣ::Float64 # rotational diffusivity
end

# Outer constructor with default values
function BrownianParticle(; ptype::Symbol = :particle,
    x_d::SVector{3, Float64}, f_d::SVector{3, Float64} = SVector(0.0, 0.0, 0.0),
    R::Float64 = 1.0, kᵦT::Float64 = 4e-21,
    γₜ::Float64 = -1.0, γᵣ::Float64 = -1.0, Dₜ::Float64 = 1.0, Dᵣ::Float64 = -1.0)


    if γₜ < 0.0
        γₜ = kᵦT / Dₜ
    end
    if Dᵣ < 0.0
        Dᵣ = 3Dₜ / R^2
    end
    if γᵣ < 0.0
        γᵣ = kᵦT / Dᵣ
    end

    return BrownianParticle(ptype, x_d, f_d, R, kᵦT, γₜ, γᵣ, Dₜ, Dᵣ)

end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cell list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export CellList

"""
Stores properties of a cell list

    CellList(; particles, L_x, L_y, cutoff)

Creates a cell list for `particles` in region/simulation with dimensions
`L_x`, `L_y`.  `cutoff` gives the approximate size of the cells, the actual size
of which is chosen so that an integer number of cells fit in each dimension.
"""
struct CellList{A <: AbstractParticle}
    particles::Array{A, 1}

    num_cells_d::Vector{Int64}
    cell_spacing_d::Vector{Float64}

    δr::Float64
    r_cut::Float64
    r_skin::Float64

    cells_to_check::Array{Array{Int64, 1}, 1}
    cells_to_particles::Array{Array{Int64, 1}, 1}
    particles_to_cells::Array{Int64, 1}
end

# Outer constructor with inferred particle type
function CellList(; particles::Array{A, 1}, L_d::SVector{3, Float64}, δr::Float64) where A <: AbstractParticle
    # Calculate properties based on `particles` and `L_d`
    δr *= particles[1].R
    r_cut = 2^(1/6) * particles[1].R  # Assuming (wrongly) it is the smallest particle radius
    r_skin = r_cut + δr

    num_cells_d = [floor(Int64, L_d[i] / r_skin) for i in 1:3]
    cell_spacing_d = [L_d[i] / num_cells_d[i] for i in 1:3]

    # Initialize the cell structure arrays
    cells_to_check = linked_grids(num_cells_d)
    cells_to_particles = [Int[] for _ in 1:prod(num_cells_d)]
    cells_to_particles[1] = collect(1:length(particles))
    particles_to_cells = ones(Int, length(particles))

    return CellList{A}(particles, num_cells_d, cell_spacing_d, δr, r_cut, r_skin, cells_to_check, cells_to_particles, particles_to_cells)
end



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Verlet list
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# export VerletList

# """
# Stores properties of a cell list

#     CellList(; particles, L_x, L_y, cutoff)

# Creates a cell list for `particles` in region/simulation with dimensions
# `L_x`, `L_y`.  `cutoff` gives the approximate size of the cells, the actual size
# of which is chosen so that an integer number of cells fit in each dimension.
# """
# struct VerletList{A <: AbstractParticle}
#     start_pid::Array{Int64, 2}
#     next_pid::Array{Int64, 1}

#     particles::Array{A, 1}

#     num_cells_x::Int64
#     num_cells_y::Int64
#     cell_spacing_x::Float64
#     cell_spacing_y::Float64

#     function CellList(; particles::Array{A, 1}, L_x::Float64, L_y::Float64, cutoff::Float64)
#         num_cells_x = trunc(Int64, L_x / cutoff)
#         num_cells_y = trunc(Int64, L_y / cutoff)
#         cell_spacing_x = L_x / num_cells_x
#         cell_spacing_y = L_y / num_cells_y

#         start_pid = -ones(Int64, num_cells_x, num_cells_y)
#         next_pid = -ones(Int64, length(particles))

#         for (n, particle) in enumerate(particles)
#             i = trunc(Int64, particle.x / cell_spacing_x) + 1
#             j = trunc(Int64, particle.y / cell_spacing_y) + 1

#             if start_pid[i, j] > 0
#                 next_pid[n] = start_pid[i, j]
#             end
#             start_pid[i, j] = n
#         end

#         new(start_pid, next_pid, particles, num_cells_x, num_cells_y, cell_spacing_x, cell_spacing_y)
#     end
# end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interactions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractInteraction
export LennardJones

abstract type AbstractInteraction end

"""
Stores properties of a Lennard-Jones interaction

    LennardJones(; particles, cell_list, ϵ, σ, cutoff, multithreaded, use_newton_3rd)
...
# Arguments
- `particles::Array{Particle, 1}`: particles to compute Lennard-Jones for
- `cell_list::CellList`: cell list for neighbors
- `ϵ::Float64`: energy scale of Lennard-Jones potential
- `σ::Float64 = -1.0`: length scale of Lennard-Jones potential.  If `σ < 0`, use diameters of interacting particles.
- `cutoff::Float64 = -1.0`: distance to truncate interaction.  If `cutoff < 0`, use `cutoff = 2^(1.0 / 6.0) * σ`, which corresponds to WCA interaction.
- `multithreaded::Bool = false`: if `true`, split particles between threads.
- `use_newton_3rd::Bool = false`: if `true`, use Newton's 3rd law to also compute force on neighbors.  Be CAREFUL with overcounting and race conditions with threads.
...
"""
struct LennardJones{T <: AbstractParticle} <: AbstractInteraction
    ϵ::Float64
    σ::Float64
    cutoff::Float64

    particles::Array{T, 1}
    cell_list::CellList

    multithreaded::Bool
    use_newton_3rd::Bool

    function LennardJones(; particles::Array{<: AbstractParticle, 1}, cell_list::CellList,
                            ϵ::Float64, σ::Float64 = -1.0, cutoff::Float64 = -1.0,
                            multithreaded::Bool = false, use_newton_3rd::Bool = false)
        new{<: AbstractParticle}(ϵ, σ, cutoff, particles, cell_list, multithreaded, use_newton_3rd)
    end
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Integrators
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export AbstractIntegrator
export Brownian

abstract type AbstractIntegrator end

"""
Stores properties of a Brownian integrator

    Brownian(; particles, dt, rotations, multithreaded)

Initialize a Brownian integrator for `particles` with timestep `dt`.  If
`rotations == true`, integrate the orientational degree of freedom.  If
`multithreaded == true`, split particles between threads.
"""
struct Brownian{A <: AbstractParticle} <: AbstractIntegrator
    dt::Float64
    particles::Array{AbstractParticle, 1}

    rotations::Bool

    multithreaded::Bool

    function Brownian(; particles::Array{<: AbstractParticle, 1}, dt::Float64, rotations::Bool = false, multithreaded::Bool = false)
        new{<: AbstractParticle}(dt, particles, rotations, multithreaded)
    end
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export Simulation

"""
Stores all information needed to run simulation

    Simulation(...)
...
# Parameters
- `descriptor::String`: a description of the simulation for easy referencing
- `L_x::Float64`: width of simulation region
- `L_y::Float64`: height of simulation region
- `L_z::Float64`: depth of simulation region
- `periodic_in_x::Bool = true`: if `true`, apply periodicity along x
- `periodic_in_y::Bool = true`: if `true`, apply periodicity along y
- `periodic_in_z::Bool = true`: if `true`, apply periodicity along z
- `particles::Array{Particle, 1}`: all particles in simulation
- `cell_lists::Array{CellList, 1}`: all cell lists used
- `verlet_lists::Array{CellList, 1}`: all cell lists used
- `interactions::Array{AbstractInteraction, 1}`: all interactions used
- `dt::Float64`: timestep
- `integrators::Array{AbstractIntegrator}`: all integrators used
- `num_steps::Int64`: total steps to run
- `save_interval::Int64`: store data periodically
...
"""
mutable struct Simulation{A <: AbstractParticle, B <: AbstractInteraction, C <: AbstractIntegrator}
    descriptor::String

    L_d::Vector{Float64}
    periodic_in_d::Vector{Bool}

    particles::Array{A, 1}
    cell_lists::Array{CellList, 1}
    interactions::Array{B, 1}

    dt::Float64
    integrators::Array{C, 1}
    num_steps::Int64

    save_interval::Int64
end

# Outer constructor with default values
function Simulation(; descriptor::String = "No description given...", 
                     L_d::Vector{Float64} = [1.0, 1.0, 0.0], 
                     periodic_in_d::Vector{Bool} = [true, true, false],
                     particles::Array{A, 1} = AbstractParticle[], 
                     cell_lists::Array{CellList, 1} = CellList[], 
                     interactions::Array{B, 1} = AbstractInteraction[],
                     dt::Float64 = 0.0, 
                     integrators::Array{C, 1} = AbstractIntegrator[], 
                     num_steps::Int64 = 0, 
                     save_interval::Int64 = 0) where {A <: AbstractParticle, B <: AbstractInteraction, C <: AbstractIntegrator}

    return Simulation{A, B, C}(descriptor, L_d, periodic_in_d, particles, cell_lists, interactions, dt, integrators, num_steps, save_interval)
end



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Useful functions and macros
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export pbc

"""
    pbc(position; period)

Returns a new position after applying periodic boundary conditions.  The
periodicity is given by `L`. 
"""
@inline Base.@pure function pbc(x::Float64, L::Float64)
    return x - L * floor(x/L + 0.5)
end