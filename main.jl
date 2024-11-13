using ActiveMatter

# Initialise the particle positions with a volume fraction of 0.3
num_particles = 100
ϕ = 0.3
particle_positions, L_d = initialise_positions_fixed_density(num_particles, ϕ);

# Initialise the Brownian particles
brownian_particles = BrownianParticle[]
for x_d ∈ particle_positions
    push!(brownian_particles, BrownianParticle(x_d = x_d))
end

# Initialise the verlet and cell lists


# sim = Simulation()
