"""

"""

export initialise_positions_fixed_density

function initialise_positions_fixed_density(num_particles::Int64, ϕ::Float64 = 0.3)
    # Determine the number of balls per row and column
    N = num_particles
    num_per_side = ceil(Int, sqrt(N))
    
    # Calculate the spacing between balls
    side_length = sqrt(N*π / ϕ)/2 
    spacing = side_length / num_per_side
    
    # Initialize arrays for coordinates
    particle_positions = SVector{3, Float64}[]
    
    # Generate coordinates
    for i in 0:num_per_side-1
        for j in 0:num_per_side-1
            x_coords = -side_length/2 + i * spacing + spacing / 2
            y_coords = -side_length/2 + j * spacing + spacing / 2
            push!(particle_positions, SVector(x_coords, y_coords, 0.0))
        end
    end

    return particle_positions, [side_length, side_length, 0.0]
end
