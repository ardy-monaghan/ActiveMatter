function linked_grids(num_cells_d)
    num_cells = num_cells_d[1]
    cells_to_check = [zeros(Int, 9) for _ in 1:num_cells^2]

    for j in 0:num_cells-1, i in 0:num_cells-1
        idx = i + j * num_cells + 1
        cells_to_check[idx][1] = idx
        cells_to_check[idx][2] = mod(i-1, num_cells) + mod(j-1, num_cells) * num_cells + 1
        cells_to_check[idx][3] = mod(i-1, num_cells) + j * num_cells + 1
        cells_to_check[idx][4] = mod(i-1, num_cells) + mod(j+1, num_cells) * num_cells + 1
        cells_to_check[idx][5] = i + mod(j+1, num_cells) * num_cells + 1
        cells_to_check[idx][6] = mod(i+1, num_cells) + mod(j+1, num_cells) * num_cells + 1
        cells_to_check[idx][7] = mod(i+1, num_cells) + j * num_cells + 1
        cells_to_check[idx][8] = mod(i+1, num_cells) + mod(j-1, num_cells) * num_cells + 1
        cells_to_check[idx][9] = i + mod(j-1, num_cells) * num_cells + 1
    end

    return cells_to_check
end

# Function to assign particles to grids based on position
function particle_grid!(p, x, y, grid_particles, current_grid)
    num_particles = p.N
    grid_length = p.grid_length
    num_cells_d = p.num_cells_d
    L_half = p.L / 2

    for pᵢ in 1:num_particles
        grid_x = Int(floor((pbc(x[pᵢ], p.L) + L_half) / grid_length) % num_cells_d) + 1
        grid_y = Int(floor((pbc(y[pᵢ], p.L) + L_half) / grid_length) % num_cells_d) + 1
        new_grid_idx = (grid_x - 1) * num_cells_d + grid_y

        # Get old grid index and update only if changed
        old_grid_idx = current_grid[pᵢ]
        if new_grid_idx != old_grid_idx
            # Move particle to new grid
            push!(grid_particles[new_grid_idx], pᵢ)
            deleteat!(grid_particles[old_grid_idx], findfirst(==(pᵢ), grid_particles[old_grid_idx]))

            # Update particle's current grid
            current_grid[pᵢ] = new_grid_idx
        end
    end
end