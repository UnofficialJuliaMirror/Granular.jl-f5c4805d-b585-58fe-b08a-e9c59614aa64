## Functions for creating grain packings

export regularPacking!
"""

    regularPacking!(simulation, n, r_min, r_max[, padding_factor,
                    size_distribution, size_distribution_parameter, seed])

Create a grid-based regular packing with grain numbers along each axis specified
by the `n` vector.

# Arguments
* `simulation::Simulation`: simulation object where the grains are inserted,
    preferably not containing prior grains.
* `n::Vector{Integer}`: 2-element vector determining number of grains along the
    `x` and `y` axes.
* `r_min::Real`: minimum desired grain radius.
* `r_max::Real`: maximum desired grain radius.
* `padding_factor::Real`: percentage-wise padding around each grain to allow for
    random perturbations to grain position.
* `size_distribution::String`: grain-size distribution to sample. Valid values
    are "powerlaw" and "uniform".
* `size_distribution_parameter::Real`: parameter to pass to the grain-size
    distribution generating function.
* `seed::Integer`: seed value to the pseudo-random number generator.
"""
function regularPacking!(simulation::Simulation,
                         n::Vector{Int},
                         r_min::Real,
                         r_max::Real;
                         padding_factor::Real=.1,
                         size_distribution::String="powerlaw",
                         size_distribution_parameter::Real=-1.8,
                         seed::Integer=1)

    r_rand = 0.
    pos = zeros(2)
    h = .5   # disc tickness
    dx = r_max * 2. * (1. + padding_factor)  # cell size
    dx_padding = r_max * 2. * padding_factor
    srand(seed)

    for iy in 1:n[2]
        for ix in 1:n[1]

            if size_distribution == "powerlaw"
                r_rand = Granular.randpower(1, size_distribution_parameter,
                                            r_min, r_max)
            elseif size_distribution == "uniform"
                r_rand = rand()*(r_max - r_min) + r_min
            end

            # Determine position from grid index and sample randomly from within
            # padding
            pos = [ix*dx - .5*dx, iy*dx - .5*dx] .+
                rand(2) .* dx_padding .- .5*dx_padding

            addGrainCylindrical!(simulation, pos, r_rand, h, verbose=false)
        end
    end

end

"""
Return random point in spherical annulus between `(r_i + r_j)` and `2.*(r_i +
r_j)` around `x_i`.  Note: there is slightly higher point density towards (r_i +
r_j).
"""
function generateNeighboringPoint(x_i::Vector, r_i::Real,
                                  r_max::Real, r_min::Real;
                                  padding::Real=0.)

    if padding > 0.
        r_j = r_min + (rand()*0.5 + 0.5)*(r_max - r_min)
    else
        r_j = r_min + rand()*(r_max - r_min)  # between r_min and r_max
    end
    #r_j = r_min + rand()*(r_i - r_min)  # between r_min and r_i
    #R = rand() * (r_i + r_j) * max_padding_factor + 2. * (r_i + r_j)
    R = r_i + r_j + padding
    T = rand() * 2. * pi
    return [x_i[1] + R * sin(T), x_i[2] + R * cos(T)], r_j
end

function generateRandomDirection()
    return rand() * 2. * pi
end

function getPositionDistancedFromPoint(T::Real, x_i::Vector, dist::Real)
    return [x_i[1] + dist * sin(T), x_i[2] + dist * cos(T)]
end

export irregularPacking!
"""
    irregularPacking!(simulation[, radius_max, radius_min, sample_limit,
                      thickness, seed, plot_during_packing, verbose)

Generate a dense disc packing in 2D using Poisson disc sampling with O(N)
complexity, as described by [Robert Bridson (2007) "Fast Poisson disk sampling
in arbitrary dimensions"](https://doi.org/10.1145/1278780.1278807). The
`simulation` can be empty or already contain grains. However, an
`simulation.ocean` or `simulation.atmosphere` grid is required.

# Arguments
* `simulation::Simulation`: simulation object where grains are inserted.
* `radius_max::Real`: largest grain radius to use.
* `radius_min::Real`: smallest grain radius to use.
* `sample_limit::Integer=30`: number of points to sample around each grain
    before giving up.
* `seed::Integer`: seed value to the pseudo-random number generator.
* `plot_during_packing::Bool=false`: produce successive plots as the packing is
    generated. Requires gnuplot (default).
* `verbose::Bool=true`: show diagnostic information to stdout.
"""
function irregularPacking!(simulation::Simulation;
                           radius_max::Real=.1,
                           radius_min::Real=.1,
                           sample_limit::Integer=30,
                           padding_factor::Real=2.,
                           binary_radius_search::Bool=false,
                           thickness::Real=1.,
                           seed::Integer=1,
                           plot_during_packing::Bool=false,
                           verbose::Bool=true)
    srand(seed)

    active_list = Int[]  # list of points to originate search from

    # Step 0: Use existing `grid` (ocean or atmosphere) for contact search
    if typeof(simulation.ocean.input_file) != Bool
        grid = simulation.ocean
    elseif typeof(simulation.atmosphere.input_file) != Bool
        grid = simulation.atmosphere
    else
        error("irregularPacking requires an ocean or atmosphere grid")
    end
    # save grid boundaries
    sw, se, ne, nw = getGridCornerCoordinates(grid.xq, grid.yq)
    width_x = se[1] - sw[1]  # assume regular grid
    width_y = nw[2] - sw[2]  # assume regular grid

    # Step 1: If grid is empty: select random initial sample and save its index
    # to the background grid. Otherwise mark all existing grains as active
    np_init = length(simulation.grains)
    if isempty(simulation.grains)
        r = rand()*(radius_max - radius_min) + radius_min
        x0 = rand(2).*[width_x, width_y] + sw
        addGrainCylindrical!(simulation, x0, r, thickness)
        sortGrainsInGrid!(simulation, grid)
        push!(active_list, 1)
    else
        for i=1:length(simulation.grains)
            push!(active_list, i)
        end
    end

    # Step 2: While the active list is not empty, choose a random index `i` from
    # it.  Generate up to `sample_limit` points chosen uniformly from the
    # distance `(r_i+r_j)` around `x_i`.
    # For each point in turn, check if it is within distance r of existing
    # samples (using the background grid to only test nearby samples). If a
    # point is adequately far from existing samples, emit it as the next sample
    # and add it to the active list. If after `sample_limit` attempts no such
    # point is found, instead remove `i` from the active list.
    i = 0; j = 0;
    x_active = zeros(2); x_candidate = zeros(2);
    r_active = 0.; r_candidate = 0.; T = 0.
    n = 0
    neighbor_found = false

    while !isempty(active_list)

        # Draw a random grain from the list of active grains
        i = active_list[rand(1:length(active_list))]

        x_active = simulation.grains[i].lin_pos
        r_active = simulation.grains[i].contact_radius

        # Did the algoritm find a neighbor to the current active grain `i`?
        neighbor_found = false

        for j=1:sample_limit

            # Generate a candidate point
            if binary_radius_search
                # Generate a point positioned at r_active + radius_max from the
                # position x_active.
                T = generateRandomDirection()
                r_candidate = radius_max
                x_candidate = getPositionDistancedFromPoint(T, x_active,
                                                            r_active + r_candidate)
            else
                if j <= 2
                    x_candidate, r_candidate = generateNeighboringPoint(
                                                   x_active,
                                                   r_active,
                                                   radius_max,
                                                   radius_min,
                                                   padding=padding_factor*radius_max)
                else
                    x_candidate, r_candidate = generateNeighboringPoint(
                                                   x_active,
                                                   r_active,
                                                   radius_max,
                                                   radius_min)
                end
            end

            # Make sure that the point is within the current grid cell
            if !(isPointInGrid(grid, x_candidate))
                continue  # skip this candidate
            end

            # If the binary_radius_search is selected, try to adjust the radius
            # to a value as large as possible
            if binary_radius_search

                radius_not_found = true

                # first test the maximum radius. If unsuccessful, iteratively
                # find the optimal radius using binary searches
                if !checkForContacts(simulation, grid, x_candidate, r_candidate)

                    # 1. Set L to min and R to max of sampling range
                    r_L = radius_min
                    r_R = radius_max

                    # size of radius sampling step
                    dr = (r_R - r_L)/25.

                    while radius_not_found

                        # 2. If L > R, the search terminates as unsuccessful
                        if r_L > r_R
                            radius_not_found = false
                            break
                        end

                        # 3. Set r to the middle of the current range
                        r_candidate = (r_L + r_R)/2.0
                        x_candidate = getPositionDistancedFromPoint(T, x_active,
                                        r_active + r_candidate)

                        # 4. If r < target, set L to r+dr and go to step 2
                        if checkForContacts(simulation, grid, x_candidate,
                                            r_candidate)
                            r_L = r_candidate + dr

                        else # 5. If r > target, set R to r-dr and go to step 2
                            r_R = r_candidate - dr
                        end
                    end
                end

            end

            # if the grain candidate doesn't overlap with any other grains,
            # add it and mark it as active
            if checkForContacts(simulation, grid, x_candidate, r_candidate)
                #println("Added grain from parent $i")
                addGrainCylindrical!(simulation, x_candidate, r_candidate,
                                     thickness, verbose=false)
                sortGrainsInGrid!(simulation, grid)
                push!(active_list, length(simulation.grains))
                break
            end

            if j == sample_limit
                # If no neighbors were found, delete the grain `i` from the list
                # of active grains
                filter!(f->f≠i, active_list)
            end
        end
        if verbose
            println("Active points: $(length(active_list))")
            #println(active_list)
        end

        if plot_during_packing
            n += 1
            filepostfix = @sprintf("packing.%05d.png", n)
            plotGrains(simulation, filetype=filepostfix, show_figure=false)
        end

    end  # end while !isempty(active_list)

    if verbose
        info("Generated $(length(simulation.grains) - np_init) points")
    end
end
