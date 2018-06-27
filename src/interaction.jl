## Interaction functions
import Compat
using Compat.LinearAlgebra

export interact!
"""
    interact!(simulation::Simulation)

Resolve mechanical interaction between all particle pairs.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
"""
function interact!(simulation::Simulation)
    for i=1:length(simulation.grains)
        for ic=1:simulation.Nc_max

            j = simulation.grains[i].contacts[ic]

            if i > j  # skip i > j and j == 0
                continue
            end

            interactGrains!(simulation, i, j, ic)
        end
    end

    for i=1:length(simulation.grains)
        @inbounds simulation.grains[i].granular_stress = 
            simulation.grains[i].force/
            simulation.grains[i].horizontal_surface_area
    end
    nothing
end


"""
    interactWalls!(sim)

Find and resolve interactions between the dynamic walls (`simulation.walls`) and
the grains.  The contact model uses linear elasticity, with stiffness based on
the grain Young's modulus `grian.E` or elastic stiffness `grain.k_n`.  The
interaction is frictionless in the tangential direction.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains and
    dynamic walls.
"""
function interactWalls!(sim::Simulation)

    orientation::Float64 = 0.0
    δ_n::Float64 = 0.0
    k_n::Float64 = 0.0
    γ_n::Float64 = 0.0
    vel_n::Float64 = 0.0
    force_n::Float64 = 0.0

    for iw=1:length(sim.walls)
        for i=1:length(sim.grains)

            orientation = sign(dot(sim.walls[iw].normal,
                                   sim.grains[i].lin_pos) - sim.walls[iw].pos)

            # get overlap distance by projecting grain position onto wall-normal
            # vector. Overlap when δ_n < 0.0
            δ_n = norm(dot(sim.walls[iw].normal[1:2],
                           sim.grains[i].lin_pos[1:2]) -
                       sim.walls[iw].pos) - sim.grains[i].contact_radius

            vel_n = dot(sim.walls[iw].normal[1:2], sim.grains[i].lin_vel[1:2])

            if δ_n < 0.
                if sim.grains[i].youngs_modulus > 0.
                    k_n = sim.grains[i].youngs_modulus * sim.grains[i].thickness
                else
                    k_n = sim.grains[i].contact_stiffness_normal
                end

                γ_n = sim.walls[iw].contact_viscosity_normal

                if k_n ≈ 0. && γ_n ≈ 0.  # No interaction
                    force_n = 0.

                elseif k_n > 0. && γ_n ≈ 0.  # Elastic (Hookean)
                    force_n = k_n * δ_n

                elseif k_n > 0. && γ_n > 0.  # Elastic-viscous (Kelvin-Voigt)
                    force_n = k_n * δ_n + γ_n * vel_n

                    # Make sure that the visous component doesn't dominate over
                    # elasticity
                    if force_n > 0.
                        force_n = 0.
                    end

                else
                    error("unknown contact_normal_rheology " *
                          "(k_n = $k_n, γ_n = $γ_n")
                end

                sim.grains[i].force += -force_n .* sim.walls[iw].normal .*
                                               orientation
                sim.walls[iw].force += force_n * orientation
            end
        end
    end
end

export interactGrains!
"""
    interactGrains!(simulation::Simulation, i::Int, j::Int, ic::Int)

Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the grain 
`pressure` field of mean compressive stress on the grain sides.

The function uses the macroscopic contact-stiffness parameterization based on 
Young's modulus and Poisson's ratio if Young's modulus is a positive value.
"""
function interactGrains!(simulation::Simulation, i::Int, j::Int, ic::Int)

    if !simulation.grains[i].enabled || !simulation.grains[j].enabled
        return
    end

    force_n = 0.  # Contact-normal force
    force_t = 0.  # Contact-parallel (tangential) force

    # Inter-position vector, use stored value which is corrected for boundary
    # periodicity
    p = simulation.grains[i].position_vector[ic][1:2]
    dist = norm(p)

    r_i = simulation.grains[i].contact_radius
    r_j = simulation.grains[j].contact_radius

    # Floe distance; <0: compression, >0: tension
    δ_n = dist - (r_i + r_j)

    # Local axes
    n = p/dist
    t = [-n[2], n[1]]

    # Contact kinematics (2d)
    vel_lin = simulation.grains[i].lin_vel[1:2] -
        simulation.grains[j].lin_vel[1:2]
    vel_n = dot(vel_lin, n)

    if !simulation.grains[i].rotating && !simulation.grains[j].rotating 
        rotation = false
    else
        rotation = true
    end
    
    if rotation
        vel_t = dot(vel_lin, t) -
            harmonicMean(r_i, r_j)*(simulation.grains[i].ang_vel[3] +
                                    simulation.grains[j].ang_vel[3])

        # Correct old tangential displacement for contact rotation and add new
        δ_t0 = simulation.grains[i].contact_parallel_displacement[ic][1:2]
        δ_t = dot(t, δ_t0 - (n*dot(n, δ_t0))) + vel_t*simulation.time_step

        # Determine the contact rotation (2d)
        θ_t = simulation.grains[i].contact_rotation[ic][3] +
            (simulation.grains[j].ang_vel[3] - 
             simulation.grains[i].ang_vel[3]) * simulation.time_step
    end

    # Effective radius
    R_ij = harmonicMean(r_i, r_j)

    # Determine which ice floe is the thinnest
    h_min, idx_thinnest = findmin([simulation.grains[i].thickness,
                                   simulation.grains[j].thickness])
    idx_thinnest == 1 ? idx_thinnest = i : idx_thinnest = j
    
    # Contact length along z
    Lz_ij = h_min

    # Contact area
    A_ij = R_ij*Lz_ij

    # Contact mechanical parameters
    if simulation.grains[i].youngs_modulus > 0. &&
        simulation.grains[j].youngs_modulus > 0.

        E = harmonicMean(simulation.grains[i].youngs_modulus,
                         simulation.grains[j].youngs_modulus)
        ν = harmonicMean(simulation.grains[i].poissons_ratio,
                         simulation.grains[j].poissons_ratio)

        # Effective normal and tangential stiffness
        k_n = E * A_ij/R_ij
        #k_t = k_n*ν   # Kneib et al 2016
        if rotation
            k_t = k_n * 2. * (1. - ν^2.) / ((2. - ν) * (1. + ν))  # Obermayr 2011
        end

    else  # Micromechanical parameterization: k_n and k_t set explicitly
        k_n = harmonicMean(simulation.grains[i].contact_stiffness_normal,
                           simulation.grains[j].contact_stiffness_normal)

        if rotation
            k_t = harmonicMean(simulation.grains[i].contact_stiffness_tangential,
                               simulation.grains[j].contact_stiffness_tangential)
        end
    end

    γ_n = harmonicMean(simulation.grains[i].contact_viscosity_normal,
                       simulation.grains[j].contact_viscosity_normal)

    if rotation
        γ_t = harmonicMean(simulation.grains[i].contact_viscosity_tangential,
                           simulation.grains[j].contact_viscosity_tangential)

        μ_d_minimum = min(simulation.grains[i].contact_dynamic_friction,
                          simulation.grains[j].contact_dynamic_friction)
    end

    # Determine contact forces
    if k_n ≈ 0. && γ_n ≈ 0.  # No interaction
        force_n = 0.

    elseif k_n > 0. && γ_n ≈ 0.  # Elastic (Hookean)
        force_n = -k_n*δ_n

    elseif k_n > 0. && γ_n > 0.  # Elastic-viscous (Kelvin-Voigt)
        force_n = -k_n*δ_n - γ_n*vel_n
        if force_n < 0.
            force_n = 0.
        end

    else
        error("unknown contact_normal_rheology (k_n = $k_n, γ_n = $γ_n,
              E = $E, A_ij = $A_ij, R_ij = $R_ij)
              ")
    end

    # Determine which grain is the weakest
    compressive_strength = simulation.grains[i].compressive_strength * 
                           sqrt(Lz_ij)
    compressive_strength_j = simulation.grains[j].compressive_strength * 
                             sqrt(Lz_ij)
    idx_weakest = i
    if compressive_strength_j < compressive_strength
        compressive_strength = compressive_strength_j
        idx_weakest = j
    end

    # Determine the compressive strength in Pa by the contact thickness and the
    # minimum compressive strength
    compressive_strength = min(simulation.grains[i].compressive_strength,
                               simulation.grains[j].compressive_strength)/
                           sqrt(Lz_ij)

    # Limit compressive stress if the prefactor is set to a positive value
    if δ_n <= 0.0 && compressive_strength > 0. &&
        norm(force_n) >= compressive_strength*A_ij

        # Register that compressive failure has occurred for this contact
        simulation.grains[i].compressive_failure[ic] = 1
    end

    # Grain-pair moment of inertia [m^4]
    if rotation
        I_ij = 2.0/3.0*R_ij^3*min(simulation.grains[i].thickness,
                                  simulation.grains[j].thickness)
    end

    # Contact tensile strength increases linearly with contact age until
    # tensile stress exceeds tensile strength.
    tensile_strength = min(simulation.grains[i].contact_age[ic]*
                           simulation.grains[i].strength_heal_rate,
                           simulation.grains[i].tensile_strength)
    if rotation
        shear_strength = min(simulation.grains[i].contact_age[ic]*
                             simulation.grains[i].strength_heal_rate,
                             simulation.grains[i].shear_strength)
        M_t = 0.0
        if tensile_strength > 0.0
            # Determine bending momentum on contact [N*m],
            # (converting k_n to E to bar(k_n))
            M_t = (k_n*R_ij/(A_ij*(simulation.grains[i].contact_radius +
                               simulation.grains[j].contact_radius)))*I_ij*θ_t
        end
    end

    # Reset contact age (breaking bond) if bond strength is exceeded
    if rotation
        if δ_n >= 0.0 && norm(force_n)/A_ij + norm(M_t)*R_ij/I_ij > tensile_strength
            force_n = 0.
            force_t = 0.
            simulation.grains[i].contacts[ic] = 0  # remove contact
            simulation.grains[i].n_contacts -= 1
            simulation.grains[j].n_contacts -= 1
        end
    else
        if δ_n >= 0.0 && norm(force_n)/A_ij > tensile_strength
            force_n = 0.
            force_t = 0.
            simulation.grains[i].contacts[ic] = 0  # remove contact
            simulation.grains[i].n_contacts -= 1
            simulation.grains[j].n_contacts -= 1
        end
    end

    # Limit compressive stress if the prefactor is set to a positive value
    if δ_n <= 0.0 && compressive_strength > 0. &&
        norm(force_n) >= compressive_strength

        # Determine the overlap distance where yeild stress is reached
        δ_n_yield = -compressive_strength*A_ij/k_n

        # Determine the excess overlap distance relative to yield
        Δr = norm(δ_n) - norm(δ_n_yield)
        
        # Adjust radius and thickness of the weakest grain
        simulation.grains[idx_weakest].contact_radius -= Δr
        simulation.grains[idx_weakest].areal_radius -= Δr
        simulation.grains[idx_weakest].thickness += 1.0/(π*Δr)
    end

    if rotation
        if k_t ≈ 0. && γ_t ≈ 0.
            # do nothing

        elseif k_t ≈ 0. && γ_t > 0.
            force_t = norm(γ_t * vel_t)

            # Coulomb slip
            if force_t > μ_d_minimum*norm(force_n)
                force_t = μ_d_minimum*norm(force_n)

                # Nguyen et al 2009 "Thermomechanical modelling of friction effects
                # in granular flows using the discrete element method"
                E_shear = norm(force_t)*norm(vel_t)*simulation.time_step

                # Assume equal thermal properties
                simulation.grains[i].thermal_energy += 0.5*E_shear
                simulation.grains[j].thermal_energy += 0.5*E_shear
            end
            if vel_t > 0.
                force_t = -force_t
            end

        elseif k_t > 0.

            force_t = -k_t*δ_t - γ_t*vel_t

            # Coulomb slip
            if norm(force_t) > μ_d_minimum*norm(force_n)
                force_t = μ_d_minimum*norm(force_n)*force_t/norm(force_t)
                δ_t = (-force_t - γ_t*vel_t)/k_t

                # Nguyen et al 2009 "Thermomechanical modelling of friction effects
                # in granular flows using the discrete element method"
                E_shear = norm(force_t)*norm(vel_t)*simulation.time_step

                # Assume equal thermal properties
                simulation.grains[i].thermal_energy += 0.5*E_shear
                simulation.grains[j].thermal_energy += 0.5*E_shear
            end

        else
            error("unknown contact_tangential_rheology (k_t = $k_t, γ_t = $γ_t")
        end
    end

    # Break bond under extension through bending failure
    if δ_n < 0.0 && tensile_strength > 0.0 && rotation && 
        shear_strength > 0.0 && norm(force_t)/A_ij > shear_strength

        force_n = 0.
        force_t = 0.
        simulation.grains[i].contacts[ic] = 0  # remove contact
        simulation.grains[i].n_contacts -= 1
        simulation.grains[j].n_contacts -= 1
    end

    simulation.grains[i].contact_age[ic] += simulation.time_step

    if rotation
        simulation.grains[i].contact_parallel_displacement[ic] = vecTo3d(δ_t.*t)
        simulation.grains[i].contact_rotation[ic] = [0., 0., θ_t]

        simulation.grains[i].force += vecTo3d(force_n.*n + force_t.*t);
        simulation.grains[j].force -= vecTo3d(force_n.*n + force_t.*t);

        simulation.grains[i].torque[3] += -force_t*R_ij + M_t
        simulation.grains[j].torque[3] += -force_t*R_ij - M_t
    else
        simulation.grains[i].force += force_n*n;
        simulation.grains[j].force -= force_n*n;
    end

    simulation.grains[i].pressure += 
        force_n/simulation.grains[i].side_surface_area;
    simulation.grains[j].pressure += 
        force_n/simulation.grains[j].side_surface_area;
    nothing
end
