#!/usr/bin/env julia
using Compat.Test
import Granular

verbose = true
debug = true
if debug
    import PyPlot
end

function plot_interaction(sim::Granular.Simulation, output::String)
    time = Float64[]
    force_n_1 = Float64[]
    force_n_2 = Float64[]
    force_t_1 = Float64[]
    force_t_2 = Float64[]
    torque_1 = Float64[]
    torque_2 = Float64[]
    compressive_failure = Bool[]
    while sim.time < sim.time_total
        Granular.run!(sim, verbose=verbose, single_step=true)
        append!(time, sim.time)
        append!(compressive_failure, sim.grains[1].compressive_failure[1])
        append!(force_n_1, sim.grains[1].force[1])
        append!(force_n_2, sim.grains[2].force[1])
        append!(force_t_1, sim.grains[1].force[2])
        append!(force_t_2, sim.grains[2].force[2])
        append!(torque_1, sim.grains[1].torque[3])
        append!(torque_2, sim.grains[2].torque[3])
    end
    PyPlot.subplot(3,1,1)
    PyPlot.plot(time, force_n_1, "-b", label="1")
    PyPlot.plot(time, force_n_2, "--y", label="2")
    PyPlot.legend(loc="upper right")
    PyPlot.ylabel("Normal force [N]")
    PyPlot.subplot(3,1,2)
    PyPlot.plot(time, force_t_1, "-b", label="1")
    PyPlot.plot(time, force_t_2, "--y", label="2")
    PyPlot.legend(loc="upper right")
    PyPlot.ylabel("Tangential force [N]")
    PyPlot.subplot(3,1,3)
    PyPlot.plot(time, torque_1, "-b", label="1")
    PyPlot.plot(time, torque_2, "--y", label="2")
    PyPlot.legend(loc="upper right")
    PyPlot.ylabel("Torque [Nm]")
    PyPlot.xlabel("Time [s]")
    PyPlot.tight_layout()
    PyPlot.savefig(output)
end

Compat.@info "Testing compressive failure"
sim = Granular.createSimulation()
Granular.addGrainCylindrical!(sim, [0.,0.], 1., 0.5,
                              compressive_strength=1285e3,
                              lin_vel=[1., 0.], fixed=true)
Granular.addGrainCylindrical!(sim, [2.,0.], 1., 0.5,
                              compressive_strength=1285e3,
                              fixed=true)
@test count(x->x==true, sim.grains[1].compressive_failure) == 0
Granular.setTimeStep!(sim, verbose=verbose)
Granular.setTotalTime!(sim, 1.0)

if debug
    Granular.setOutputFileInterval!(sim, 0.01)
    plot_interaction(sim, "compressive_failure.pdf")
    #Granular.render(sim)
else
    Granular.run!(sim, verbose=verbose)
end

@test sim.grains[1].compressive_failure[1] == true
@test count(x->x==true, sim.grains[1].compressive_failure) == 1
@test sim.grains[1].force[1] < 0.0
@test sim.grains[1].force[2] ≈ 0.0
@test sim.grains[2].force[1] > 0.0
@test sim.grains[2].force[2] ≈ 0.0
