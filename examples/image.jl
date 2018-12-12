#!/usr/bin/env julia

import Granular
import FileIO
import Colors
using Random

const verbose = true

const img_file = "aadcroft.png"

img = FileIO.load(img_file)

# resize the image if it is too large, preceed with lopass to avoid antialias
max_pixels = 100^2
if size(img, 1)*size(img, 2) > max_pixels
    cp(img_file, "backup-" * img_file, force=true)
    run(`convert $(img_file) -resize "$(max_pixels)@>" $(img_file)`)
    img = FileIO.load(img_file)
end

const img_bw = Colors.Gray.(img)

const forcing = "gyres"
#const forcing = "down"
#const forcing = "convergent"
#const forcing = "sandpile"

const dx = 1.
const dy = dx

const nx = size(img_bw, 2) + 1
const ny = size(img_bw, 1) + 1

Lx = nx*dx
if forcing == "sandpile"
    Lx *= 1.5
end
const Ly = ny*dy

const youngs_modulus = 2e7
const tensile_strength = 0e3
const h = 0.5

sim = Granular.createSimulation(id="image")

@info "nx = $nx, ny = $ny"

for iy=1:size(img_bw, 1)
    for ix=1:size(img_bw, 2)

        x = ix*dx - dx
        if forcing == "sandpile"
            x += Lx/6.0
        end
        y = Ly - (iy*dy - dy)
        r = 0.5*dx*((1.0 - Float64(img_bw[iy, ix])))

        if r > 0.1*dx
            Granular.addGrainCylindrical!(sim, [x + dx, y - dy], r, h,
                                          tensile_strength=tensile_strength,
                                          youngs_modulus=youngs_modulus,
                                          verbose=verbose)
        end
    end
end

# set ocean forcing
sim.ocean = Granular.createRegularOceanGrid([nx, ny, 1], [Lx, Ly, 1.0],
                                          name="image_ocean")

if forcing == "gyres"
    epsilon = 0.25  # amplitude of periodic oscillations
    t = 0.0
    a = epsilon*sin(2.0*pi*t)
    b = 1.0 - 2.0*epsilon*sin(2.0*pi*t)
    for i=1:size(sim.ocean.u, 1)
        for j=1:size(sim.ocean.u, 2)

            x = sim.ocean.xq[i, j]/(Lx*0.5)  # x in [0;2]
            y = sim.ocean.yq[i, j]/Ly       # y in [0;1]

            f = a*x^2.0 + b*x
            df_dx = 2.0*a*x + b

            sim.ocean.u[i, j, 1, 1] = -pi/10.0*sin(pi*f)*cos(pi*y) * 4e1
            sim.ocean.v[i, j, 1, 1] = pi/10.0*cos(pi*f)*sin(pi*y)*df_dx * 4e1
        end
    end

elseif forcing == "down" || forcing == "sandpile"
    Random.seed!(1)
    sim.ocean.u[:, :, 1, 1] = (rand(nx+1, ny+1) - 0.5)*0.1
    sim.ocean.v[:, :, 1, 1] = -Ly/5.0

elseif forcing == "convergent"
    Random.seed!(1)
    sim.ocean.u[:, :, 1, 1] = (rand(nx+1, ny+1) - 0.5)*0.1
    for j=1:size(sim.ocean.u, 2)
        sim.ocean.v[:, j, 1, 1] = -(j/ny - 0.5)*10.0
    end

else
    error("Forcing not understood")
end

# Initialize confining walls, which are ice floes that are fixed in space
r = dx/4.

## N-S wall segments
for y in range(r, stop=Ly-r, length=Int(round((Ly - 2.0*r)/(r*2))))
    Granular.addGrainCylindrical!(sim, [r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    Granular.addGrainCylindrical!(sim, [Lx-r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end

## E-W wall segments
for x in range(3.0*r, stop=Lx-3.0*r, length=Int(round((Lx - 6.0*r)/(r*2))))
    Granular.addGrainCylindrical!(sim, [x, r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    Granular.addGrainCylindrical!(sim, [x, Ly-r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end


# Finalize setup and start simulation
Granular.setTimeStep!(sim, verbose=verbose)

Granular.setTotalTime!(sim, 5.0)
Granular.setOutputFileInterval!(sim, .1)

Granular.removeSimulationFiles(sim)
Granular.run!(sim, verbose=verbose)
