#!/usr/bin/env julia
using Compat

if is_apple() 
    using Homebrew
    Homebrew.update()
    Homebrew.add("gnuplot")
    Homebrew.add("imagemagick")
end

if is_windows() 
    using WinRPM
    WinRPM.update()
    WinRPM.install("gnuplot")
    WinRPM.install("imagemagick")
end

