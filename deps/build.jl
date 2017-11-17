
#!/usr/bin/env julia
using BinDeps
@BinDeps.setup

import Compat
import Compat.Sys

gnuplot = library_dependency("gnuplot")
imagemagick = library_dependency("imagemagick", aliases = ["ImageMagick"])

# Debian derivatives: https://www.debian.org/distrib/packages#search_packages
provides(AptGet, "gnuplot", gnuplot, os = :Linux)
provides(AptGet, "imagemagick", imagemagick, os = :Linux)

# RHEL derivatives: http://rpm.pbone.net/index.php3/stat/2/simple/2
provides(Yum, "gnuplot", gnuplot, os = :Linux)
provides(Yum, Dict("ImageMagick" => imagemagick), os = :Linux)

# Arch: https://www.archlinux.org/packages/
provides(Pacman, "gnuplot", gnuplot, os = :Linux)
provides(Pacman, "imagemagick", imagemagick, os = :Linux)

# Mac: http://formulae.brew.sh/
if Compat.Sys.isapple()
    using Homebrew
    provides(Homebrew.HB, "gnuplot", gnuplot, os = :Darwin)
    provides(Homebrew.HB, "imagemagick", imagemagick, os = :Darwin)
end

# Windows: http://software.opensuse.org/search
if Compat.Sys.iswindows()
    using WinRPM
    provides(WinRPM.RPM, "gnuplot", gnuplot, os = :Windows)
    provides(WinRPM.RPM, "imagemagick", imagemagick, os = :Windows)
end

@BinDeps.install Dict([
                       (:gnuplot => :gnuplot),
                       (:imagemagick => :imagemagick),
                      ])
