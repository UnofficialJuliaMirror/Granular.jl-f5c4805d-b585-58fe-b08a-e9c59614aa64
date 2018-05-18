#!/bin/bash

# Optionally use this script to launch profiling.jl with different Julia 
# optimization levels.  Defaults are --optimize=2, --math-mode=ieee.

set -e

declare -a arr=(\
    "--procs 1 --optimize=1 --math-mode=ieee" \
    "--procs 1 --optimize=3 --math-mode=ieee" \
    "--procs 1 --optimize=1 --math-mode=fast" \
    "--procs 1 --optimize=3 --math-mode=fast")

for flags in "${arr[@]}"; do
    for nthreads in 1 2 4; do
        echo "JULIA_NUM_THREADS=$nthreads julia --color=yes $flags profiling.jl"
        JULIA_NUM_THREADS=$nthreads julia --color=yes $flags profiling.jl
        mv profiling-cpu.pdf "profiling-cpu.${flags}.${nthreads}threads.pdf"
        mv profiling-memory-usage.pdf \
            "profiling-memory-usage.${flags}.${nthreads}threads.pdf"
    done
done
