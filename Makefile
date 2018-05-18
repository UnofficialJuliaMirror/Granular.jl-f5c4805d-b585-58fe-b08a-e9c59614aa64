default: test

.PHONY: test
test: \
	test-julia-0.6 \
	test-julia-0.6-2threads \
	test-julia-0.7 \
	test-julia-0.7-2threads \


.PHONY: test-julia-0.6
test-julia-0.6:
	julia --color=yes -e 'Pkg.test("Granular")' \
		&& notify-send Granular.jl tests completed successfully on Julia 0.6 \
		|| notify-send Granular.jl failed on Julia 0.6

.PHONY: test-julia-0.6-2threads
test-julia-0.6-2threads:
	JULIA_NUM_THREADS=2 julia --color=yes --optimize=3 --math-mode=fast \
					  -e 'Pkg.test("Granular")' \
		&& notify-send Granular.jl tests with two threads completed successfully on Julia 0.6 \
		|| notify-send Granular.jl tests with two threads failed on Julia 0.6

.PHONY: test-julia-0.7
test-julia-0.7:
	julia-0.7 --color=yes -e 'import Pkg; Pkg.test("Granular")' \
		&& notify-send Granular.jl tests completed successfully on Julia 0.7 \
		|| notify-send Granular.jl failed on Julia 0.7

.PHONY: test-julia-0.7-2threads
test-julia-0.7-2threads:
	JULIA_NUM_THREADS=2 julia-0.7 --color=yes --optimize=3 --math-mode=fast \
					  -e 'Pkg.test("Granular")' \
		&& notify-send Granular.jl tests with two threads completed successfully on Julia 0.7 \
		|| notify-send Granular.jl tests with two threads failed on Julia 0.7


.PHONY: docs
docs:
	cd docs && julia --color=yes make.jl
	open docs/build/index.html

