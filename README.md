# EAGODomainReduction.jl
Domain Reduction Procedures in Global Optimization

[![Build Status](https://travis-ci.org/MatthewStuber/EAGODomainReduction.jl.svg?branch=master)](https://travis-ci.org/MatthewStuber/EAGODomainReduction.jl)
[![Coverage Status](https://coveralls.io/repos/github/MatthewStuber/EAGODomainReduction.jl/badge.svg?branch=master)](https://coveralls.io/github/MatthewStuber/EAGODomainReduction.jl?branch=master)
[![codecov.io](http://codecov.io/github/MatthewStuber/EAGODomainReduction.jl/coverage.svg?branch=master)](http://codecov.io/github/MatthewStuber/EAGODomainReduction.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/latest)

## Authors

 [Matthew Wilhelm](httppsor.uconn.eduour-team), Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)

## Installation

```julia
julia> Pkg.add("EAGODomainReduction.jl")
```

## Capabilities

**EAGODomainReduction.jl** provides a series of subroutines for tightening the domains of subproblems
solved in global optimization (potentially to in-feasibility). Currently, it supports routines for
nonconvex nonlinear programs:
- **Interval Contractor Propagation:** A forward-backward interval contractor using **IntervalArthimetic.jl**
and **IntervalContractors.jl** for the operator library.
- **Duality-Based Bound Tightening:** Provides algorithms for tightening domains based on duality of solutions
found for subproblems.
- **Standard Range Reduction:** Contracts subproblem domain via linear-relaxations generated using McCormick
relaxations.
- **Implicit Subroutine support:** Supports domain reduction of reduced space lower-bound problems defined through
relaxation of implicit functions by fixed-point methods.

The routine are used extensively in the [`EAGO.jl`](https://github.com/MatthewStuber/EAGO.jl) package solver.
Please see the example files for usage cases.

## Future Work

- Update the interval-constraint propagation algorithm to incorporate propagation heuristics in Vu2008.
- Incorporate control-flow syntax support into constraint propagation algorithm.
- Incorporate improvements to probing and optimality based-bound tightening.
- Add support for Mixed-Integer NLP.

## Related Packages
- [**EAGO.jl**](https://github.com/MatthewStuber/EAGO.jl): A package containing global and robust solvers based mainly on McCormick relaxations.
This package supports a JuMP and MathProgBase interface.
- [**IntervalConstraintProgramming.jl**](https://github.com/JuliaIntervals/IntervalConstraintProgramming.jl): Provides algorithms that furnish bounds
on constraints defined by expressions. The constraint propagation routine in **EAGODomainReduction.jl** can generate tape objects that are
reusable for generically-defined functions. In addition, we use a `Vector{Interval}` storage object that allows for in-place mutation of intervals.
- [**IntervalContractors.jl**](https://github.com/JuliaIntervals/IntervalContractors.jl): Provides a library of reverse interval contractors.

## References
- Benhamou, F., & Older, W.J. (1997). Applying interval arithmetic to real, integer, and boolean constraints. The Journal of Logic Programming, 32, 1–24.
- Caprara, A., & Locatelli, M. (2010). Global optimization problems and domain reduction strategies. Mathematical Programming, 125, 123–137.
- Gleixner, A.M., Berthold, T., Müller, B., & Weltge, S. (2016). Three enhancements for optimization-based bound tightening. ZIB Report, 15–16.
- Ryoo, H.S., & Sahinidis, N.V. (1996). A branch-and-reduce approach to global optimization. Journal of Global Optimization, 8, 107–139.
- Schichl, H., & Neumaier, A. (2005). Interval analysis on directed acyclic graphs for global optimization. Journal of Global Optimization, 33, 541–562.
- Tawarmalani, M., & Sahinidis, N.V. (2005). A polyhedral branch-and-cut approach to global optimization. Mathematical Programming, 103, 225–249.
- Vu, X., Schichl, H., & Sam-Haroud, D. (2009). Interval propagation and search on directed acyclic
graphs for numerical constraint solving. Journal of Global Optimization, 45, 499–531.
