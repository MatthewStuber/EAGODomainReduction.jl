#__precompile__()

module EAGODomainReduction

using MathProgBase
using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using IntervalContractors

# import basic functions to overload
import Base: +, -, /, *, ^, exp, exp2, exp10, log, log2, log10, sign, step, abs,
             sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh,
             atanh, one, zero, sqrt, min, max, convert, promote_rule

# replace with using once IntervalContractor updated
import IntervalContractors: plus_rev, minus_rev, mul_rev, power_rev, sqr_rev,
                            sqrt_rev, abs_rev, sin_rev, cos_rev, tan_rev,
                            asin_rev, log_rev, exp_rev, inv_rev

# exports main FBBT functions
export Generate_Tape, Generate_TapeList, Generate_Fixed_Tape, Generate_Fixed_TapeList,
       FFBT_Refine, SetConstraint!, SetVarBounds!, GetVarBounds, ForwardPass!,
       ReversePass!, DAGContractor!, NodeFinder, Tape, TapeList, getDAG

# exports supplemental reverse contractors (push to IntervalContractor Library)
export tanh_rev, div_revDR, acos_rev, atan_rev, sinh_rev, cosh_rev, tanh_rev,
       asinh_rev, acosh_rev, atanh_rev, min_rev, max_rev, step_rev, sign_rev,
       exp2_rev, exp10_rev, log2_rev, log10_rev, one_rev, zero_rev

# exports OBBT functions
export Variable_DR!, Variable_DR_Imp!,
       STD_Linear_RR!,
       STD_LP_Probe! #LP_contractor, poorLP_contractor

# includes utility subroutines
include("src/utils.jl")

# includes subroutines for optimality-based bound tightening
include("src/OBBT/Variable_Dual.jl")
include("src/OBBT/Linear_RR.jl")
include("src/OBBT/STD_LP_Probe.jl")
include("src/OBBT/poorLP.jl")

# includes subroutines for feasibility-based bound tightening
include("src/FBBT/NodeFinder.jl")
include("src/FBBT/DAGprop.jl")
include("src/FBBT/Utilities.jl")

#=
function __init__()
end
=#

end # module
