module TestDAGcntr

using Compat
using Compat.Test
using IntervalArithmetic
using EAGODomainReduction

single_expr = :(x[1]+x[2]^2-exp(x[1]))
tape_out = Generate_Tape(single_expr,2,-2,4)
dag_out = getDAG(single_expr,[1,2])
@test dag_out[1] == Int64[1; 2; 3; 4; 5; 6; 7]
@test dag_out[3][4] == Int64[4; 5]

#=
generates the tape list for the constraints -1 <= (x[1]+x[2]^2-exp(x[1]) <= 2
and 2 <= :(x[1]*x[2]-x[5]) <= 4 and runs the contractor 6 times in place on the
variable box X = [1,20]^5
=#
mult_expr = [:(x[1]+x[2]^2)
             :(x[1]+x[5])
            ]
tapelist_out = Generate_TapeList(mult_expr,5,[-1.0,2.0],[2.0,3.0])
X = [Interval(0.5,10.0) for i=1:5]
DAGContractor!(X,tapelist_out,6)
@test 0.5-1E-4 <= X[1].lo <= 0.5+1E-4
@test 0.5-1E-4 <= X[2].lo <= 0.5+1E-4
@test 0.5-1E-4 <= X[3].lo <= 0.5+1E-4
@test 0.5-1E-4 <= X[4].lo <= 0.5+1E-4
@test 0.5-1E-4 <= X[5].lo <= 0.5+1E-4
@test 1.75-1E-4 <= X[1].hi <= 1.75+1E-4
@test 1.22475-1E-4 <= X[2].hi <= 1.22475+1E-4
@test 10-1E-4 <= X[3].hi <= 10+1E-4
@test 10-1E-4 <= X[4].hi <= 10+1E-4
@test 2.5-1E-4 <= X[5].hi <= 2.5+1E-4

end
