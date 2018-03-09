"""
    STD_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard probing using linear underestimator forms via
McCormick relaxations, sparse calculations, and contracts `X::Vector{Interval{Float64}}` in place.
The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function STD_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

    # constructs LP relaxation
    Xlo::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    Xhi::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = (Xlo + Xhi)/2.0
    x_SMC::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],
                                                                    x0[i],
                                                                    seed_g(Float64,opt[1].numVar,i),
                                                                    seed_g(Float64,opt[1].numVar,i),
                                                                    X[i],
                                                                    false,
                                                                    X,
                                                                    x0) for i=1:opt[1].numVar]

    # probes upper bound
    f::SMCg{opt[1].numVar,Float64} = opt[1].f(x_SMC)
    f_cv::Float64 = f.cv
    c = opt[1].g(x_SMC)
    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    cval = vcat(c_cv[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-c_cc[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
    dcdx::Array{Float64,2} = vcat(dcdx_cv[opt[1].gU_loc,:],dcdx_cc[opt[1].gL_loc,:])
    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([-x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs = vcat(rhs1[opt[1].gU_loc]+opt[1].gU[opt[1].gU_loc]-c_cv[opt[1].gU_loc],
               rhs2[opt[1].gL_loc]-opt[1].gL[opt[1].gL_loc]+c_cc[opt[1].gL_loc])

    if (opt[1].numConstr>0)
        temp_model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, Xlo, Xhi, opt[1].solver.LP_solver)
    else
        temp_model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], zeros(opt[1].numVar,opt[1].numVar), '<', zeros(opt[1].numVar), Xlo, Xhi, opt[1].solver.LP_solver)
    end
    for i=1:opt[1].numVar

        # probes upper bound
        setvarLB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        setvarUB!(temp_model,Xhi)
        result = solvelp(temp_model)
        val = result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar])
        mult::Vector{Float64} = result.attrs[:redcost]
        mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        Variable_DR!(X,mult_lo,mult_hi,LBD,UBD)

        # probes lower bound
        setvarLB!(temp_model,Xlo)
        setvarUB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        result = solvelp(temp_model)
        val = result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar])
        mult = result.attrs[:redcost]
        mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        Variable_DR!(X,mult_lo,mult_hi,val,UBD)

    end
end

"""
    Imp_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard probing using linear underestimator forms via
McCormick relaxations, sparse calculations, and contracts `X::Vector{Interval{Float64}}`
in place only tightening the following intervals in the vector `X[(nx+1):(nx+np)]`.
The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function Imp_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

    # constructs LP relaxation
    Xlo::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    Xhi::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = (Xlo + Xhi)/2.0
    x_SMC::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],
                                                                    x0[i],
                                                                    seed_g(Float64,opt[1].numVar,i),
                                                                    seed_g(Float64,opt[1].numVar,i),
                                                                    X[i],
                                                                    false,
                                                                    X,
                                                                    x0) for i=1:opt[1].numVar]

    # probes upper bound
    f::SMCg{opt[1].numVar,Float64} = opt[1].f(x_SMC)
    f_cv::Float64 = f.cv
    c = opt[1].g(x_SMC)
    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    cval = vcat(c_cv[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-c_cc[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
    dcdx::Array{Float64,2} = vcat(dcdx_cv[opt[1].gU_loc,:],dcdx_cc[opt[1].gL_loc,:])
    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([-x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs = vcat(rhs1[opt[1].gU_loc]+opt[1].gU[opt[1].gU_loc]-c_cv[opt[1].gU_loc],
               rhs2[opt[1].gL_loc]-opt[1].gL[opt[1].gL_loc]+c_cc[opt[1].gL_loc])

    if (opt[1].numConstr>0)
        temp_model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, Xlo, Xhi, opt[1].solver.LP_solver)
    else
        temp_model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], zeros(opt[1].numVar,opt[1].numVar), '<', zeros(opt[1].numVar), Xlo, Xhi, opt[1].solver.LP_solver)
    end
    for i=1:opt[1].numVar

        # probes upper bound
        setvarLB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        setvarUB!(temp_model,Xhi)
        result = solvelp(temp_model)
        val = result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar])
        mult::Vector{Float64} = result.attrs[:redcost]
        mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        Variable_DR!(X,mult_lo,mult_hi,LBD,UBD)

        # probes lower bound
        setvarLB!(temp_model,Xlo)
        setvarUB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        result = solvelp(temp_model)
        val = result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar])
        mult = result.attrs[:redcost]
        mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        Variable_DR!(X,mult_lo,mult_hi,val,UBD)

    end
end
