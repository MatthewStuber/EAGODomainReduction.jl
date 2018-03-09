"""
    STD_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard range reduction using linear underestimator forms via
McCormick relaxations and contracts `X::Vector{Interval{Float64}}` in place.
The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function STD_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)
    # gets prior McCormick relaxation setting
    mu_temp::Int64 = copy(EAGOSmoothMcCormickGrad.MC_param.mu)
    set_diff_relax(0)

    # sets variable bounds
    l::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    u::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = mid.(X)

    # computes relaxation of g and f
    x_mc::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],x0[i],seed_g(opt[1].numVar,i),seed_g(opt[1].numVar,i),X[i],false,X,x0) for i=1:opt[1].numVar]
    f_mc::SMCg{opt[1].numVar,Float64} = opt[1].f(x_mc)
    c::Vector{SMCg{opt[1].numVar,Float64}} = opt[1].g(x_mc)

    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    cval = vcat(c_cv[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-c_cc[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
    dcdx::Array{Float64,2} = vcat([f_mc.cv_grad[i] for i=1:opt[1].numVar]',
                                  dcdx_cv[opt[1].gU_loc,:],
                                  dcdx_cc[opt[1].gL_loc,:])

    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([x0[j]*dcdx_cc[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs = vcat([UBD + f_mc.cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])],
               rhs1[opt[1].gU_loc]+opt[1].gU[opt[1].gU_loc]-c_cv[opt[1].gU_loc],
               rhs2[opt[1].gL_loc]-opt[1].gL[opt[1].gL_loc]+c_cc[opt[1].gL_loc])

    # solve nx lower & upper bounding problems
    for i=1:opt[1].numVar
        # solve lower bounding problems
        c_obj::Vector{Float64} = [j == i ? 1.0 : 0.0 for j=1:opt[1].numVar]
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag::Symbol = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            l[i] = MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
        # solve upper bounding problems
        c_obj[i] = -1.0
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            u[i] = -MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax(mu_temp)

    X[:] = [Interval(l[i],u[i]) for i=1:opt[1].numVar]
    return true
end

"""
    Imp_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard range reduction using linear underestimator forms via
McCormick relaxations and contracts `X::Vector{Interval{Float64}}` in place. However,
it only contracts the p variables. The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc`, `solver.LP_solver`, and `solver.Implicit_Opts.nx`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function Imp_Linear_RR!(X::Vector{Interval{Float64}},opt, UBD::Float64)
    # gets prior McCormick relaxation setting
    mu_temp::Int64 = copy(EAGOSmoothMcCormickGrad.MC_param.mu)
    set_diff_relax(0)

    # sets variable bounds
    nx::Int64 = opt[1].solver.Implicit_Options.nx
    np::Int64 = opt[1].numVar - nx

    l::Vector{Float64} = [X[i].lo for i=(nx+1):(nx+np)]
    u::Vector{Float64} = [X[i].hi for i=(nx+1):(nx+np)]
    x0::Vector{Float64} = (l+u)/2.0

    # computes relaxation of g and f
    p_mc::Vector{SMCg{np,Float64}} = [SMCg{np,Float64}(x0[i],x0[i],seed_g(opt[1].numVar,i),seed_g(opt[1].numVar,i),X[i],false,X,x0) for i=1:opt[1].numVar]
    out,z,x_mc = GenExpansionParams(opt[1].solver.Implicit_Options.h,
                                    opt[1].solver.Implicit_Options.hj,
                                    X[1:nx],
                                    X[(nx+1):(nx+np)],
                                    p0,
                                    opt[1].solver.Implicit_Options.Param)
    f_mc::SMCg{np,Float64} = opt[1].solver.ImplicitOpts.f(x_mc[end],p_SMC)
    c::Vector{SMCg{np,Float64}} = opt[1].solver.ImplicitOpts.g(x_mc[end],p_SMC)

    # UPDATE gL_loc & gU_loc
    gL_temp = filter(x->(x>nx),opt[1].gL_loc)
    gU_temp = filter(x->(x>nx),opt[1].gU_loc)

    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:np]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:np]
    cval = vcat(c_cv[gU_temp-nx]-opt[1].gU[gU_temp-nx],-c_cc[gL_temp-nx]+opt[1].gL[gL_temp-nx])
    dcdx::Array{Float64,2} = vcat([f_mc.cv_grad[i] for i=1:np]',
                                  dcdx_cv[opt[1].gU_loc,:],
                                  dcdx_cc[opt[1].gL_loc,:])

    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:np]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([x0[j]*dcdx_cc[i,j] for j=1:np]) for i=1:opt[1].numConstr]
    rhs = vcat([UBD + f_mc.cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])],
               rhs1[gU_temp-nx]+opt[1].gU[gU_temp]-c_cv[gU_temp-nx],
               rhs2[gL_temp-nx]-opt[1].gL[gL_temp]+c_cc[gL_temp-nx])

    # solve nx lower & upper bounding problems
    for i=1:opt[1].numVar
        # solve lower bounding problems
        c_obj::Vector{Float64} = [j == i ? 1.0 : 0.0 for j=1:np]
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag::Symbol = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            l[i] = MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
        # solve upper bounding problems
        c_obj[i] = -1.0
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            u[i] = -MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax(mu_temp)

    X[(nx+1):(nx+np)] = [Interval(l[i],u[i]) for i=1:np]
    return true
end

"""
    STDSparse_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard range reduction using linear underestimator forms via
McCormick relaxations and contracts `X::Vector{Interval{Float64}}` in place using
sparse routines to calculate affine relaxation of the problem. The `opt` object
has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function STDSparse_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)
    # gets prior McCormick relaxation setting
    mu_temp::Int64 = copy(EAGOSmoothMcCormickGrad.MC_param.mu)
    set_diff_relax(0)

    # sets variable bounds
    l::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    u::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = mid.(X)

    # computes relaxation of g and f
    x_mc::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],x0[i],seed_g(opt[1].numVar,i),seed_g(opt[1].numVar,i),X[i],false,X,x0) for i=1:opt[1].numVar]
    f_mc::SMCg{opt[1].numVar,Float64} = opt[1].f(x_mc)
    c::Vector{SMCg{opt[1].numVar,Float64}} = opt[1].g(x_mc)

    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:opt[1].numVar]
    cval = vcat(c_cv[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-c_cc[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
    dcdx::Array{Float64,2} = vcat([f_mc.cv_grad[i] for i=1:opt[1].numVar]',
                                  dcdx_cv[opt[1].gU_loc,:],
                                  dcdx_cc[opt[1].gL_loc,:])

    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([x0[j]*dcdx_cc[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]
    rhs = vcat([UBD + f_mc.cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])],
               rhs1[opt[1].gU_loc]+opt[1].gU[opt[1].gU_loc]-c_cv[opt[1].gU_loc],
               rhs2[opt[1].gL_loc]-opt[1].gL[opt[1].gL_loc]+c_cc[opt[1].gL_loc])

    # solve nx lower & upper bounding problems
    for i=1:opt[1].numVar
        # solve lower bounding problems
        c_obj::Vector{Float64} = [j == i ? 1.0 : 0.0 for j=1:opt[1].numVar]
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag::Symbol = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            l[i] = MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
        # solve upper bounding problems
        c_obj[i] = -1.0
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            u[i] = -MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax(mu_temp)

    X[:] = [Interval(l[i],u[i]) for i=1:opt[1].numVar]
    return true
end

"""
    ImpSparse_Linear_RR!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard range reduction using linear underestimator forms via
McCormick relaxations and contracts `X::Vector{Interval{Float64}}` in place. This
However, it only contracts the p variables and employs sparse arithmetic. The
`opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc`, `solver.LP_solver`, and `solver.Implicit_Opts.nx`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function ImpSparse_Linear_RR!(X::Vector{Interval{Float64}},opt, UBD::Float64)
    # gets prior McCormick relaxation setting
    mu_temp::Int64 = copy(EAGOSmoothMcCormickGrad.MC_param.mu)
    set_diff_relax(0)

    # sets variable bounds
    nx::Int64 = opt[1].solver.Implicit_Options.nx
    np::Int64 = opt[1].numVar - nx

    l::Vector{Float64} = [X[i].lo for i=(nx+1):(nx+np)]
    u::Vector{Float64} = [X[i].hi for i=(nx+1):(nx+np)]
    x0::Vector{Float64} = (l+u)/2.0

    # computes relaxation of g and f
    p_mc::Vector{SMCg{np,Float64}} = [SMCg{np,Float64}(x0[i],x0[i],seed_g(opt[1].numVar,i),seed_g(opt[1].numVar,i),X[i],false,X,x0) for i=1:opt[1].numVar]
    out,z,x_mc = GenExpansionParams(opt[1].solver.Implicit_Options.h,
                                    opt[1].solver.Implicit_Options.hj,
                                    X[1:nx],
                                    X[(nx+1):(nx+np)],
                                    p0,
                                    opt[1].solver.Implicit_Options.Param)
    f_mc::SMCg{np,Float64} = opt[1].solver.ImplicitOpts.f(x_mc[end],p_SMC)
    c::Vector{SMCg{np,Float64}} = opt[1].solver.ImplicitOpts.g(x_mc[end],p_SMC)

    # UPDATE gL_loc & gU_loc
    gL_temp = filter(x->(x>nx),opt[1].gL_loc)
    gU_temp = filter(x->(x>nx),opt[1].gU_loc)

    c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
    c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
    dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:np]
    dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:np]
    cval = vcat(c_cv[gU_temp-nx]-opt[1].gU[gU_temp-nx],-c_cc[gL_temp-nx]+opt[1].gL[gL_temp-nx])
    dcdx::Array{Float64,2} = vcat([f_mc.cv_grad[i] for i=1:np]',
                                  dcdx_cv[opt[1].gU_loc,:],
                                  dcdx_cc[opt[1].gL_loc,:])

    rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:np]) for i=1:opt[1].numConstr]
    rhs2::Vector{Float64} = [sum([x0[j]*dcdx_cc[i,j] for j=1:np]) for i=1:opt[1].numConstr]
    rhs = vcat([UBD + f_mc.cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])],
               rhs1[gU_temp-nx]+opt[1].gU[gU_temp]-c_cv[gU_temp-nx],
               rhs2[gL_temp-nx]-opt[1].gL[gL_temp]+c_cc[gL_temp-nx])

    # solve nx lower & upper bounding problems
    for i=1:opt[1].numVar
        # solve lower bounding problems
        c_obj::Vector{Float64} = [j == i ? 1.0 : 0.0 for j=1:np]
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag::Symbol = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            l[i] = MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
        # solve upper bounding problems
        c_obj[i] = -1.0
        model = buildlp(c_obj, dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        MathProgBase.optimize!(model)
        stat_flag = MathProgBase.status(model)
        if (stat_flag == :Optimal)
            u[i] = -MathProgBase.getobjval(model)
        elseif (stat_flag == :Infeasible)
            feas = false
            return feas
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax(mu_temp)

    X[(nx+1):(nx+np)] = [Interval(l[i],u[i]) for i=1:np]
    return true
end
