function mul_revDR(a::Interval,b::Interval,c::Interval) # DONE
    if (in(0,c) && ~in(0,b))
        c = c ∩ (a / b)
    elseif (~in(0,c) && in(0,b))
        b = b ∩ (a / c)
    elseif (~in(0,c) && ~in(0,b))
        b = b ∩ (a / c)
        c = c ∩ (a / b)
    end
    return a, b, c
end
mul_revDR(a,b,c) = div_revDR(promote(a,b,c)...)

function exp_revDR(y::Interval, x::Interval) # DONE
    y_new = y ∩ (0..∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end
exp_revDR(a,b) = exp_revDR(promote(a,b)...)

function acos_rev(y::Interval, x::Interval) # DONE
    y_new = y ∩ (0.0..2*pi)
    x_new = x ∩ cos(y_new)

    return y_new, x_new
end
acos_rev(a,b) = acos_rev(promote(a,b)...)

function atan_rev(y::Interval, x::Interval) # DONE
    y_new = y ∩ (-pi/2.0..pi/2.0)
    x_new = x ∩ tan(y_new)

    return y_new, x_new
end
atan_rev(a,b) = atan_rev(promote(a,b)...)

function sinh_rev(y::Interval,x::Interval) # DONE
    x_new = x ∩ asinh(y_new)

    return y, x_new
end
sinh_rev(a,b) = sinh_rev(promote(a,b)...)

function cosh_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (1.0..∞)
    x_new = x ∩ acosh(y_new)

    return y_new, x_new
end
cosh_rev(a,b) = cosh_rev(promote(a,b)...)

function tanh_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ ((-1.0)..(1.0))
    x_new = x ∩ atanh(y_new)

    return y_new, x_new
end
tanh_rev(a,b) = tanh_rev(promote(a,b)...)

function asinh_rev(y::Interval,x::Interval) # DONE
    x_new = x ∩ sinh(y_new)

    return y, x_new
end
asinh_rev(a,b) = asinh_rev(promote(a,b)...)

function acosh_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0.0..∞)
    x_new = x ∩ cosh(y_new)

    return y_new, x_new
end
acosh_rev(a,b) = acosh_rev(promote(a,b)...)

function atanh_rev(y::Interval,x::Interval) # DONE
    x_new = x ∩ tanh(y_new)

    return y, x_new
end
atanh_rev(a,b) = atanh_rev(promote(a,b)...)

function step_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0..1)

    return y_new, x
end
step_rev(a,b) = step_rev(promote(a,b)...)

function sign_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ ((-1.0)..(1.0))

    return y_new, x
end
sign_rev(a,b) = sign_rev(promote(a,b)...)

function exp2_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0..∞)
    x_new = x ∩ log10(y_new)

    return y_new, x_new
end
exp2_rev(a,b) = exp2_rev(promote(a,b)...)

function exp10_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0..∞)
    x_new = x ∩ log10(y_new)

    return y_new, x_new
end
exp10_rev(a,b) = exp10_rev(promote(a,b)...)

function log2_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0..∞)
    x_new = x ∩ exp2(y_new)

    return y_new, x_new
end
log2_rev(a,b) = log2_rev(promote(a,b)...)

function log10_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ (0..∞)
    x_new = x ∩ exp10(y_new)

    return y_new, x_new
end
log10_rev(a,b) = log10_rev(promote(a,b)...)

function one_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ ((1.0)..(1.0))

    return y_new, x
end
one_rev(a,b) = one_rev(promote(a,b)...)

function zero_rev(y::Interval,x::Interval) # DONE
    y_new = y ∩ ((0.0)..(0.0))

    return y_new, x
end
zero_rev(a,b) = zero_rev(promote(a,b)...)

function min_rev(z::Interval,x::Interval,y::Interval)
    y_new = y ∩ (0..∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end
function max_rev(z::Interval,x::Interval,y::Interval)
    y_new = y ∩ (0..∞)
    x_new = x ∩ log(y_new)

    return y_new, x_new
end
