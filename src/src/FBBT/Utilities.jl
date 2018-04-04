function mul_revDR(a::Interval{T},b::Interval{T},c::Interval{T}) where {T<:AbstractFloat}
    if (in(zero(T),c) && ~in(zero(T),b))
        c = c ∩ (a / b)
    elseif (~in(zero(T),c) && in(zero(T),b))
        b = b ∩ (a / c)
    elseif (~in(zero(T),c) && ~in(zero(T),b))
        b = b ∩ (a / c)
        c = c ∩ (a / b)
    end
    return a, b, c
end
function mul_revDR(a::MCInterval{T},b::MCInterval{T},c::MCInterval{T}) where {T<:AbstractFloat}
    if (in(zero(T),c) && ~in(zero(T),b))
        c = c ∩ (a / b)
    elseif (~in(zero(T),c) && in(zero(T),b))
        b = b ∩ (a / c)
    elseif (~in(zero(T),c) && ~in(zero(T),b))
        b = b ∩ (a / c)
        c = c ∩ (a / b)
    end
    return a, b, c
end

function exp_revDR(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
    y_new = y ∩ Interval{T}(zero(T),Inf)
    x_new = x ∩ log(y_new)
    return y_new, x_new
end
function exp_revDR(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
    y_new = y ∩ MCInterval{T}(zero(T),Inf)
    x_new = x ∩ log(y_new)
    return y_new, x_new
end

function acos_rev(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),2.0*pi)
        x_new = x ∩ cos(y_new)

        return y_new, x_new
end
function acos_rev(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),2.0*pi)
        x_new = x ∩ cos(y_new)

        return y_new, x_new
end

function atan_rev(y::Interval{T}, x::Interval{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(-pi/2.0,pi/2.0)
        x_new = x ∩ tan(y_new)

        return y_new, x_new
end
function atan_rev(y::MCInterval{T}, x::MCInterval{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(-pi/2.0,pi/2.0)
        x_new = x ∩ tan(y_new)

        return y_new, x_new
end

for V in (Interval,MCInterval)
    @eval function sinh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        x_new = x ∩ asinh(y)

        return y, x_new
    end

    @eval function cosh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(one(T),∞)
        x_new = x ∩ acosh(y_new)

        return y_new, x_new
    end

    @eval function tanh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(-one(T),one(T))
        x_new = x ∩ atanh(y_new)

        return y_new, x_new
    end

    @eval function asinh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        x_new = x ∩ sinh(y)

        return y, x_new
    end

    @eval function acosh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ cosh(y_new)

        return y_new, x_new
    end

    @eval function atanh_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        x_new = x ∩ tanh(y)

        return y, x_new
    end

    @eval function step_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),one(T))

        return y_new, x
    end

    @eval function sign_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(-one(T),one(T))

        return y_new, x
    end

    @eval function exp2_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ log2(y_new)

        return y_new, x_new
    end

    @eval function exp10_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ log10(y_new)

        return y_new, x_new
    end

    @eval function log2_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ exp2(y_new)

        return y_new, x_new
    end

    @eval function log10_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ exp10(y_new)

        return y_new, x_new
    end

    @eval function one_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(one(T))

        return y_new, x
    end

    @eval function zero_rev(y::$V{T},x::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T))

        return y_new, x
    end

    @eval function min_rev(z::$V{T},x::$V{T},y::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ log(y_new)

        return y_new, x_new
    end
    @eval function max_rev(z::$V{T},x::$V{T},y::$V{T}) where {T<:AbstractFloat}
        y_new = y ∩ $(V{T})(zero(T),∞)
        x_new = x ∩ log(y_new)

        return y_new, x_new
    end
end

mul_revDR(a,b,c) = div_revDR(promote(a,b,c)...)
exp_revDR(a,b) = exp_revDR(promote(a,b)...)
sinh_rev(a,b) = sinh_rev(promote(a,b)...)
cosh_rev(a,b) = cosh_rev(promote(a,b)...)
tanh_rev(a,b) = tanh_rev(promote(a,b)...)
asinh_rev(a,b) = asinh_rev(promote(a,b)...)
acosh_rev(a,b) = acosh_rev(promote(a,b)...)
atanh_rev(a,b) = atanh_rev(promote(a,b)...)
step_rev(a,b) = step_rev(promote(a,b)...)
sign_rev(a,b) = sign_rev(promote(a,b)...)
exp2_rev(a,b) = exp2_rev(promote(a,b)...)
exp10_rev(a,b) = exp10_rev(promote(a,b)...)
log2_rev(a,b) = log2_rev(promote(a,b)...)
log10_rev(a,b) = log10_rev(promote(a,b)...)
one_rev(a,b) = one_rev(promote(a,b)...)
zero_rev(a,b) = zero_rev(promote(a,b)...)
