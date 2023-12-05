import Base: inv, *, +, - 

struct Projectivity{Q}
    P::Q
    exists::Bool
end

Projectivity(P::Q) where Q = Projectivity(P, true)
Projectivity(b::Bool) = Projectivity(nothing, b)

function unit_normalize(P::Projectivity)
    if P.exists
        return Projectivity(P.P/norm(P.P))
    else
        return Projectivity(false)
    end
end

#Overload operators for the Projectivity struct

function Base.zero(P::Projectivity{Q}) where Q
    return !P.exists
end

function +(P::Projectivity, Q::Projectivity)
    if P.exists && Q.exists
        return Projectivity(P.P + Q.P)
    else
        return Projectivity(false)
    end
end

function -(P::Projectivity, Q::Projectivity)
    if P.exists && Q.exists
        return Projectivity(P.P - Q.P)
    else
        return Projectivity(false)
    end
end

function *(P::Projectivity, Q::Projectivity)
    if P.exists && Q.exists
        return Projectivity(P.P*Q.P)
    else
        return Projectivity(false)
    end

end

function *(s::T , P::Projectivity) where T<:AbstractFloat
    if P.exists && !iszero(s)
        return Projectivity(s*P.P)
    else
        return Projectivity(false)
    end
end

function *(P::Projectivity, s::T) where T<:AbstractFloat
    *(s,P)
end

function inv(P::Projectivity)
    if P.exists
        return Projectivity(inv(P.P))
    else
        return Projectivity(false)
    end
end

