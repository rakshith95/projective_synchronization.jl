import Base: inv, *, +, -, /

struct Projectivity{Q}
    P::Q
    exists::Bool
end

Projectivity(P::Q) where Q = Projectivity(P, true)
Projectivity(b::Bool) = Projectivity(nothing, b)

function unwrap!(Q::Projectivity, P::AbstractMatrix{T}) where T
    if Q.exists
        P[:,:] = @views Q.P
        return true
    else
        return false
    end
end

function unwrap!(Z′::AbstractMatrix{F}, Z::AbstractMatrix{Projectivity}) where F
    n = size(Z,1)
    dims=0
    T = missing
    for i=1:n
        for j=1:n
            if Z[i,j].exists
                dims = size(Z[i,j].P, 1)
                T = eltype(Z[i,j].P)
                break
            end
        end
    end
    if dims==0
        return false
    end

    for i=1:dims:(n*dims)-dims+1
        for j=1:dims:(n*dims)-dims+1
            if Z[div(i,dims)+1,div(j,dims)+1].exists
                Z′[i:i+dims-1, j:j+dims-1] = @views Z[div(i,dims)+1,div(j,dims)+1].P
            end
        end
    end
end

function unit_normalize(P::Projectivity)
    if P.exists
        return Projectivity(P.P/norm(P.P))
    else
        return Projectivity(false)
    end
end

function unit_normalize!(P::Projectivity)
    if P.exists
        return P.P[:,:] = P.P/norm(P.P)
    end
end

function unit_normalize(P::Projectivity, metric)
    if P.exists
        if occursin("det", lowercase(metric))
            d = det(P.P)
            n = size(P.P,1)
            if d>0
                return Projectivity(P.P/(d^(1/n)))
            elseif d<0
                return Projectivity(P.P/(Complex(d)^(1/n)))
            end
        end
    else
        return Projectivity(false)
    end
end

function unit_normalize!(P::Projectivity, metric)
    if P.exists
        if occursin("det", lowercase(metric))
            d = det(P.P)
            n = size(P.P,1)
            if d>0
                P.P[:,:] = P.P/(d^(1/n))
            elseif d<0
                P.P[:,:] = P.P/(Complex(d)^(1/n))
            end
        end
    end
end

function unit_normalize(v::SVector{N,T}) where {N,T}
    return v/norm(v)
end

function unit_normalize(v::AbstractVector{T}) where T
    return v/norm(v)
end

function unit_normalize!(v::AbstractVector{T}) where T
    v[:] = v/norm(v)
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

function *(A::AbstractMatrix, P::Projectivity)
    if P.exists
        return A*P.P
    else
        return A
    end
end

function *(P::Projectivity, A::AbstractMatrix)
    if P.exists
        return P.P*A
    else
        return A
    end
end

function /(A::Projectivity, B::Projectivity)
    if A.exists && B.exists
        return Projectivity(A.P / B.P)
    else
        return Projectivity(false)
    end
end

function inv(P::Projectivity)
    if P.exists
        return Projectivity(inv(P.P))
    else
        return Projectivity(false)
    end
end

function norm(P::Projectivity)
    if P.exists
        return norm(P.P)
    else
        return 0.0
    end
end