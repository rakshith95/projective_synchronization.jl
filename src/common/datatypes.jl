struct Projectivity{Q}
    P::Q
    exists::Bool
end

Projectivity(P::Q) where Q = Projectivity(P, true)
Projectivity(b::Bool) = Projectivity(nothing, b)

function Base.zero(P::Projectivity{Q}) where Q
    return !P.exists
end


# S = SMatrix{2,2, Projectivity}([ Projectivity(R) Projectivity(false); Projectivity(3*R) Projectivity(4*R)])
# SP = sparse(S)
