function angular_distance(a::AbstractVector{T}, b::AbstractVector{T}) where T<:AbstractFloat
    a = a/norm(a)
    b = b/norm(b)
    dotval = dot(a,b)
    if isapprox(dotval, 1.0)
        dotval = 1.0
    elseif isapprox(dotval, -1.0)
        dotval = -1.0
    end
    θ = acos(dotval)
    return min(θ, π-θ)
end
angular_distance(A::Projectivity, B::Projectivity) = A.exists && B.exists ? angular_distance(vec(A.P), vec(B.P)) : Inf
angular_distance(a::Bool, b::Bool) = a && b ? 1.0 : 0.0

function normalized_euclidean_distance(a::AbstractVector, b::AbstractVector; p=2)
    a = a/norm(a)
    b = b/norm(b)
    return norm(a-b, p)
end
normalized_euclidean_distance(A::Projectivity, B::Projectivity) = A.exists && B.exists ? normalized_euclidean_distance(vec(A.P), vec(B.P)) : Inf
normalized_euclidean_distance(a::Bool, b::Bool) = a && b ? 1.0 : 0.0

function orthogonal_projection_distance(a::AbstractVector, b::AbstractVector)
    #Orthogonal projection distance 
    
    a = SVector{length(a), eltype(a)}(a/norm(a))
    b = SVector{length(b), eltype(b)}(b/norm(b))
    return norm((SMatrix{length(a),length(a), eltype(a)}(I) - a*a')*b  ) 
end
orthogonal_projection_distance(A::Projectivity, B::Projectivity) = A.exists && B.exists ? orthogonal_projection_distance(vec(A.P), vec(B.P)) : Inf
orthogonal_projection_distance(a::Bool, b::Bool) = a && b ? 1.0 : 0.0