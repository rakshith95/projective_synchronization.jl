function angular_distance(a::AbstractVector, b::AbstractVector)
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

function normalized_euclidean_distance(a::AbstractVector, b::AbstractVector; p=2)
    a = a/norm(a)
    b = b/norm(b)
    return norm(a-b, p)
end

function orthogonal_projection_distance(a::AbstractVector, b::AbstractVector)
    #Orthogonal projection distance 
    
    a = SVector{length(a), eltype(a)}(a/norm(a))
    b = SVector{length(b), eltype(b)}(b/norm(b))
    return norm((SMatrix{length(a),length(a), eltype(a)}(I) - a*a')*b  )
    
end