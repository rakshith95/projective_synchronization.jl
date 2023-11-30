function anglular_distance(a::AbstractVector, b::AbstractVector)
    a = a/norm(a)
    b = b/norm(b)

    θ = acos(dot(a,b))
    return θ
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