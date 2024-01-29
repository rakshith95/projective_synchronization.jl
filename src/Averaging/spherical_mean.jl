function spherical_mean(M::AbstractArray;  initialize=nothing, c₀=nothing, max_iterations=1e3, δ=1e-6, σ=1e-4)
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M = M ./ norm.(eachcol(M))'
    if isnothing(c₀)
        if isnothing(initialize)
            c₀ = mean.(eachrow(M))
        else
            c₀ = initialize(M)
        end
    end
    c = c₀
    # if the mean is already a point in the sphere  do nothing
    # this happens with a single point or coincident points 
    if !isapprox(norm(c₀),1.0)
        it=0
        while it < max_iterations
            c_prev = c
            # If current iterate is equal to any of the points, then disturb it by some amount
            if any(norm.(eachcol((M .- vec(c_prev)))) .< 1e-6)
                Δ=rand(Distributions.Normal(0, σ),  size(M)[1])
                c_prev = unit_normalize(c_prev + Δ)
            end    

            c = SVector{length(c_prev)}(zero(c_prev))
            for i in collect(1:size(M,2))
                c = c + view(M,:,i)/ sqrt( 1 - (c_prev'*view(M,:,i))^2 )
            end
            if isapprox(norm(c),0.0)
                c = c_prev
                break 
            end
            c = unit_normalize(c)
            it += 1
            if norm(c-c_prev) < δ 
                break
            end
        end
    end
    return c
end

function spherical_mean(M::AbstractArray, weights::AbstractVector{T};  initialize=nothing, c₀=nothing, max_iterations=1e3, δ=1e-5, σ=1e-4) where T<:AbstractFloat
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M = M ./ norm.(eachcol(M))'
    if isnothing(c₀)
        if isnothing(initialize)
            if !isapprox(sum(weights),1.0)
                weights = weights/sum(weights)
            end
            c₀ = WeightedSum(M, weights)
        else
            c₀ = initialize(M)
        end
    end
    c = c₀
    # if the mean is already a point in the sphere  do nothing
    # this happens with a single point or coincident points 
    it=0
    if !isapprox(norm(c₀),1.0)
        while it < max_iterations
            c_prev = c            
            # If current iterate is equal to any of the points, then disturb it by some amount
            if any(norm.(eachcol((M .- vec(c_prev)))) .< 1e-6)
            # if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
                Δ=rand(Distributions.Normal(0, σ),  size(M)[1])
                c_prev = unit_normalize(c_prev + Δ)
            end    
            c = SVector{length(c_prev)}(zero(c_prev))
            for i in collect(1:size(M,2))
                c = c + (weights[i]*view(M,:,i))/ sqrt( 1 - (c_prev'*view(M,:,i))^2 )
            end
            if isapprox(norm(c),0.0)
                c = c_prev
                break 
            end
            c = unit_normalize(c)
            it += 1
            if norm(c-c_prev) < δ 
                break
            end
        end
    end
    return c
end

function WeightedSum(M::AbstractArray{T}, weights::AbstractVector{T}) where T<:AbstractFloat
    return SVector{size(M,1),T}(sum(weights.*eachcol(M)))
end

function RotateUnitInDirection!(x::AbstractVector{T}, dir::AbstractVector{T}) where T<:AbstractFloat
    θ = norm(dir)
    # If θ=0 , then do nothing. Otherwise,...
    if !isapprox(θ,0.0)
        c = cos(θ)
        s = sin(θ)
        dirUnit = dir/θ
        x[:] = c*x + s*dirUnit
    end
end

function sphere_to_tangent(x::AbstractVector, v::AbstractVector)
    dim = length(x)
    # Take a point v to tangent hyperplane of x 
    c = v'*x #cosθ = Dot prod
    vPerp = v - c*x
    s = norm(vPerp) #sinθ
    if isapprox(s,0)
        return zeros(dim)
    else
        θ = atan(s,c)
        return (θ/s)*vPerp # Take to tangent space. (How does this relate to methodology described in paper?)
    end
end


function weighted_spherical_mean_A1(M::AbstractArray{T}, weights::SVector{N,Float64}; p=1, initialize=nothing, max_iterations=1e3, δ=1e-8) where {N,T}
    # Algorithm A1 in "Spherical Averages and Applications to Spherical Splines and Interpolation" by Samuel R Buss et al. 
    
    #Make columns unit vectors
    M = M ./ norm.(eachcol(M))'
    
    num_pts = size(M,2)
    if !isapprox(sum(weights),1.0)
        weights = weights/sum(weights)
    end

    #Initialize average
    if isnothing(initialize)
        # Try with Weiszfeld
        if p==1
            xVec = weiszfeld(M, weights)
        elseif p==2
            xVec = WeightedSum(M, weights)
            #Project it to unit sphere
            unit_normalize!(xVec)
        end
    else
        xVec = initialize(M, weights)
    end
    
    it=0
    localPoints = zero(M)
    while it <= max_iterations
        xPrev = copy(xVec)
        for i=1:num_pts
            localPoints[:,i] = sphere_to_tangent(xVec, view(M,:,i))
        end
        if p==1
            xDisp = weiszfeld(localPoints, weights, normalize=false)
        elseif p==2
            xDisp = WeightedSum(localPoints, weights) #Avg points in tangent space
        end
        RotateUnitInDirection!(xVec, xDisp) # Take back to sphere by rotating xVec in the direction of xDisp. (How does this relate to methodology described in paper?)
        # unit_normalize!(xVec) #Project to sphere if there are any errs.
        if norm(xVec - xPrev) <= δ
            break
        end
        it += 1
    end
    return xVec
end            

weighted_spherical_mean_A1(M::AbstractArray{T}) where T = weighted_spherical_mean_A1(M, SVector{size(M,2), Float64}(ones(size(M,2))) )
weighted_spherical_mean_A1_L1(M::AbstractArray{T}, weights::SVector{N,Float64}; initialize=nothing, max_iterations=1e3, δ=1e-8) where {N,T} = weighted_spherical_mean_A1(M,weights, p=1, initialize=initialize, max_iterations=max_iterations, δ=δ)
weighted_spherical_mean_A1_L1(M::AbstractArray{T}; initialize=nothing, max_iterations=1e3, δ=1e-8) where T = weighted_spherical_mean_A1(M, SVector{size(M,2), Float64}(ones(size(M,2))) , p=1, initialize=initialize, max_iterations=max_iterations, δ=δ)
weighted_spherical_mean_A1_L2(M::AbstractArray{T}, weights::SVector{N,Float64}; initialize=nothing, max_iterations=1e3, δ=1e-8) where {N,T} = weighted_spherical_mean_A1(M,weights, p=2, initialize=initialize, max_iterations=max_iterations, δ=δ)
weighted_spherical_mean_A1_L2(M::AbstractArray{T}; initialize=nothing, max_iterations=1e3, δ=1e-8) where T = weighted_spherical_mean_A1(M, SVector{size(M,2), Float64}(ones(size(M,2))) , p=2, initialize=initialize, max_iterations=max_iterations, δ=δ)