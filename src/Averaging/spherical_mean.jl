function spherical_mean(M::AbstractArray;  initialize=nothing, c₀=nothing, max_iterations=1e3, δ=1e-6, σ=1e-4)
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M[:,:] = M ./ norm.(eachcol(M))'
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
            if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
                Δ=rand(Distributions.Normal(0, σ),  size(M)[1])
                c_prev = c_prev + Δ
            end    

            c = zeros(length(c_prev))
            for i in collect(1:size(M,2))
                c = c + M[:,i]/ sqrt( 1 - (c_prev'*M[:,i])^2 )
            end
            if isapprox(norm(c),0.0)
                c = c_prev
                break 
            end
            c = c/norm(c)
            it += 1
            if norm(c-c_prev) < δ 
                break
            end
        end
    end
    return c
end

function spherical_mean(M::AbstractArray, weights::AbstractVector{T};  initialize=nothing, c₀=nothing, max_iterations=1e3, δ=1e-6, σ=1e-4 ) where T<:AbstractFloat
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M[:,:] = M ./ norm.(eachcol(M))'
    if isnothing(c₀)
        if isnothing(initialize)
            if !isapprox(sum(weights),1.0)
                weights[:] = weights/sum(weights)
            end
            c₀ = WeightedSum(M, weights)
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
            if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
                Δ=rand(Distributions.Normal(0, σ),  size(M)[1])
                c_prev = c_prev + Δ
            end    
            
            c = zeros(length(c_prev))
            for i in collect(1:size(M,2))
                c = c + (weights[i]*M[:,i])/ sqrt( 1 - (c_prev'*M[:,i])^2 )
            end
            if isapprox(norm(c),0.0)
                c = c_prev
                break 
            end
            c = c/norm(c)
            it += 1
            if norm(c-c_prev) < δ 
                break
            end
        end
    end
    return c
end

function WeightedSum(M::AbstractArray{T}, weights::AbstractVector{T}) where T<:AbstractFloat
    return sum(weights.*eachcol(M))
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


function weighted_spherical_mean_A1(M::AbstractArray{T}; p=1, weights=nothing, initialize=nothing, max_iterations=1e3, δ=1e-8) where T
    # Algorithm A1 in "Spherical Averages and Applications to Spherical Splines and Interpolation" by Samuel R Buss et al. 
    
    #Make columns unit vectors
    M[:,:] = M ./ norm.(eachcol(M))'
    
    num_pts = size(M,2)
    if isnothing(weights)
        weights = SVector{num_pts, T}(ones(num_pts)*(1/num_pts))
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
            localPoints[:,i] = sphere_to_tangent(xVec, M[:,i])
        end
        if p==1
            xDisp = weiszfeld(localPoints, weights, normalize=false)
        elseif p==2
            xDisp = WeightedSum(localPoints, weights) #Avg points in tangent space
        end
        RotateUnitInDirection!(xVec, xDisp) # Take back to sphere by rotating xVec in the direction of xDisp. (How does this relate to methodology described in paper?)
        unit_normalize!(xVec) #Project to sphere if there are any errs.
        if norm(xVec - xPrev) <= δ
            break
        end
        it += 1
    end
    return xVec
end            