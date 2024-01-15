function spherical_mean(M::AbstractArray;  initialize=nothing, max_iterations=1e3, δ=1e-6)
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M[:,:] = M ./ norm.(eachcol(M))'
    if isnothing(initialize)
        c₀ = mean.(eachrow(M))
    else
        c₀ = initialize(M)
    end

    c = c₀
    # if the mean is already a point in the sphere  do nothing
    # this happens with a single point or coincident points 
    if !isapprox(norm(c₀),1.0)
        it=0
        while it < max_iterations
            c_prev = c
            c = zeros(length(c_prev))
            for i in collect(1:size(M,2))
                if (c_prev'*M[:,i]) < 1.0
                    if isapprox(abs(c_prev'*M[:,i]), 1.0)
                        c = c + M[:,i]
                    else 
                        c = c + M[:,i]/ sqrt( 1 - (c_prev'*M[:,i])^2 )
                    end
                end
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

function weighted_spherical_mean_A1(M::AbstractArray{T}; weights=nothing, initialize=nothing, max_iterations=1e3, δ=1e-8) where T
    # Algorithm A1 in "Spherical Averages and Applications to Spherical Splines and Interpolation" by Samuel R Buss et al. 
    
    #Make columns unit vectors
    M[:,:] = M ./ norm.(eachcol(M))'
    
    num_pts = size(M,2)
    dim = size(M,1)
    if isnothing(weights)
        weights = SVector{num_pts, T}(ones(num_pts)*(1/num_pts))
    end

    #Initialize average
    if isnothing(initialize)
        xVec = WeightedSum(M, weights)
    else
        xVec = initialize(M, weights)
    end
    #Project it to unit sphere
    unit_normalize!(xVec)
    
    it=0
    localPoints = copy(M)
    while it <= max_iterations
        xPrev = copy(xVec)
        for i=1:num_pts
            c = M[:,i]'*xVec #cosθ = Dot prod
            vPerp = M[:,i] - c*xVec
            s = norm(vPerp) #sinθ
            if isapprox(s,0)
                localPoints[:,i] = zeros(dim)
            else
                θ = atan(s,c)
                localPoints[:,i] = (θ/s)*vPerp # Take to tangent space. (How does this relate to methodology described in paper?)
            end
        end
        xDisp = WeightedSum(localPoints, weights) #Avg points in tangent space
        RotateUnitInDirection!(xVec, xDisp) # Take back to sphere by rotating xVec in the direction of xDisp. (How does this relate to methodology described in paper?)
        unit_normalize!(xVec) #Project to sphere if there are any errs.
        if norm(xVec - xPrev) <= δ
            break
        end
        it += 1
    end
    return xVec
end            