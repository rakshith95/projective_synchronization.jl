function spherical_mean(M::AbstractArray;  initialize=nothing, max_iterations=1e3, δ=1e-10)
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    M = M ./ norm.(eachcol(M))'
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
                if c_prev'*M[:,i] < 1.0
                    c = c + M[:,i]/sqrt( 1 - (c_prev'*M[:,i])^2 )
                end
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
