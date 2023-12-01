function spherical_mean(M::AbstractArray;  initialize=nothing, max_iterations=1e3, δ=1e-10)
    # Columns of M are unit vectors which need to be averaged 
    # initialize
    if isnothing(initialize)
        c₀ = mean.(eachrow(M./norm(eachcol(M)')))
    else
        c₀ = initialize(M)
    end
    # Make unit vector
    if !isapprox(norm(c₀),1.0)
        c = c₀/norm(c₀)
    else
        c = c₀
    end
    it=0
    while it < max_iterations
        c_prev = c
        c = zeros(length(c_prev))
        for i in collect(1:size(M)[2])
            if dot(c, M[:,i]) < 1.0
                c = c + M[:,i]/sqrt(1 - dot(c,M[:,i])^2 )
            end
        end
        c = c/norm(c)
        it += 1
        if norm(c-c_prev) < δ 
            break
        end
    end

    return c
end

# A = SMatrix{3,3,Float64}(rand(3,3))
# normalize_column(x) = x/norm(x)
# B = mapslices(normalize_column, A, dims=1)
# spherical_mean(A)

