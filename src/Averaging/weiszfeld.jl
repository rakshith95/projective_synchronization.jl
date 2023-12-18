function weiszfeld(M::AbstractMatrix; max_iterations=1e3, initialize=nothing, δ=1e-10, σ=1e-3, Δ=rand(Distributions.Normal(0, σ), size(M)[1]) )
    # As described in L1 rotation averaging using the Weiszfeld algorithm: Hartley et al.
    #Initialize
    #Make cols of M unit vectors
    M = M ./ norm.(eachcol(M))'

    if isnothing(initialize)
        c₀ = mean.(eachrow(M))
    else
        c₀ = initialize(M)
    end
    
    it=0
    c = c₀
    while it < max_iterations
        c_prev = c
        # If current iterate is equal to any of the points, then disturb it by some amount
        if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
            c_prev = c_prev + Δ
        end

        # ∇ = zeros(eltype(c_prev), length(c_prev))
        # λ = 0
        # for i in collect(1:size(M,2))
            # ∇ += (M[:,i] - c_prev)/norm(M[:,i]-c_prev,1)
            # λ += norm(M[:,i] - c_prev,1)^(-1)
        # end  
        # c = c_prev + ∇/λ

        N = zeros(eltype(c_prev), length(c_prev))
        D = 0
        for i in collect(1:size(M,2))
            N += M[:,i] / norm(M[:,i]-c_prev, 1)
            D += 1/norm(M[:,i]-c_prev, 1)
        end
        c = N/D
        c = c/norm(c)
        it+=1
        if norm(c-c_prev,1) < δ
            break
        end
    end
    return c
end