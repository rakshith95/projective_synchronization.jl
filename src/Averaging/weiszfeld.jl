function weiszfeld(M::AbstractMatrix; normalize=true, max_iterations=1e3, initialize=nothing, c₀=nothing, δ=1e-5, σ=1e-4 )
    # As described in L1 rotation averaging using the Weiszfeld algorithm: Hartley et al.
    #Initialize
    #Make cols of M unit vectors
    if normalize
        M = M ./ norm.(eachcol(M))'
    end

    if isnothing(c₀)
        if isnothing(initialize)
            c₀ = mean.(eachrow(M))
        else
            c₀ = initialize(M)
        end
    end
    
    it=0
    c = c₀
    if !isapprox(norm(c), 1.0)
        while it < max_iterations
            c_prev = c
            # If current iterate is equal to any of the points, then disturb it by some amount
            if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
                Δ=rand(Distributions.Normal(0, σ), size(M)[1])
                c_prev = c_prev + Δ
            end

            # ∇ = zeros(eltype(c_prev), length(c_prev))
            # λ = 0
            # for i in collect(1:size(M,2))
                # ∇ += (M[:,i] - c_prev)/norm(M[:,i]-c_prev,1)
                # λ += norm(M[:,i] - c_prev,1)^(-1)
            # end  
            # c = c_prev + ∇/λ

            N = SVector{length(c_prev), eltype(c_prev)}(zero(c_prev) )
            D = 0
            for i in collect(1:size(M,2))
                N += @views M[:,i] / norm(M[:,i]-c_prev, 1)
                D += 1/norm( view(M,:,i)-c_prev, 1)
            end
            c = N/D
            it+=1
            if norm(c-c_prev,1)/norm(c_prev,1) < δ
                break
            end
        end
    end
    if normalize
        return c/norm(c)
    else
        return c
    end

end


function weiszfeld(M::AbstractMatrix, weights::AbstractVector{T}; normalize=true, max_iterations=1e2, initialize=nothing, c₀=nothing, δ=1e-5, σ=1e-4 ) where T<:AbstractFloat
    # As described in L1 rotation averaging using the Weiszfeld algorithm: Hartley et al.
    #Initialize
    #Make cols of M unit vectors
    if normalize
        M = M ./ norm.(eachcol(M))'
    end
    
    if isnothing(c₀)
        if isnothing(initialize)
            if !isapprox(sum(weights),1.0)
                weights = weights/sum(weights)
            end
            c₀ = WeightedSum(M, weights)
            # c₀ = mean.(eachrow(M))
        else
            c₀ = initialize(M)
        end
    end
    
    it = 0
    c = c₀
    if !isapprox(norm(c),1.0)
        while it < max_iterations
            c_prev = c
            # If current iterate is equal to any of the points, then disturb it by some amount
            if any(isapprox.(norm.(eachcol((M .- vec(c_prev)))),0))
                Δ=rand(Distributions.Normal(0, σ), size(M)[1])
                c_prev = c_prev + Δ
            end
    
            ∇ = SVector{length(c_prev), eltype(c_prev)}(zero(c_prev))
            λ = 0
            for i in collect(1:size(M,2))
                ∇ += weights[i]*((view(M,:,i) - c_prev)/norm(view(M,:,i)-c_prev,1))
                λ += weights[i]*(norm(view(M,:,i) - c_prev,1))^(-1)
            end  
            c = c_prev + ∇/λ
            it+=1
            if norm(c-c_prev,1)/norm(c_prev,1) < δ
                break
            end
        end
    end
    if normalize
        return c/norm(c)
    else
        return c
    end
end