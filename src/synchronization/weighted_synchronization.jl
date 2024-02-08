function  iteratively_weighted_synchronization(Z::AbstractArray{Projectivity}, synchronization_method::String; averaging_max_it=10, weight_function=cauchy, c=c_cauchy, h=h_robust, error_measure=angular_distance, max_it=100, δ=1e-6)
    exit_loop = false
    iter = 0
    n = size(Z,1)
    dims = 0
    for i=1:n
        for j=1:i
            if Z[i,j].exists
                dims = size(Z[i,j].P,1)
                break
            end
        end
    end
    weights = SMatrix{n,n,Float64}(ones(n,n))
    X_prev = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dims,dims, Float64}(I))], n)) 
    # Q = MMatrix{dims*dims, n, Float64}(zeros(dims*dims, n))
    # Ẑ = SizedMatrix{n,n, Projectivity}(repeat([Projectivity(SMatrix{dims,dims,Float64}(I) )], n, n)) 
    Ẑ = SparseMatrixCSC{Projectivity, Integer}(repeat([Projectivity( SMatrix{dims,dims,Float64}(I) )],n,n)) # Relative projectivities
    prev_weights = zero(weights)
    while !exit_loop && iter < max_it
        if occursin("spectral", lowercase(synchronization_method))
            X = projectivity_synch_spectral(copy(Z), weights)
        else
            if iszero(iter)
                X = iterative_projective_synchronization(copy(Z); averaging_method=synchronization_method, weights=weights, max_iterations=averaging_max_it)
            else
                X = iterative_projective_synchronization(copy(Z); X₀=X_prev, averaging_method=synchronization_method, weights=weights, max_iterations=averaging_max_it)
            end
        end
        # Q_avg = compute_Q(Q, X, X_prev, average)
        compute_Z!(Ẑ, X)

        # How to deal with the fact that error(a, b) !=  error(inv(a), inv(b)). Discuss further
        # E = error_measure.(Z,Ẑ)
        # E = (E_UT + E_LT')/2
        E_UT = error_measure.(UpperTriangular(Z),UpperTriangular(Ẑ))
        E_LT = error_measure.(LowerTriangular(Z),LowerTriangular(Ẑ))
        E = min.(E_UT, E_LT')
        E = E + E'

        s = StatsBase.mad(E[.!isinf.(E)])
        if iszero(s)
            s = std(E[.!isinf.(E)])
        end
        if iszero(s)
            return X
        end
        weights = SMatrix{n,n,Float64}(weight_function.(E/(h*c*s)) )
        iter+=1
        #Check break logic
        # println(mean(compute_err(X, X_prev, SMatrix{4,4,Float64}(I), error_measure)), "\t", δ )
        # rel_weights_diff = abs.(weights - prev_weights)./prev_weights
        # if any(isfinite.(rel_weights_diff)) && median(rel_weights_diff[isfinite.(rel_weights_diff)]) < δ
        # println(rad2deg(mean((compute_err(X, X_prev, SMatrix{4,4,Float64}(I), error_measure)))))
        if mean(compute_err(X, X_prev, SMatrix{4,4,Float64}(I), error_measure)) < δ
            X_prev = X
            exit_loop = true
        end
        X_prev = X
        prev_weights = weights
    end
    # println(iter)
    return X_prev
end
function iteratively_weighted_averaging(M::AbstractArray{T}, average; weight_function=cauchy, c=c_cauchy, h=h_robust, error_measure=angular_distance, max_it=100, δ=1e-5) where T<:AbstractFloat
    exit_loop = false
    iter = 0
    n = size(M,2)
    if n==1
        return M[:,1]
    end
    weights = SVector{n,Float64}(ones(n))
    x_prev = zero(view(M,:,1))
    while !exit_loop
        if iter >= max_it
            break
        end
        x = average(M, weights)
        dists = Vector{Float64}(undef, n)
        for i=1:n
            dists[i] = error_measure(x, view(M,:,i))
        end
        s = StatsBase.mad(dists)
        if iszero(s)
            s = std(dists)
        end
        weights = SVector{n, Float64}(weight_function.(dists/(h*c*s)))
        iter+=1
        if norm(x-x_prev)/norm(x_prev) < δ || iszero(s)
            exit_loop = true
        end
        x_prev = x
    end
    return x_prev
end

# Y = iteratively_weighted_synchronization(copy(Z), "sphere", error_measure=angular_distance, max_it=100, averaging_max_it=5, δ=deg2rad(0.1));
# Q = MMatrix{4*4, 25, Float64}(zeros(16, 25));
# Q_avg = compute_Q(Q, Y, X_gt, spherical_mean);
# rad2deg(median(compute_err(X_gt, Y, Q_avg, angular_distance)))


#  
# X_nr = iterative_projective_synchronization(copy(Z); averaging_method="sphere");
# X_nr = iteratively_weighted_synchronization(copy(Z), "sphere", error_measure=angular_distance, max_it=30, averaging_max_it=30, δ=deg2rad(0.1));
# Q_avg_nr = compute_Q(Q, X_nr, X_gt, weiszfeld);
# rad2deg(median(compute_err(X_gt, X_nr, Q_avg_nr, angular_distance)))

# X_r_inner = iterative_projective_synchronization(copy(Z); averaging_method="robust-sphere");
# Q_avg_r_innter = compute_Q(Q, X_r_inner, X_gt, weiszfeld);
# rad2deg(median(compute_err(X_gt, X_r_inner, Q_avg_r_innter, angular_distance)))
