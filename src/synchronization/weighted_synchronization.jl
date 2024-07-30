function compute_weights(Z::AbstractMatrix, Ẑ::AbstractMatrix;error_measure=angular_distance, weight_function=cauchy, c=c_cauchy, h=h_robust)
    E_UT = error_measure.(UpperTriangular(Z),UpperTriangular(Ẑ))
    E_LT = error_measure.(LowerTriangular(Z),LowerTriangular(Ẑ))
    E = min.(E_UT, E_LT')
    E = E + E'
    s = StatsBase.mad(E[.!isinf.(E)])
    if iszero(s)
        s = std(E[.!isinf.(E)])
    end
    if iszero(s)
        return false
    end
    wts = weight_function.(E/(h*c*s)) 
end

function  iteratively_weighted_synchronization(Z::AbstractArray{Projectivity}, synchronization_method::String;weights=nothing, X₀=nothing, averaging_max_it=10, weight_function=cauchy, c=c_cauchy, h=h_robust, error_measure=angular_distance, max_it=100, δ_irls=1e-6, kwargs...)
    max_it_init = get(kwargs, :averaging_max_it_init, 100)
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
    if isnothing(weights)
        if n*dims<=100
            wts = SMatrix{n,n,Float64}(ones(n,n))
        else
            wts = ones(n,n)
        end
    else
        wts = weights
    end
    X_prev = sparse(repeat([Projectivity(SMatrix{dims,dims, Float64}(I))], n)) 
    Ẑ = SparseMatrixCSC{Projectivity, Integer}(repeat([Projectivity(false)],n,n)) # Relative projectivities
    while !exit_loop && iter < max_it
        if iszero(iter)
            if isnothing(X₀)
                X = iterative_projective_synchronization(copy(Z); averaging_method=synchronization_method, weights=wts, max_iterations=max_it_init , kwargs...)
            else
                X = iterative_projective_synchronization(copy(Z); X₀=copy(X₀), averaging_method=synchronization_method, weights=wts, max_iterations=max_it_init , kwargs...)
            end
        else
            X = iterative_projective_synchronization(copy(Z); X₀=copy(X_prev), min_updates=0, averaging_method=synchronization_method, weights=wts, max_iterations=averaging_max_it, kwargs...)
        end
        compute_Z!(Ẑ, X)
        # How to deal with the fact that error(a, b) !=  error(inv(a), inv(b)). Discuss further
        wts_prev = wts
        wts = compute_weights(Z, Ẑ, error_measure=error_measure, weight_function=weight_function, c=c, h=h)
        if typeof(wts) == Bool
            return X, wts_prev
        end
        iter+=1
        #Check break logic
        if mean(compute_err(X, X_prev, SMatrix{4,4,Float64}(I), error_measure)) <= δ_irls
            X_prev = X
            exit_loop = true
        end
        X_prev = X
    end
    # println(iter)
    return X_prev, wts
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

# Y, wts = iteratively_weighted_synchronization(Z, "sphere", weight_function=welsch, c=c_welsch, error_measure=angular_distance, max_it=20, averaging_max_it=10, δ_irls=deg2rad(0.1), anchor="nothing", update="start-R-all-random");

# Ẑ = SparseMatrixCSC{Projectivity, Integer}(repeat([Projectivity(SMatrix{4,4,Float64}(I) )],25,25)) # Relative projectivities
# compute_Z!(Ẑ, X_gt);
# E_UT = angular_distance.(UpperTriangular(Z),UpperTriangular(Ẑ))
# E_LT = angular_distance.(LowerTriangular(Z),LowerTriangular(Ẑ))
# E = max.(E_UT, E_LT')
# E = E + E'
# s = StatsBase.mad(E[.!isinf.(E)])
# s = std(E[.!isinf.(E)])
# err_weights = SMatrix{25,25,Float64}(cauchy.(E/(1.0*c_cauchy*s)) )


# i = rand(1:25)
# j = rand(1:25)
# for ind in outies
#     i,j = ind[1], ind[2]
#     if Z[i,j].exists
#         println(wts[i,j])
#         # if !(ind in outies) && !(CartesianIndex(j,i) in outies)
#             # print(ind,"\t")
#             # println(wts[ind])
#         # end
#     end
# end
# for ind in findall(wts .<= 1.0 )#outies
#     i,j = ind[1], ind[2]
#     if Z[i,j].exists
#         # println(weights[i,j])
#         if !(ind in outies) && !(CartesianIndex(j,i) in outies)
#             # println(ind in outies)
#             print(ind,"\t")
#             println(wts[ind])
#         end
#     end
# end
