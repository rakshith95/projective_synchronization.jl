function create_synthetic(σ;error=orthogonal_projection_distance,average=spherical_mean, averaging_methods=["sphere", "dyadic", "least-squares-orthogonal", "gaia", "crossproduct", "euclidean", "weiszfeld" ], kwargs...)
    normalize_matrix = get(kwargs, :normalize_matrix, false)
    dims = get(kwargs, :dimension, 4)
    n = get(kwargs, :frames, 25)
    ρ = get(kwargs, :holes_density, 0.5) 
    Ρ = get(kwargs, :outliers, 0.0)

    X_gt = SizedVector{n, Projectivity}(repeat([Projectivity(false)],n)) # Ground truth nodes with dimxdim matrices
    for i=1:n
        Xᵢ = rand(dims,dims)
        X_gt[i] = Projectivity(SMatrix{dims,dims,Float64}(Xᵢ/norm(Xᵢ)))
    end

    Z = SizedMatrix{n,n,Projectivity}(repeat([Projectivity(false)],n,n)) # Relative projectivities
    for i=1:n
        for j=1:i
            if i==j
                continue
            end
            Z_ij = X_gt[i]*inv(X_gt[j]) + Projectivity(rand(Distributions.Normal(0, σ), dims, dims)) # Add noise
            Z[i,j] =  Z_ij
            Z[j,i] = inv(Z_ij) # Symmetric block is inverse
        end
    end
    
    # Make holes
    A = sprand(n, n, ρ)
    A[A.!=0] .= 1.0
    A = sparse(ones(n,n)) - A
    A = triu(A,1) + triu(A,1)'
    Z = Z.*A
    
    # outliers
    num_UT =  length(findall(triu(A,1) .!= 0))
    num_outliers = Int(round(Ρ*num_UT))
    UT_outliers = StatsBase.sample( findall(triu(A,1).!=0), num_outliers, replace=false  )
    for ind in UT_outliers
        Z[ind] = Projectivity(rand(dims,dims))
        Z[CartesianIndex(reverse(ind.I))] = inv(Z[ind])
    end

    if normalize_matrix
        Z = unit_normalize.(Z)
    end
    X_sol_spectral = projectivity_synch_spectral(copy(Z))
    # X_solᵢ*Q = Xᵢ
    # Either take Q = Xᵢ for the anchor node i, OR
    # Take avg of Qᵢ = inv(X_solᵢ)*Xᵢ ∀ i ∈ 1..n

    Q = MMatrix{dims*dims, n, Float64}(zeros(dims*dims, n))
    for i=1:n
        Q[:,i] = vec((inv(X_sol_spectral[i])*X_gt[i]).P)
    end
    Q_avg = SMatrix{dims,dims,Float64}(reshape(average(Q),dims,dims))
    err = zeros(n)
    for i=1:n
        xgtᵢ = vec(@views X_gt[i].P)
        xsolᵢ = vec(@views X_sol_spectral[i] * Q_avg)
        err[i] = error(xgtᵢ, xsolᵢ)
    end
    
    for method in averaging_methods
        X_sol_iterative = iterative_projective_synchronization(copy(Z);X₀=X_sol_spectral, averaging_method=method,kwargs...)
        for i=1:n
            Q[:,i] = vec((inv(X_sol_iterative[i])*X_gt[i]).P)
        end
        Q_avg = SMatrix{dims,dims,Float64}(reshape(average(Q),dims,dims))
        err_method = zeros(n)
        for i=1:n
            xgtᵢ = vec(@views X_gt[i].P)
            xsolᵢ = vec(@views X_sol_iterative[i] * Q_avg)
            err_method[i] = error(xgtᵢ, xsolᵢ)
        end
        err = hcat(err, err_method)
    end
    return err

end

# Err = create_synthetic(0.0, error=angular_distance, outliers=0.1,  frames=15);
# rad2deg.(mean.(eachcol(Err)))