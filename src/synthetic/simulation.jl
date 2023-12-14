function create_synthetic(σ;kwargs...)
    normalize_matrix = get(kwargs, :normalize_matrix, false)
    dim = get(kwargs, :dimension, 4)
    n = get(kwargs, :frames, 25)
    ρ = get(kwargs, :holes_density, 0.5) 

    X_gt = SizedVector{n, Projectivity}(repeat([Projectivity(false)],n)) # Ground truth nodes with dimxdim matrices
    for i=1:n
        Xᵢ = rand(dim,dim)
        X_gt[i] = Projectivity(SMatrix{dim,dim,Float64}(Xᵢ/norm(Xᵢ)))
    end

    Z = SizedMatrix{n,n,Projectivity}(repeat([Projectivity(false)],n,n)) # Relative projectivities
    for i=1:n
        for j=1:i
            Z_ij = X_gt[i]*inv(X_gt[j]) + Projectivity(rand(Distributions.Normal(0, σ), dim, dim)) # Add noise
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
    
    if normalize_matrix
        Z = unit_normalize.(Z)
    end
    X_sol = projectivity_synch_spectral(copy(Z))
    # X_sol = iterative_projective_synchronization(Z; kwargs...)
    return X_gt, X_sol
    # X_solᵢ*Q = Xᵢ
    # Either take Q = Xᵢ for the anchor node i, OR
    # Take avg of Qᵢ = inv(X_solᵢ)*Xᵢ

end

# X_gt, X_ans = create_synthetic(0.0, frames=15);

# Q =  inv(unit_normalize(X_ans[1]))*X_gt[1]
# Q2 = inv(unit_normalize(X_ans[2]))*X_gt[2]
# 
# 
# unit_normalize(Q).P
# unit_normalize(Q2).P
# 