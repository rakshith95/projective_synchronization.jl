function angular_noise!(x::AbstractVector{T}, θ::T) where T
    dim = length(x)
    B = missing
    while true
        B = rand(dim,dim)
        if !iszero(det(B))
            break
        end
    end

    # Take the bases into tangent space
    for i=1:dim
        B[:,i] = sphere_to_tangent(x, view(B,:,i))
    end

    # Get a random vector as linear combination of bases
    coords = rand(dim)
    v = vec(sum(coords'.*B, dims=2))
    unit_normalize!(v)
    v = θ*v

    # Take vector from tangent space to sphere
    RotateUnitInDirection!(x, v)
end

function angular_noise(x::AbstractVector{T}, θ::T) where T
    x_cpy = copy(x)
    angular_noise!(x_cpy, θ)
    return x_cpy
end

function angular_noise!(Z::Projectivity, θ::T) where T<:AbstractFloat
    if !isapprox(norm(Z.P), 1.0)
        unit_normalize!(Z)
    end
    angular_noise!( vec(Z.P) , θ) 
end

function compute_Q(Q, X_sol, X_gt, average_method)
    n = length(X_gt)
    dims = size(X_gt[1].P,1)
    for i=1:n
        Q[:,i] = vec((inv(X_sol[i])*X_gt[i]).P)
    end
    for j=2:size(Q,2)
        if dot(view(Q,:,1), view(Q,:,j)) < 0
            Q[:,j] = -Q[:,j]
        end
    end
    return SMatrix{dims,dims,Float64}(reshape(average_method(Q),dims,dims))
end
    
function compute_err(X_gt, X_sol, Q_avg, error_measure)
    n = length(X_gt)
    err = zeros(n)
    for i=1:n
        xgtᵢ = vec(@views X_gt[i].P)
        xsolᵢ = vec(@views X_sol[i] * Q_avg)
        err[i] = error_measure(xgtᵢ, xsolᵢ)
    end
    return err
end

function compute_Z!(Z::AbstractArray{Projectivity}, X::AbstractVector{Projectivity};dims=4, noise_type="elemental", σ=0.0)
    n = size(Z,1)
    for i=1:n
        for j=1:i
            if i==j
                continue
            end
            if occursin("elemental", noise_type)
                Z[i,j] = X[i]*inv(X[j]) + Projectivity(rand(Distributions.Normal(0, σ), dims, dims)) # Add noise
            elseif occursin("angular", noise_type)
                Z[i,j] = unit_normalize(X[i]*inv(X[j]))
                θ = abs(rand(Distributions.Normal(0, σ)))
                angular_noise!(Z[i,j], θ)
            end
            Z[j,i] = inv(Z[i,j]) # Symmetric block is inverse
        end
    end
end


function create_synthetic(σ;noise_type="elemental_gaussian", error=orthogonal_projection_distance,average=spherical_mean, averaging_methods=["sphere"], kwargs...)
    normalize_matrix = get(kwargs, :normalize_matrix, false)
    dims = get(kwargs, :dimension, 4)
    n = get(kwargs, :frames, 25)
    ρ = get(kwargs, :holes_density, 0.5) 
    Ρ = get(kwargs, :outliers, 0.0)
    
    X_gt = SizedVector{n, Projectivity}(repeat([Projectivity(false)],n)) # Ground truth nodes with dimxdim matrices
    for i=1:n
        Xᵢ = rand(dims,dims)
        X_gt[i] = Projectivity(MMatrix{dims,dims,Float64}(Xᵢ/norm(Xᵢ)))
    end

    # Z = SizedMatrix{n,n,Projectivity}(repeat([Projectivity(false)],n,n)) # Relative projectivities
    Z = SparseMatrixCSC{Projectivity, Integer}(repeat([Projectivity(false)],n,n)) # Relative projectivities
    compute_Z!(Z, X_gt;dims=dims, noise_type=noise_type, σ=σ)
    
    # Make holes
    A = missing
    while true 
        A = sprand(n, n, ρ)
        A[A.!=0] .= 1.0
        A = sparse(ones(n,n)) - A
        A = triu(A,1) + triu(A,1)'
        G = Graph(A)
        if Graphs.is_connected(G)
            break
        end
    end
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
    # return X_gt, Z
    #Global Methods:
    # X_solᵢ*Q = Xᵢ
    # Either take Q = Xᵢ for the anchor node i, OR
    # Take avg of Qᵢ = inv(X_solᵢ)*Xᵢ ∀ i ∈ 1..n
    
    Q = MMatrix{dims*dims, n, Float64}(zeros(dims*dims, n))
    #1. Spanning Tree 
    X_sol_spanningTree = spanning_tree_synchronization(copy(Z))
    Q_avg_spanningTree = compute_Q(Q, X_sol_spanningTree, X_gt, average)
    err = compute_err(X_gt, X_sol_spanningTree, Q_avg_spanningTree, error)

    #2. Spectral 
    # X_sol_spectral = projectivity_synch_spectral(copy(Z))
    # Q_avg_spectral = compute_Q(Q, X_sol_spectral, X_gt, average)
    # err = hcat(err, compute_err(X_gt, X_sol_spectral, Q_avg_spectral, error))
 
    #3. Robust Spectral
    # X_sol_robust_spectral = iteratively_weighted_synchronization(copy(Z), "spectral", max_it=30, δ=1e-2)
    # Q_avg_robust_spectral = compute_Q(Q, X_sol_robust_spectral, X_gt, average)
    # err = hcat(err, compute_err(X_gt, X_sol_robust_spectral, Q_avg_robust_spectral, error))

    # Iterative methods
    for method in averaging_methods
        if occursin("irls", lowercase(method))
            X_sol_iterative = iteratively_weighted_synchronization(copy(Z), method, error_measure=error, max_it=100, averaging_max_it=5, δ=deg2rad(0.1))
        else
            X_sol_iterative = iterative_projective_synchronization(copy(Z);X₀=X_sol_spanningTree, averaging_method=method,kwargs...)
        end            
        Q_avg_method = compute_Q(Q, X_sol_iterative, X_gt, average)
        err = hcat(err, compute_err(X_gt, X_sol_iterative, Q_avg_method, error))
    end

    return err
    
end

# avg_methods = ["sphere", "sphere-robust", "sphere-irls" ];
# avg_methods = ["sphere", "weiszfeld", "robust-sphere", "robust-A1L2", "robust-weiszfeld" ]
# Err = create_synthetic(0.1, average=spherical_mean , averaging_methods=avg_methods, error=angular_distance, outliers=0.0,  frames=25);
# rad2deg.(mean.(eachcol(Err)))

# X_gt, Z = create_synthetic(0.1, outliers=0.3);