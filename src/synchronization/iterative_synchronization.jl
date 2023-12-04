function iterative_projective_synchronization(Z::AbstractMatrix{Projectivity};kwargs...) 
    #Z contains nxn relative projective quantities (4x4 matrices) for n nodes
    
    method = get(kwargs, :averaging_method, "least_squares")
    dim = get(kwargs, :dimension, 4)
    max_it = get(kwargs, :max_iterations, 1000)
    max_updates = get(kwargs, :max_updates, 70)
    min_updates = get(kwargs, :min_updates, 20)
    δ = get(kwargs, :δ, 1e-12)
    
    if occursin("euclidean", lowercase(method))
        average = average_ls_euclidean
    elseif occursin("crossproduct", lowercase(method))
        average = average_ls_crossProduct
    elseif occursin("dyad", lowercase(method))
        average = average_dyadic
    elseif occursin("sphere", lowercase(method)) || occursin("spherical", lowercase(method))
        average = average_sphere
    elseif occursin("gaia", lowercase(method))
        average = average_lsg
    elseif occursin("weiszfeld", lowercase(method))
        average = average_weiszfeld
    elseif occursin("orthogonal", lowercase(method))
        average = average_ls_orthogonal
    end

    n = size(Z,1)
    Adj = MMatrix{n,n, Bool}(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    
    G = Graph(Adj)
    steady = zeros(Bool,n )
    updated = zeros(Int, n)

    # Initialize nodes to have identity dimxdim matrices (this will be overridden for all but the anchor)
    X = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dim,dim, Float64}(I))], n)) 
    # Set anchor node according to maximum closeness centrality
    C = closeness_centrality(G)
    _,anchor = findmax(C)
    updated[anchor] = 1
    
    iter=0
    while iter < max_it
        iter += 1
        j = rand(1:n)
        if updated[j]  > max_updates || steady[j]
            continue
        end
        
        N = neighbors(G, j)
        if iszero(updated[N])
            continue
        end

        oldX = X[j]
        X[j] = Projectivity(average(j, N[findall(!zero, updated[N])], Z, X))
        updated[j] += 1
        steady[j] = updated[j] > min_updates && norm(oldX.P - X[j].P)/norm(oldX.P) < δ
    
        if all(steady) || max_it || all(updated .> max_updates)
            break
        end

    end
    
    return X
end   

function average_ls_euclidean(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(Z[j,i].P*X[i].P)
        M = [M h]
    end
    M = M[:,2:end]

    avg = reshape(ls_euclidean(M), n,n)
    return SMatrix{n,n,Float64}(avg)
end


function average_ls_crossProduct(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(Z[j,i].P*X[i].P)
        M = [M h]
    end
    M = M[:,2:end]

    avg = reshape(ls_crossProduct(M), n,n)
    return SMatrix{n,n,Float64}(avg)
end

function average_sphere(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(Z[j,i].P*X[i].P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    a = M[:,1]
    for i=2:size(M,2)
        if dot(a,M[:,i]) < 0
            M[:,i] = -M[:,i]
        end
    end

    avg = reshape(spherical_mean(M),n,n)
    return SMatrix{n,n,Float64}(avg)
end

function average_weiszfeld(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(Z[j,i].P*X[i].P)
        M = [M h]
    end
    M = M[:,2:end]
    
    avg = reshape(weiszfeld(M),n,n)
    return SMatrix{n,n,Float64}(avg)
end

function average_ls_orthogonal(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(X[i].P)
        M = [M (SMatrix{n^2,n^2, Float64}(I) - dot(h,h)) * kron(SMatrix{n,n,Float64}(I), Z[i,j].P )]
    end
    M = M[:,2:end]

    S = svd(M)
    return SMatrix{n,n,Float64}( reshape(S.V[:,end],n,n) )
end

function average_lsg(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for i in N
        h = vec(X[i].P)
        M = [M (dot(h,h)*SMatrix{n^2,n^2, Float64}(I) - h*h' ) * kron(SMatrix{n,n,Float64}(I), Z[i,j].P) ]
    end
    M = M[:,2:end]

    S = svd(M)
    return SMatrix{n,n,Float64}( reshape(S.V[:,end],n,n) )
end

function average_dyadic(j::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2, n^2)

    for i in N
        h = vec(Z[j,i].P*X[i].P)
        M = M + (h*h')/dot(h,h)
    end
    
    return SMatrix{n,n,Float64}(reshape(eigvecs(M)[:,end],n,n))
end
