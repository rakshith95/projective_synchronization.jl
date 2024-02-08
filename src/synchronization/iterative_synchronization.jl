function iterative_projective_synchronization(Z::AbstractMatrix{Projectivity};X₀=nothing,  weights=nothing, kwargs...) 
    #Z contains nxn relative projective quantities (4x4 matrices) for n nodes
    
    method = get(kwargs, :averaging_method, "sphere")
    dim = get(kwargs, :dimension, 4)
    max_it = get(kwargs, :max_iterations, 100)
    max_updates = get(kwargs, :max_updates, 100)
    min_updates = get(kwargs, :min_updates, 20)
    δ = get(kwargs, :δ, 1e-10)
    
    if occursin("robust", lowercase(method))
        if occursin("sphere", lowercase(method))
            average = robust_average_sphere
        elseif occursin("A1", uppercase(method)) 
            average = robust_average_sphere_A1
        elseif occursin("weiszfeld", lowercase(method))
            average = robust_average_weiszfeld
        end
    else
        if occursin("euclidean", lowercase(method))
            average = average_ls_euclidean
        elseif occursin("crossproduct", lowercase(method))
            average = average_ls_crossProduct
        elseif occursin("dyad", lowercase(method))
            average = average_dyadic
        elseif occursin("sphere", lowercase(method)) || occursin("spherical", lowercase(method))
            average = average_sphere
        elseif occursin("A1", uppercase(method)) || occursin("buss", lowercase(method))
            average = average_sphere_A1
        elseif occursin("gaia", lowercase(method))
            average = average_lsg
        elseif occursin("weiszfeld", lowercase(method))
            average = average_weiszfeld
        elseif occursin("orthogonal", lowercase(method))
            average = average_ls_orthogonal
        end
    end
    n = size(Z,1)
    if isnothing(weights)
        weights = ones(n, n)
    end
    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    
    G = Graph(Adj)
    steady = zeros(Bool,n )
    updated = zeros(Int, n)

    # Initialize nodes to have identity dimxdim matrices (this will be overridden for all but the anchor)
    if isnothing(X₀)
        X = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dim,dim, Float64}(I))], n)) 
    else
        X = copy(X₀)
        updated .= 1
    end
    # Set anchor node according to maximum closeness centrality
    C = closeness_centrality(G)
    _,anchor = findmax(C)
    updated[anchor] = max_updates
    steady[anchor] = true
    iter=0

    nodes = collect(1:n)
    # nodes = sortperm(C)
    # Try reverse
    
    exit_loop = false
    while !exit_loop
        if iter>=max_it
            exit_loop=true
            continue
        end
        Random.shuffle!(nodes)
        for i in nodes
            # wts =  (max_updates .- updated)./max_updates .* .!steady
            # Change to random permutation at each step. i.e. update each node once before updating any node again.
            # i = StatsBase.wsample(unit_normalize!(wts))
            if steady[i] 
                continue
            end
            N = neighbors(G, i)
            if iszero(updated[N])
                continue
            end

            oldX = X[i]
            updated_N = N[updated[N].!=0]
            X[i] = Projectivity(average(i, X, updated_N, Z, SVector{length(updated_N),Float64}(weights[i,updated_N])))
            updated[i] += 1
            steady[i] = updated[i] > min_updates && norm((oldX - X[i]).P)/norm(oldX.P) < δ
            if all(steady)  || all(updated .>= max_updates)
                exit_loop=true
                break
            end
        end
        iter += 1
    end
    println(iter)
    return X
end   

function average_ls_euclidean(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]
    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end

    return SMatrix{n,n,Float64}(reshape(ls_euclidean(M, weights), n,n))
end


function average_ls_crossProduct(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    return SMatrix{n,n,Float64}(reshape(ls_crossProduct(M), n,n))
end

function average_sphere(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(spherical_mean(M, weights),n,n))
end

function robust_average_sphere(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights=nothing) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(iteratively_weighted_averaging(M, spherical_mean),n,n))
end


function average_sphere_A1(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(weighted_spherical_mean_A1(M, weights),n,n))    
end

function robust_average_sphere_A1(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity},weights=nothing) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(iteratively_weighted_averaging(M, weighted_spherical_mean_A1_L1),n,n))    
end

function average_weiszfeld(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end

    return SMatrix{n,n,Float64}(reshape(weiszfeld(M, weights),n,n))
end

function robust_average_weiszfeld(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity},weights=nothing) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    for j=2:size(M,2)
        if dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end

    return SMatrix{n,n,Float64}(reshape(iteratively_weighted_averaging(M, weiszfeld),n,n))
end

function average_ls_orthogonal(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2, n^2)
    for j in N
        h = vec(X[j].P)
        M = vcat(M,(SMatrix{n^2,n^2, Float64}(I) - (h*h')/(h'*h) ) * kron(SMatrix{n,n,Float64}(I), Z[j,i].P ))
    end
    M = M[n^2+1:end,:]
    S = svd(M)
    return SMatrix{n,n,Float64}( reshape(S.V[:,end],n,n) )
end

function average_lsg(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2, n^2)
    for j in N
        h = vec(unit_normalize(X[j]).P)
        M = vcat(M,((h'*h)*SMatrix{n^2,n^2, Float64}(I) - h*h' ) * kron(SMatrix{n,n,Float64}(I), Z[j,i].P) )
    end
    M = M[n^2+1:end,:]
    S = svd(M)
    return SMatrix{n,n,Float64}( reshape(S.V[:,end],n,n) )
end

function average_dyadic(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}, weights::AbstractVector{Float64}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2, n^2)

    for (ct,j) in enumerate(N)
        h = vec((Z[i,j]*unit_normalize(X[j])).P)
        M = M + weights[ct]*(h*h')/(h'*h)
    end
    ev_max = eigvecs(M)[:,end]
    return SMatrix{n,n,Float64}(reshape(ev_max,n,n))
end

function average_dyadic(M::AbstractArray)
    n = size(M,1)
    D = zeros(n, n)
    for i in 1:size(M,2)
        D = D + (M[:,i]*M[:,i]')/(M[:,i]'*M[:,i])
    end
    ev_max = eigvecs(D)[:,end]
    return SVector{n, Float64}(ev_max)
end

function average_dyadic(M::AbstractArray, weights::AbstractVector{T}) where T<:AbstractFloat
    if !isapprox(sum(weights),1.0)
        weights = weights/sum(weights)
    end
    n = size(M,1)
    D = zeros(n, n)
    for i in 1:size(M,2)
        D = D + weights[i]*(M[:,i]*M[:,i]')/(M[:,i]'*M[:,i])
    end
    ev_max = eigvecs(D)[:,end]
    return SVector{n, Float64}(ev_max)
end