function iterative_projective_synchronization(Z::AbstractMatrix{Projectivity};X₀=nothing,  kwargs...) 
    #Z contains nxn relative projective quantities (4x4 matrices) for n nodes
    
    method = get(kwargs, :averaging_method, "sphere")
    dim = get(kwargs, :dimension, 4)
    max_it = get(kwargs, :max_iterations, 1000)
    max_updates = get(kwargs, :max_updates, 70)
    min_updates = get(kwargs, :min_updates, 20)
    δ = get(kwargs, :δ, 1e-10)
    
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

    n = size(Z,1)
    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    
    G = Graph(Adj)
    steady = zeros(Bool,n )
    updated = zeros(Int, n)

    # Initialize nodes to have identity dimxdim matrices (this will be overridden for all but the anchor)
    if isnothing(X₀)
        X = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dim,dim, Float64}(I))], n)) 
    else
        X = X₀
    end
    # Set anchor node according to maximum closeness centrality
    C = closeness_centrality(G)
    _,anchor = findmax(C)
    updated[anchor] = 1
    steady[anchor] = true
    iter=0

    nodes = collect(1:n)
    exit_loop = false
    while !exit_loop
        Random.shuffle!(nodes)
        for i in nodes
            iter += 1
            # wts =  (max_updates .- updated)./max_updates .* .!steady
            # Change to random permutation at each step. i.e. update each node once before updating any node again.
            # i = StatsBase.wsample(unit_normalize!(wts))
            if i==anchor || steady[i] 
                continue
            end
            N = neighbors(G, i)
            if iszero(updated[N])
                continue
            end

            oldX = X[i]
            X[i] = Projectivity(average(i, X, N[updated[N].!=0], Z))
            updated[i] += 1
            steady[i] = updated[i] > min_updates && norm((oldX - X[i]).P)/norm(oldX.P) < δ

            if all(steady)  || all(updated .> max_updates)
                exit_loop=true
                break
            end
        end
    end
    return X
end   

function average_ls_euclidean(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]
    for j=2:size(M,2)
        if dot(M[:,1],M[:,j]) < 0
            M[:,j] = -M[:,j]
        end
    end

    return SMatrix{n,n,Float64}(reshape(ls_euclidean(M), n,n))
end


function average_ls_crossProduct(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    return SMatrix{n,n,Float64}(reshape(ls_crossProduct(M), n,n))
end

function average_sphere(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(M[:,1],M[:,j]) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(spherical_mean(M),n,n))
end

function average_sphere_A1(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}; algorithm="lagrangian") where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]

    #identity antipodal points
    for j=2:size(M,2)
        if dot(M[:,1],M[:,j]) < 0
            M[:,j] = -M[:,j]
        end
    end
    return SMatrix{n,n,Float64}(reshape(weighted_spherical_mean_A1(M),n,n))    
end


function average_weiszfeld(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2)
    for j in N
        h = vec((Z[i,j]*X[j]).P)
        M = [M h]
    end
    M = M[:,2:end]
    for j=2:size(M,2)
        if dot(M[:,1],M[:,j]) < 0
            M[:,j] = -M[:,j]
        end
    end

    return SMatrix{n,n,Float64}(reshape(weiszfeld(M),n,n))
end

function average_ls_orthogonal(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
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

function average_lsg(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
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

function average_dyadic(i::Int, X::AbstractVector{Projectivity}, N::AbstractVector{T}, Z::AbstractMatrix{Projectivity}) where T
    n = size(X[1].P)[1]
    M = zeros(n^2, n^2)

    for j in N
        h = vec((Z[i,j]*unit_normalize(X[j])).P)
        M = M + (h*h')/(h'*h)
    end
    ev_max = eigvecs(M)[:,end]
    return SMatrix{n,n,Float64}(reshape(ev_max,n,n))
end

function average_dyadic(M::AbstractArray)
    n = size(M,1)
    D = zeros(n, n)
    for i in size(M,2)
        h = vec(M[:,i])
        D = D + (h*h')/(h'*h)
    end
    ev_max = eigvecs(D)[:,end]
    return SVector{n, Float64}(ev_max)
end