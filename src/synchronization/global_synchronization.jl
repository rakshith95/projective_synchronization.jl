function projectivity_synch_spectral(Z::AbstractMatrix{Projectivity})
    # Based on spectral solution for synchronization in GL(d) : https://iris.unitn.it/retrieve/handle/11572/272068/360005/ijcv20.pdf

    n = size(Z,1)
    dims = size(Z[findfirst(!zero, Z)].P, 1)    
    X = SizedVector{n, Projectivity}(repeat([Projectivity(false )], n)) 

    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    for i=1:n
        for j=1:i
            Z[i,j] = unit_normalize(Z[i,j], "determinant")
            Z[j,i] = inv(Z[i,j])
        end
    end
    
    iD = diagm(1 ./sum.(eachcol(Adj)))
    Z′ = sparse(zeros(Complex, n*dims, n*dims))
    unwrap!(Z′, Z)
    S = kron(iD, SMatrix{dims,dims,ComplexF64}(I))*Z′
    # M, will contain the blocks of the nodes 
    λ, M = eigs(S, nev=4)
    M = M[:, end:-1:1] # Reverse order for highest to lowest
    for k=1:n
        Xₖ = M[(k-1)*dims + 1:(k-1)*dims+dims,:]
        #Force last block to I 
        Xₖ = Xₖ/M[(n-1)*dims + 1:(n-1)*dims+dims,:]
        c = real(view(Xₖ,:)) \ imag(view(Xₖ,:))
        scale = exp(-im*atan(c))
        X[k] = Projectivity(real(scale.*Xₖ))
    end

    return X
end

function projectivity_synch_spectral(Z::AbstractMatrix{Projectivity}, weights::SMatrix{N,N,T}) where {N,T<:AbstractFloat}
    # Based on spectral solution for synchronization in GL(d) : https://iris.unitn.it/retrieve/handle/11572/272068/360005/ijcv20.pdf

    n = size(Z,1)
    dims = size(Z[findfirst(!zero, Z)].P, 1)    
    X = SizedVector{n, Projectivity}(repeat([Projectivity(false )], n)) 

    # Adj[findall(!zero, Z)] .= 1
    # Adj = sparse(Adj .* weights)
    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] = weights[findall(!zero,Z)]
    
    for i=1:n
        for j=1:i
            Z[i,j] = unit_normalize(Z[i,j], "determinant")
            Z[j,i] = inv(Z[i,j])
        end
    end
    
    iD = diagm(1 ./sum.(eachcol(Adj)))
    Z′ = sparse(zeros(Complex, n*dims, n*dims))
    unwrap!(Z′, Z)
    S = kron(iD, SMatrix{dims,dims,ComplexF64}(I))*Z′
    # M, will contain the blocks of the nodes 
    λ, M = eigs(S, nev=4)
    M = M[:, end:-1:1] # Reverse order for highest to lowest
    for k=1:n
        Xₖ = M[(k-1)*dims + 1:(k-1)*dims+dims,:]
        #Force last block to I 
        Xₖ = Xₖ/M[(n-1)*dims + 1:(n-1)*dims+dims,:]
        c = real(view(Xₖ,:)) \ imag(view(Xₖ,:))
        scale = exp(-im*atan(c))
        X[k] = Projectivity(real(scale.*Xₖ))
    end

    return X
end

function spanning_tree_synchronization(Z::AbstractMatrix{Projectivity})
    n = size(Z,1)
    dims = size(Z[findfirst(!zero, Z)].P, 1)    
    X = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dims,dims,Float64}(I) )], n)) 

    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    G = Graph(Adj)
    
    # Get spanning tree
    ST = prim_mst(G)
    root_node = setdiff(collect(1:n), [v.dst for v in ST])[1]
    
    for e in ST
        src = e.src
        dst = e.dst
        if X[dst].P != SMatrix{dims,dims,Float64}(I)
            continue
        end
        
        #loop
        path = []
        while src != root_node && X[src].P == SMatrix{dims,dims,Float64}(I)
            append!(path, src)
            src = ST[src-1].src
        end
        for p in reverse(path)
            X[p] = inv(Z[ ST[p-1].src , p])*X[ST[p-1].src] 
        end
        
        X[dst] = inv(Z[e.src , dst])*X[e.src] 
    end 

    return X
end

function spanning_tree_synchronization(Z::AbstractMatrix{Projectivity}, weights::SMatrix{N,N,T}) where {N,T<:AbstractFloat}
    n = size(Z,1)
    dims = size(Z[findfirst(!zero, Z)].P, 1)    
    X = SizedVector{n, Projectivity}(repeat([Projectivity(SMatrix{dims,dims,Float64}(I) )], n)) 

    Adj = sparse(zeros(n,n))
    Adj[findall(!zero, Z)] .= 1
    G = Graph(Adj)
    
    # Get spanning tree
    ST = prim_mst(G, weights)
    root_node = setdiff(collect(1:n), [v.dst for v in ST])[1]
    
    for e in ST
        src = e.src
        dst = e.dst
        if X[dst].P != SMatrix{dims,dims,Float64}(I)
            continue
        end
        
        #loop
        path = []
        while src != root_node && X[src].P == SMatrix{dims,dims,Float64}(I)
            append!(path, src)
            src = ST[src-1].src
        end
        for p in reverse(path)
            X[p] = inv(Z[ ST[p-1].src , p])*X[ST[p-1].src] 
        end
        
        X[dst] = inv(Z[e.src , dst])*X[e.src] 
    end 

    return X
end