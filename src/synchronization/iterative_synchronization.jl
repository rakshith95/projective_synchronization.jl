function iterative_projective_synchronization(Z::AbstractMatrix{Projectivity};kwargs...) 
    #Z contains nxn relative projective quantities (4x4 matrices) for n nodes
    
    dim = get(kwargs, :dimension, 4)
    # if occursin("least-squares", method)
    n = size(Z,1)
    Adj = zeros(n,n)
    Adj[findall(!zero, Z)] .= 1
    
    G = Graph(Adj)
    steady = zeros(Bool,n )
    updated = zeros(n)

    

end   

# R = SMatrix{3,3,Float64}(rand(3,3))
# Z = SMatrix{2,2, Projectivity}([ Projectivity(R) Projectivity(false); Projectivity(false) Projectivity(4*R)])


# iterative_projective_synchronization(Z)


# function average_ls_skew(H::AbstractMatrix)
