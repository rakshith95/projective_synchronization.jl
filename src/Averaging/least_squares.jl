function ls_euclidean(n::Int, X::AbstractVector)
    # least squares minimization for hᵢ ~= c
    #X is stacked vector of some m=length(X)/n unit vectors of dimension n each
    A = kron(ones(div(length(X), n)), SMatrix{n,n,eltype(X)}(I))
    return A\X
end

function make_crossProduct_matrix(a::AbstractVector)
    n = length(a)
    B = MMatrix{div(n*(n-1),2), n, eltype(a)}(zeros(div(n*(n-1),2), n))
    x = 1
    r₀ = 0
    for i=1:n-1
        r = n-i
        Bᵢ = MMatrix{r,n,eltype(a)}(zeros(r,n))
        for j=1:r
            Bᵢ[j,i] = -a[i+j]
        end
        Bᵢ[:,i+1:end] = a[i]*MMatrix{r,r,eltype(a)}(I)
        B[x+r₀:x+r₀+r-1,:] = Bᵢ
        x = x+r₀
        r₀ = r
    end
    return B
end

function ls_crossProduct(M::AbstractMatrix)
    #Each column of M is a vectorized Projetive Matrix
    # From "Uncalibrated Dynamic Stereo Using Parallax" by Francesco Malapelle et al. : https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6703743
    hₓᵢ = make_crossProduct_matrix.(eachcol(M))
    hₓ = reduce(vcat, hₓᵢ)
    U_Σ_vt = svd(hₓ)
    V = U_Σ_vt.V
    #= Last column of V in SVD, i.e. the eigenvector corresponding to the minimum singular value minimizes ||hₓc||₂ for ||c||=1
    https://math.berkeley.edu/~hutching/teach/54-2017/svd-notes.pdf shows proof for opposite case, i.e. max ||Ax|| for ||x||=1  =#

    return V[:,end]
end

# c = ls_crossProduct(A)