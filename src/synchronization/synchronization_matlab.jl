push!(LOAD_PATH, "/home/rakshith/PoliMi/Projective Synchronization/projective-synchronization-julia/projective_synchronization.jl")
import projective_synchronization,SparseArrays, LinearAlgebra
import MAT

function run_sync_from_mat(Z_filename::String, weights_filename, output_filename::String, dimension; method="sphere", varname="Z")
    file = MAT.matopen(Z_filename)
    Z_mat = read(file, varname)
    close(file)
    weights_mat = nothing
    if !isnothing(weights_filename)
        weights_file = MAT.matopen(weights_filename)
        weights_mat = read(weights_file, "edgeWeights")
        close(weights_file)
    end
    n = div(size(Z_mat,1), dimension)
    Z = SparseArrays.SparseMatrixCSC{projective_synchronization.Projectivity, Integer}(repeat([projective_synchronization.Projectivity(false)],n,n)) # Relative projectivities
    projective_synchronization.wrap!(Z, Z_mat, dimension)
    X₀ = projective_synchronization.spanning_tree_synchronization(copy(Z))
    if occursin("irls", method)
        println("Running IRLS")
        X, wts = projective_synchronization.iteratively_weighted_synchronization(Z, method, weights=weights_mat, weight_function=projective_synchronization.cauchy, c=projective_synchronization.c_cauchy, max_it=20, averaging_max_it=15, averaging_max_it_init=75, δ_irls=deg2rad(0.1), δ=1e-12, anchor="centrality", update="start-centrality-update-all" );
        # X = projective_synchronization.spanning_tree_synchronization(copy(Z))
    else
        X = projective_synchronization.iterative_projective_synchronization(Z, averaging_method=method, weights=weights_mat, δ=1e-12, max_iterations=n*100, anchor="fixed", update="start-centrality-update-all");
    end
    X_vec = repeat([zeros(dimension, dimension)], n)
    for i=1:n
        tmp = zeros(dimension, dimension)
        projective_synchronization.unwrap!(X[i], tmp)
        X_vec[i] = tmp
    end

    #=
    Ẑ = SparseArrays.SparseMatrixCSC{projective_synchronization.Projectivity, Integer}(repeat([projective_synchronization.Projectivity(false )],n,n)) # Relative projectivities
    # compute_Z!(Ẑ, X_gt);
    projective_synchronization.compute_Z!(Ẑ, X)
    Z_new = zero(Z_mat)
    projective_synchronization.unwrap!(Z_new, Ẑ)

    file = MAT.matopen("Z_sync.mat", "w")
    write(file, "Z_sync", Z_new)
    close(file)
    
    wts = projective_synchronization.compute_weights(Z,Ẑ, dimension)
    file = MAT.matopen("wts.mat", "w")
    write(file, "wts", wts)
    close(file)
    =#
    file = MAT.matopen(output_filename, "w")
    write(file, "T", X_vec)
    close(file)
end


function main(args)
    Z_file = args[1]
    output_file = args[2]
    dimension = 4
    weights_file = nothing
    method = "sphere"
    if length(args) > 2
        method = args[3]
    end 
    if length(args) > 3
        weights_file = args[4]
    end
    if length(args) > 4
        dimension = args[5]
    end

    run_sync_from_mat(Z_file, weights_file, output_file, dimension, method=method);
end

main(ARGS)