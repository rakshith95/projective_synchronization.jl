import Pkg
# Pkg.add(path="/home/rakshith/PoliMi/Projective Synchronization/projective-synchronization-julia/projective_synchronization.jl/")
# Pkg.add("MAT")
import projective_synchronization, SparseArrays
import MAT

function run_sync_from_mat(filename::String, output_filename::String, dimension; method="sphere", varname="Z")
    file = MAT.matopen(filename)
    Z_mat = read(file, varname)
    close(file)
    n = div(size(Z_mat,1), dimension)
    Z = SparseArrays.SparseMatrixCSC{projective_synchronization.Projectivity, Integer}(repeat([projective_synchronization.Projectivity(false)],n,n)) # Relative projectivities
    println(Z[1])
    projective_synchronization.wrap!(Z, Z_mat, dimension)
    X = projective_synchronization.iterative_projective_synchronization(Z, averaging_method=method, δ=1e-6, max_iterations=200);
    X_vec = repeat([zeros(dimension, dimension)], n)
    for i=1:n
        tmp = zeros(dimension, dimension)
        unwrap!(X[i], tmp)
        X_vec[i] = tmp
    end
    file = MAT.matopen(output_filename, "w")
    write(file, "T", X_vec)
    close(file)
end


function main(args)
    Z_file = args[1]
    output_file = args[2]
    dimension = 4
    method = "sphere"
    if length(args) > 2
        dimension = args[3]
    end
    if length(args) > 3
        method = args[4]
    end 
    # output_file = "synch.mat"
    # Z_file = "/home/rakshith/PoliMi/Projective Synchronization/projective-synchronization-julia/GPSFM-code/GPSFM/Z.mat"
    run_sync_from_mat(Z_file, output_file, dimension, method=method);
end

main(ARGS)