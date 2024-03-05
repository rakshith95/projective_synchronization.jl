push!(LOAD_PATH, "/home/rakshith/PoliMi/Projective Synchronization/projective-synchronization-julia/projective_synchronization.jl")

import projective_synchronization
import MAT
import LinearAlgebra

function average_from_mat(filename::String, output_filename::String; method="robust-sphere", varname="L")
    file = MAT.matopen(filename)
    M = read(file, varname)
    close(file)
    #identity antipodal points
    for j=2:size(M,2)
        if LinearAlgebra.dot(view(M,:,1),view(M,:,j)) < 0
            M[:,j] = -M[:,j]
        end
    end
    avg=missing;
    if (occursin("sphere", method) || occursin("lagrangian", method)) 
        if occursin("robust", lowercase(method))
            avg = projective_synchronization.iteratively_weighted_averaging(M, projective_synchronization.spherical_mean)
        else
            avg = projective_synchronization.spherical_mean(M)    
        end
    elseif occursin("weiszfeld", lowercase(method))
        if occursin("robust", lowercase(method))
            avg = projective_synchronization.iteratively_weighted_averaging(M, projective_synchronization.weiszfeld)
        else
            avg = projective_synchronization.weiszfeld(M)
        end
    elseif occursin("euclidean", method)
        avg = projective_synchronization.ls_euclidean(M)
    elseif occursin("dyad", method)
        avg = projective_synchronization.average_dyadic(M)
    end
    file = MAT.matopen(output_filename, "w")
    write(file, "avg", Vector(avg))
    close(file)
end


function main(args)
    L_file = args[1]
    output_file = args[2]
    method = "sphere"
    if length(args) > 2
        method = args[3]
    end 
    average_from_mat(L_file, output_file, method=method);
end

main(ARGS)