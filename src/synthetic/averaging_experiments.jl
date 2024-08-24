function create_environment(σ;kwargs...)
    num_vecs = get(kwargs, :num_vecs, 10)
    dim = get(kwargs, :dimension, 16)
    normalize = get(kwargs, :normalize, true)
    Ρ = get(kwargs, :outliers, 0.0)
    gt_vec = rand(dim)
    if normalize
        unit_normalize!(gt_vec)
    end

    M = zeros(dim, num_vecs)
    outlier_inds = StatsBase.sample(collect(1:num_vecs), Int(round(num_vecs*Ρ)), replace=false)
    for i=1:num_vecs
        θ = abs(rand(Distributions.Normal(0, σ)))
        M[:,i] = i in outlier_inds ? rand(dim) : rotate_vector(gt_vec, θ)
    end
    return gt_vec,M
end

function get_errors(σ::Float64, average_methods::AbstractVector{String}; error_measure=angular_distance, kwargs...)
    gt_vec, M = create_environment(σ; kwargs... )
    errs = []
    for method in average_methods
        if !occursin("robust", lowercase(method))
            if (occursin("sphere", method) || occursin("lagrangian", method)) && !occursin("A1", method)
                est_lagrangian = spherical_mean(M)
                errs = [errs; rad2deg(error_measure(est_lagrangian, gt_vec)) ]
            elseif occursin("A1", method) || occursin("Buss", method)
                if occursin("l1", lowercase(method))
                    est_A1_l1 = weighted_spherical_mean_A1_L1(M, SVector{size(M,2), Float64}(ones(size(M,2))))
                    errs = [errs; rad2deg(error_measure(est_A1_l1, gt_vec)) ]
                elseif occursin("l2", lowercase(method))
                    est_A1_l2 = weighted_spherical_mean_A1_L2(M, SVector{size(M,2), Float64}(ones(size(M,2))))
                    errs = [errs; rad2deg(error_measure(est_A1_l2, gt_vec)) ]
                end 
            elseif occursin("weiszfeld", lowercase(method))
                est_weiszfeld = weiszfeld(M)
                errs = [errs; rad2deg(error_measure(est_weiszfeld, gt_vec)) ]
            elseif occursin("euclidean", method)
                est_euclidean = ls_euclidean(M)
                errs = [errs; rad2deg(error_measure(est_euclidean, gt_vec)) ]
            elseif occursin("dyad", method)
                est_dyadic = average_dyadic(M)
                errs = [errs; rad2deg(error_measure(est_dyadic, gt_vec)) ]
            end
        else
            if (occursin("sphere", method) || occursin("lagrangian", method)) && !occursin("A1", method)
                est_lagrangian = iteratively_weighted_averaging(M, spherical_mean)
                errs = [errs; rad2deg(error_measure(est_lagrangian, gt_vec)) ]
            elseif occursin("A1", method) || occursin("Buss", method)
                if occursin("l1", lowercase(method))
                    est_A1_l1 = iteratively_weighted_averaging(M, weighted_spherical_mean_A1_L1)
                    errs = [errs; rad2deg(error_measure(est_A1_l1, gt_vec)) ]
                elseif occursin("l2", lowercase(method))
                    est_A1_l2 = iteratively_weighted_averaging(M, weighted_spherical_mean_A1_L2)
                    errs = [errs; rad2deg(error_measure(est_A1_l2, gt_vec)) ]
                end 
            elseif occursin("weiszfeld", lowercase(method))
                est_weiszfeld = iteratively_weighted_averaging(M,weiszfeld)
                errs = [errs; rad2deg(error_measure(est_weiszfeld, gt_vec)) ]
            elseif occursin("euclidean", method)
                est_euclidean = iteratively_weighted_averaging(M,ls_euclidean)
                errs = [errs; rad2deg(error_measure(est_euclidean, gt_vec)) ]
            elseif occursin("dyad", method)
                est_dyadic = iteratively_weighted_averaging(M,average_dyadic)
                errs = [errs; rad2deg(error_measure(est_dyadic, gt_vec)) ]
            end
        end
    end
    return errs
end

function averaging_sensitivity(avg_methods;fixed_σ=0.1, varying_parameter="noise", num_trials=100, param_min=0.0, param_step=0.05, param_max=0.2, kwargs...)
    E = Vector{Vector{Vector{Float64}}}(undef, length(collect(param_min:param_step:param_max)))
    for (i,param) in enumerate(collect(param_min:param_step:param_max))
        Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
        if occursin("noise", lowercase(varying_parameter))
            for j in tqdm(1:num_trials)
                Eᵢ[j] = get_errors(param, avg_methods; kwargs...)
            end
            E[i] = Eᵢ   
        elseif occursin("outliers", lowercase(varying_parameter))
            for j in tqdm(1:num_trials)
                Eᵢ[j] = get_errors(fixed_σ, avg_methods; outliers=param, kwargs...)
            end
            E[i] = Eᵢ   
        end
    end
    return E
end

# avg_methods = ["sphere", "A1_L1", "A1_L2", "weiszfeld", "robust-sphere", "robust-A1L1", "robust-A1L2", "robust-weiszfeld" ];
# Errs = averaging_sensitivity(avg_methods, num_vecs=25, fixed_σ=0.1, varying_parameter="outliers", param_step=0.05, param_max=0.6, num_trials=50);
# import GLMakie
# colors = ["red","green","blue", "orange", "red", "green", "blue", "orange"];
# f,a = plot_sensitivity_curves(Errs, varying_parameter="Outliers", avg_methods, colors=colors, parameter= rad2deg.(collect(0:0.05:0.6)))
# GLMakie.save("pointOneNoise_Outliers_Frames25.png", f)
