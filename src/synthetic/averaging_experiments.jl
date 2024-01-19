function create_environment(σ;kwargs...)
    num_vecs = get(kwargs, :num_vecs, 10)
    dim = get(kwargs, :dimension, 16)
    normalize = get(kwargs, :normalize, true)

    gt_vec = rand(dim)
    if normalize
        unit_normalize!(gt_vec)
    end

    M = zeros(dim, num_vecs)
    for i=1:num_vecs
        θ = abs(rand(Distributions.Normal(0, σ)))
        M[:,i] = angular_noise(gt_vec, θ)
    end
    return gt_vec,M
end

function get_errors(σ::Float64, average_methods::AbstractVector{String}; error_measure=angular_distance, kwargs...)
    gt_vec, M = create_environment(σ; kwargs... )
    errs = []
    for method in average_methods
        if (occursin("sphere", method) || occursin("lagrangian", method)) && !occursin("A1", method)
            est_lagrangian = spherical_mean(M)
            errs = [errs; rad2deg(error_measure(est_lagrangian, gt_vec)) ]
        elseif occursin("weighted_sphere", method) || occursin("A1", method) || occursin("Buss", method)
            if occursin("l1", lowercase(method))
                est_A1_l1 = weighted_spherical_mean_A1(M, p=1)
                errs = [errs; rad2deg(error_measure(est_A1_l1, gt_vec)) ]
            elseif occursin("l2", lowercase(method))
                est_A1_l2 = weighted_spherical_mean_A1(M, p=2)
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
    end
    return errs
end

function averaging_sensitivity(avg_methods;num_trials=100, param_min=0.0, param_step=0.05, param_max=0.2, kwargs...)
    E = Vector{Vector{Vector{Float64}}}(undef, length(collect(param_min:param_step:param_max)))
    for (i,σ) in enumerate(collect(param_min:param_step:param_max))
        Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
        for j in tqdm(1:num_trials)
            Eᵢ[j] = get_errors(σ, avg_methods; kwargs...)
        end
        E[i] = Eᵢ   
    end
    return E
end

# avg_methods = ["sphere", "A1_L1", "A1_L2", "weiszfeld", "dyadic", "euclidean" ]
# Errs = averaging_sensitivity(avg_methods, param_step=0.025, param_max=0.5, num_trials=1000);
# import GLMakie
# f,a = plot_sensitivity_curves(Errs, varying_parameter="Noise (σ) degrees", avg_methods, parameter= rad2deg.(collect(0:0.025:0.5)))
# GLMakie.save("singularAveraging.png", f)
