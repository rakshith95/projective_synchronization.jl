
function synchronization_sensitivity(;num_trials=Int(1e3), σ_min=0.0, σ_range=0.1, σ_max=1.0, error=angular_distance, kwargs...)
    E = Vector{Vector{Vector{Float64}}}(undef, length(collect(σ_min:σ_range:σ_max)))
    i=1
    for σ=tqdm(σ_min:σ_range:σ_max)
        # Eᵢ = MMatrix{num_trials, num_methods, Float64}(zeros(num_trials, num_methods))
        Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
        for j=tqdm(1:num_trials)
            Eⱼ = create_synthetic(σ, error=error,kwargs...)    
            e = rad2deg.(median.(eachcol(Eⱼ)))
            Eᵢ[j] = e
        end
        E[i] = Eᵢ
        i += 1
    end
    return E
end

# methods = ["spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "euclidean", "crossproduct", "weiszfeld" ]
# E = synchronization_sensitivity(σ_range=0.05, σ_max=0.2, num_trials=100, error=angular_distance);
# p = plot_sensitivity_boxplot(E, ["spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "euclidean", "crossproduct", "weiszfeld" ], σ=collect(0:0.1:0.5))

# using JLD2
# @save "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2" E
# @load "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2"  
# @load "saved_envs/synchronization_sensitivity_noise_zero_pointOne_point5.jld2"

# import GLMakie
# f,a = plot_sensitivity_curves(E, ["spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "crossproduct", "euclidean",   "weiszfeld" ], σ=collect(0:0.05:0.2))
# GLMakie.save("fig.png", f)