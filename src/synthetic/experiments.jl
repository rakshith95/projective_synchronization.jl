function synchronization_sensitivity(;methods=["spanning-tree", "spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "crossproduct", "euclidean", "weiszfeld" ],vary_parameter="noise", noise_type="elemental", σ_fixed=0.0, num_trials=Int(1e3), parameter_min=0.0, parameter_range=0.1, parameter_max=1.0, error=angular_distance)
    E = Vector{Vector{Vector{Float64}}}(undef, length(collect(parameter_min:parameter_range:parameter_max)))
    if occursin("noise", lowercase(vary_parameter))
        i=1
        for σ=tqdm(parameter_min:parameter_range:parameter_max)
            Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Eⱼ = create_synthetic(σ, noise_type=noise_type, averaging_methods=methods, error=error)    
                if error==angular_distance
                    e = rad2deg.(median.(eachcol(Eⱼ)))
                else
                    e = mean.(eachcol(Eⱼ))
                end
                Eᵢ[j] = e
            end
            E[i] = Eᵢ
            i += 1
        end
        return E
    elseif occursin("density", lowercase(vary_parameter)) && !occursin("outliers", lowercase(vary_parameter))
        i=1
        for ρ=tqdm(parameter_min:parameter_range:parameter_max)
            Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Eⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, holes_density=ρ)    
                if error==angular_distance
                    e = rad2deg.(median.(eachcol(Eⱼ)))
                else
                    e = mean.(eachcol(Eⱼ))
                end 
                Eᵢ[j] = e
            end
            E[i] = Eᵢ
            i += 1
        end
        return E
    elseif occursin("frames", lowercase(vary_parameter))
        i=1
        for n=tqdm(parameter_min:parameter_range:parameter_max)
            Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Eⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, frames=n)    
                if error==angular_distance
                    e = rad2deg.(median.(eachcol(Eⱼ)))
                else
                    e = mean.(eachcol(Eⱼ))
                end
                Eᵢ[j] = e
            end
            E[i] = Eᵢ
            i += 1
        end
        return E
    elseif occursin("outliers", lowercase(vary_parameter))
        i=1
        for Ρ=tqdm(parameter_min:parameter_range:parameter_max)
            Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=1:num_trials
                Eⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, outliers=Ρ)    
                if error==angular_distance
                    e = rad2deg.(median.(eachcol(Eⱼ)))
                else
                    e = mean.(eachcol(Eⱼ))
                end
                Eᵢ[j] = e
            end
            E[i] = Eᵢ
            i += 1
        end
        return E    
    end

end

# avg_methods = ["sphere", "A1"]#, "dyadic", "least-squares-orthogonal", "gaia", "euclidean", "weiszfeld" ]
# avg_methods = ["sphere", "A1", "dyadic", "euclidean", "weiszfeld" ]
# all_methods = ["spanning-tree"; "spectral"; avg_methods]
# E_noise = synchronization_sensitivity(methods=avg_methods, noise_type="angular", parameter_min=0.0, parameter_range=0.05, parameter_max=0.2, num_trials=100, error=angular_distance);
# E_outliers_pointOne = synchronization_sensitivity(methods=avg_methods, vary_parameter="outliers", σ_fixed=0.1 ,parameter_min=0.0, parameter_range=0.1, parameter_max=0.4, num_trials=100, error=angular_distance);
# E_frames_pointOne = synchronization_sensitivity(methods=avg_methods, vary_parameter="frames", σ_fixed=0.1 ,parameter_min=10, parameter_range=15, parameter_max=50, num_trials=100, error=angular_distance);
# p = plot_sensitivity_boxplot(E_noise, ["spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "euclidean", "weiszfeld" ], σ=collect(0:0.05:0.2))

# import GLMakie
# f,a = plot_sensitivity_curves(E_noise, varying_parameter="Noise (σ)", all_methods, parameter=collect(0:0.05:0.2))
# GLMakie.save("outliers_pointOne.png", f)

# M = stack(E[1])'
# median.(eachcol(M))



#=
using JLD2
@save "varying_frames_fixedNoise_pointFive.jld2" E_frames_pointFive
@save "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2" E
@load "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2"  
@load "saved_envs/synchronization_sensitivity_noise_zero_pointOne_point5.jld2"
=#