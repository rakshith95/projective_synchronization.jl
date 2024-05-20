function synchronization_sensitivity(;methods=["spanning-tree", "spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "crossproduct", "euclidean", "weiszfeld" ],vary_parameter="noise", noise_type="elemental", σ_fixed=0.0, num_trials=Int(1e3), parameter_min=0.0, parameter_range=0.1, parameter_max=1.0, error=angular_distance)
    E = Vector{Vector{Vector{Float64}}}(undef, length(collect(parameter_min:parameter_range:parameter_max)))
    if occursin("noise", lowercase(vary_parameter))
        i=1
        for σ=tqdm(parameter_min:parameter_range:parameter_max)
            Eᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Eⱼ = create_synthetic(σ, noise_type=noise_type, averaging_methods=methods, error=error, anchor="centrality", update="all-random")    
                if error==angular_distance
                    e = rad2deg.(mean.(eachcol(Eⱼ)))
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
                Eⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, holes_density=ρ, anchor="centrality", update="all-random")  
                if error==angular_distance
                    e = rad2deg.(mean.(eachcol(Eⱼ)))
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
                Eⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, frames=n, anchor="centrality", update="all-random")
                if error==angular_distance
                    e = rad2deg.(mean.(eachcol(Eⱼ)))
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
            for j=tqdm(1:num_trials)
                Eⱼ = create_synthetic(σ_fixed, holes_density=0.8, noise_type=noise_type, averaging_methods=methods, error=error, outliers=Ρ, anchor="centrality", update="all-random")    
                if error==angular_distance
                    e = rad2deg.(mean.(eachcol(Eⱼ)))
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

function synchronization_timer(vary_parameter;methods=["spanning-tree", "spectral", "sphere", "dyadic", "euclidean", "weiszfeld" ], noise_type="angular", σ_fixed=0.1, num_trials=Int(1e3), parameter_min=0.0, parameter_range=0.1, parameter_max=1.0, error=angular_distance)
    T = Vector{Vector{Vector{Float64}}}(undef, length(collect(parameter_min:parameter_range:parameter_max)))
    
    if occursin("density", lowercase(vary_parameter)) && !occursin("outliers", lowercase(vary_parameter))
        i=1
        for ρ=tqdm(parameter_min:parameter_range:parameter_max)
            Tᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Tᵢ[j] = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, holes_density=ρ, anchor="centrality", update="all-random")
            end
            T[i] = Tᵢ
            i += 1
        end
        return T
    elseif occursin("frames", lowercase(vary_parameter))
        i=1
        for n=tqdm(parameter_min:parameter_range:parameter_max)
            Tᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=tqdm(1:num_trials)
                Tⱼ = create_synthetic(σ_fixed, noise_type=noise_type, averaging_methods=methods, error=error, frames=n, anchor="centrality", update="all-random")
                Tᵢ[j] = Tⱼ
            end
            T[i] = Tᵢ
            i += 1
        end
        return T
    elseif occursin("noise", lowercase(vary_parameter))
        i=1
        for σ=parameter_min:parameter_range:parameter_max
            Tᵢ = Vector{Vector{Float64}}(undef, num_trials)
            for j=1:num_trials
                Tⱼ = create_synthetic(σ, noise_type=noise_type, averaging_methods=methods, error=error, anchor="centrality", update="all-random")
                Tᵢ[j] = Tⱼ
            end
            T[i] = Tᵢ
            i += 1
        end
        return T
    end
end
    

# avg_methods = ["sphere", "sphere-init", "weiszfeld", "weiszfeld-init"];
# all_methods = ["spanning-tree"; "spectral"; avg_methods];
# E_noise_initST = synchronization_sensitivity(methods=avg_methods, noise_type="angular", parameter_min=0.0, parameter_range=0.025, parameter_max=0.2, num_trials=200, error=angular_distance);
# T_noise_initST = synchronization_timer("noise"; methods=avg_methods, parameter_min=0, parameter_range=0.025, parameter_max=0.2, num_trials=100, error=angular_distance);
# E_outliers_Zero_initST = synchronization_sensitivity(methods=avg_methods, noise_type="angular", vary_parameter="outliers", σ_fixed=0.0 ,parameter_min=0.0, parameter_range=0.1, parameter_max=0.6, num_trials=100, error=angular_distance);
# E_outliers_pointOne = synchronization_sensitivity(methods=avg_methods, noise_type="angular", vary_parameter="outliers", σ_fixed=0.1, parameter_min=0.0, parameter_range=0.1, parameter_max=0.6, num_trials=1000, error=angular_distance);




# all_methods = ["spanning-tree"; "spectral"; avg_methods];
# T_holesDensity_pointOne = synchronization_timer("density";methods=avg_methods, σ_fixed=0.1 ,parameter_min=0.0, parameter_range=0.05, parameter_max=0.95, num_trials=1000, error=angular_distance);
# avg_methods = ["dyadic", "sphere", "weiszfeld"];
# T_frames_pointOne = synchronization_timer("frames"; methods=avg_methods,  σ_fixed=0.1 ,parameter_min=10, parameter_range=15, parameter_max=100, num_trials=1000, error=angular_distance);
# Ts_matrix = stack(stack.(T_frames_pointOne)');
# Ts_matrix = dropdims(Ts_matrix, dims = tuple(findall(size(Ts_matrix) .== 1)...));
# file = MAT.matopen("times_new_frames_final.mat", "w")
# write(file, "T", Ts_matrix)
# close(file)





# E_noise = synchronization_sensitivity(methods=avg_methods, noise_type="angular", parameter_min=0.0, parameter_range=0.025, parameter_max=0.2, num_trials=1000, error=angular_distance);
# E_holesDensity_pointOne = synchronization_sensitivity(;methods=avg_methods, noise_type="angular", vary_parameter="dens    ity", σ_fixed=0.1 ,parameter_min=0.0, parameter_range=0.05, parameter_max=0.95, num_trials=1000, error=angular_distance);
# E_outliers_Zero = synchronization_sensitivity(methods=avg_methods, noise_type="angular", vary_parameter="outliers", σ_fixed=0.0 ,parameter_min=0.0, parameter_range=0.1, parameter_max=0.6, num_trials=1000, error=angular_distance);
# E_outliers_pointOne = synchronization_sensitivity(methods=avg_methods, noise_type="angular", vary_parameter="outliers", σ_fixed=0.1, parameter_min=0.0, parameter_range=0.1, parameter_max=0.6, num_trials=1000, error=angular_distance);
# E_frames_pointOne = synchronization_sensitivity(methods=avg_methods, noise_type="angular", vary_parameter="frames", σ_fixed=0.1 ,parameter_min=10, parameter_range=15, parameter_max=100, num_trials=1000, error=angular_distance);
# p = plot_sensitivity_boxplot(E_noise, ["spectral", "sphere", "dyadic", "least-squares-orthogonal", "gaia", "euclidean", "weiszfeld" ], σ=collect(0:0.05:0.2))

# import GLMakie
# all_methods[4] = "A1-L1"
# colors = ["green", "red", "purple", "yellow", "black", "blue" ];
# f,a = plot_sensitivity_curves(E_outliers_Zero_initST, varying_parameter="Noise Density", all_methods, parameter=collect(0.0:0.1:0.6),)
# GLMakie.save("Holes_outliers.png", f)

# M = stack(E[1])'
# median.(eachcol(M))

# Errs_matrix = stack(stack.(E_noise_initST)');
# Errs_matrix = dropdims(Errs_matrix, dims = tuple(findall(size(Errs_matrix) .== 1)...));
# file = MAT.matopen("Noise_initST_0_point025_point2_100Trials.mat", "w")
# write(file, "E", Errs_matrix)
# close(file)

# Ts_matrix = stack(stack.(T_noise_initST)');
# Ts_matrix = dropdims(Ts_matrix, dims = tuple(findall(size(Ts_matrix) .== 1)...));
# file = MAT.matopen("times_initST_100trials.mat", "w")
# write(file, "T", Ts_matrix)
# close(file)


#=
using JLD2
@save "varying_frames_fixedNoise_pointFive.jld2" E_frames_pointFive
@save "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2" E
@load "synchronization_sensitivity_noise_zero_to_zeroPointtwo.jld2"  
@load "saved_envs/synchronization_sensitivity_noise_zero_pointOne_point5.jld2"
=#


