# import GLMakie

function plot_sensitivity_boxplot(Errs::Vector{Vector{Vector{T}}}, methods::AbstractVector{String}; σ=collect(0:0.1:0.5)) where {T <: AbstractFloat}
    num_trials = length(Errs[1])
    if length(Errs[1][1]) != length(methods)
        display("WARNING: Might have more/less methods than expected")
    end
    X_axis = Vector{String}([])
    Errs_matrix = stack(stack.(Errs)')

    for σᵢ in σ
        X_axis = [X_axis; repeat([string(σᵢ)], num_trials)]
    end
    boxes = [box(y=reduce(vcat,Errs_matrix[:,i,:,:]), x=X_axis, name=method,  marker_size=2, fillcolor="white" ) for (i,method) in enumerate(methods)]
    return PlotlyJS.plot(boxes, Layout(legend=attr(x=0.8,y=1), plot_bgcolor="rgb(255,255,255)", yaxis_title="θ deviation (degrees)", xaxis_title="Noise(σ)",  ticklen=2, background=false,  boxmode="group"))
end

function plot_sensitivity_curves(Errs::Vector{Vector{Vector{T}}}, methods::Vector{String}; colors=nothing, varying_parameter="Noise (σ)", parameter=collect(0:0.1:0.5)) where {T <: AbstractFloat}
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1,1], title="Sensitivity Curve")
    
    if length(Errs[1][1]) != length(methods)
        display("WARNING: Might have more/less methods than expected")
    end
    Errs_matrix = stack(stack.(Errs)')
    lin=[]
    for (i,method) in enumerate(methods)
        median_i = mean.(eachcol(Errs_matrix[:,i,1,:]) )
        if isnothing(colors)
            lin = [lin;GLMakie.lines!(ax, parameter, median_i, markersize=rand()*10)]
        else
            lin = [lin;GLMakie.lines!(ax, parameter, median_i, color=colors[i])]
        end
    end
    GLMakie.Legend(fig[1,2], lin, methods)
    ax.xlabel = varying_parameter
    ax.ylabel = "θ deviation (degrees)" 
    return fig,ax
end