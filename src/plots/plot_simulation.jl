using PlotlyJS

function plot_sensitivity_boxplot(Errs::Vector{Vector{Vector{T}}}, methods::AbstractVector{String}; σ=collect(0:0.1:0.5)) where {T <: AbstractFloat}
    num_trials = length(E[1])
    if length(Errs[1][1]) != length(methods)
        display("WARNING: Might have more/less methods than expected")
    end
    X_axis = Vector{String}([])
    Y = Vector{Vector{}}
    Errs_matrix = stack(stack.(Errs)')

    for σᵢ in σ
        X_axis = [X_axis; repeat([string(σᵢ)], num_trials)]
    end
    boxes = [box(y=reduce(vcat,Errs_matrix[:,i,:,:]), x=X_axis, name=method,  marker_size=2, fillcolor="white" ) for (i,method) in enumerate(methods)]
    return PlotlyJS.plot(boxes, Layout(legend=attr(x=0.8,y=1), plot_bgcolor="rgb(255,255,255)", yaxis_title="θ deviation (degrees)", xaxis_title="Noise(σ)",  ticklen=2, background=false,  boxmode="group"))
end