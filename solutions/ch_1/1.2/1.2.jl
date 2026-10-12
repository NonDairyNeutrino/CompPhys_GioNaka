using Plots

function euler(fn :: Number, dfdt :: Number, dt :: Number)
    fn + dfdt * dt
end

function main()
    # intialize
    VELOCITY :: Float64             = 40
    INITIAL_POSITION :: Float64    = 0.0
    time_step_vector :: Vector{Float64} = Float64[1,2,5] #= logrange(0.001, 1.0, 4) =#
    # allocate
    time_matrix :: Vector{Vector{Float64}} = (step -> 0.0:step:10.0).(time_step_vector) .|> collect
    # there should be a velocity for every moment in time, so copy the time and overwrite it
    position_matrix :: Vector{Vector{Float64}} = deepcopy(time_matrix)

    for i in 1:length(position_matrix)
        position_matrix[i][1] = INITIAL_POSITION
    end

    for j in 1:length(time_step_vector)
        time_step = time_step_vector[j]
        for i in 2:length(position_matrix[j])
            position_matrix[j][i] = euler(position_matrix[j][i-1], VELOCITY, time_step)
        end
    end

    p = plot(
             time_matrix,
             position_matrix,
             xlabel = "t (s)",
             ylabel = "x (m)",
             labels = "dt = " .* string.(transpose(time_step_vector)),
             marker = "o",
             # size   = (1920, 1080)
         )
    savefig(p, "1.2.position.png")
end

main()
