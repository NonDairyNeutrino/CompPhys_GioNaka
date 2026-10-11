using Plots

function euler(yn :: Number, dydt :: Number, dt :: Number)
    yn + dydt * dt
end

function main()
    # intialize
    GRAVITY :: Float64             = -9.8
    INITIAL_VELOCITY :: Float64    = 0.0
    time_step_vector :: Vector{Float64} = logrange(0.001, 1.0, 4)
    # allocate
    time_matrix :: Vector{Vector{Float64}} = (step -> 0.0:step:10.0).(time_step_vector) .|> collect
    # there should be a velocity for every moment in time, so copy the time and overwrite it
    velocity_matrix :: Vector{Vector{Float64}} = deepcopy(time_matrix)

    for i in 1:length(velocity_matrix)
        velocity_matrix[i][1] = INITIAL_VELOCITY
    end

    @show length.(velocity_matrix)
    for j in 1:length(time_step_vector)
        time_step = time_step_vector[j]
        for i in 2:length(velocity_matrix[j])
            velocity_matrix[j][i] = euler(velocity_matrix[j][i-1], GRAVITY, time_step)
        end
        # println("velocity sequence: ", velocity_matrix[:, j])
    end

    p = plot(time_matrix, velocity_matrix)
    @show savefig(p, "velocity.png")
end

main()
