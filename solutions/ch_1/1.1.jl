using Plots

function euler(yn :: Number, dydt :: Number, dt :: Number)
    yn + dydt * dt
end

function main()
    GRAVITY :: Float64             = -9.8
    INITIAL_VELOCITY :: Float64    = 0.0
    # time_step   :: Float64         = 0.1
    time_step_vector :: Vector{Float64} = [0.001, 0.01, 0.1, 1.0]
    time_matrix :: Matrix{Float64} = zeros(Float64, 1, length(time_step_vector))
    velocity_matrix :: Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}} = Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}(undef, size(time_matrix, 2))
    velocity_matrix[1,:] .= INITIAL_VELOCITY

    for j in 1:length(time_step_vector)
        time_step = time_step_vector[j]
        time_range = 0.0:time_step:10.0
        println("using time range: ", time_range)
        time_matrix[:, j] .= Ref(time_range)
        for i in 2:length(time_range)
            velocities[i, j] = euler(velocities[i-1, j], GRAVITY, time_step)
        end
        println("velocity sequence: ", velocity_matrix[:, j])
    end

    # p = plot(time_domain, velocities)
    # @show savefig(p, "velocity.png")
end

main()
