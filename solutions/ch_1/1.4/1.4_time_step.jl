using Plots

function euler(yn :: Number, dydt :: Number, dt :: Number)
    yn + dydt * dt
end

function allocate(time_step_vector :: Vector{Float64}, time_final :: Float64)
    time_matrix :: Vector{Vector{Float64}} = (step -> 0.0:step:10.0).(time_step_vector) .|> collect
    # there should be a velocity for every moment in time, so copy the time and overwrite it
    count_a :: Vector{Vector{Float64}} = deepcopy(time_matrix)
    count_b :: Vector{Vector{Float64}} = deepcopy(time_matrix)
    return (time_matrix, count_a, count_b)
end

function initialize!(
                     count_a :: Vector{Vector{Float64}}, initial_count_a :: Float64,
                     count_b :: Vector{Vector{Float64}}, initial_count_b :: Float64
                    )
    # set the initial value in each count sequence
    for i in 1:length(count_a)
        count_a[i][1] = initial_count_a
        count_b[i][1] = initial_count_b
    end
end

function evolve(count_a_old :: Float64, count_b_old :: Float64, time_step :: Float64)
    count_a_new = euler(count_a_old, decay_rate_a(count_a_old), time_step)
    count_b_new = euler(count_b_old, decay_rate_b(count_a_old, count_b_old), time_step)
    return (count_a_new, count_b_new)
end

function plot_decay(time_matrix, count_a, count_b, time_step_vector)
    p = plot(
             time_matrix / DECAY_A,
             [count_a, count_b] / INITIAL_COUNT_A,
             xlabel = "t/τ_A",
             # xticks = (0:(TIME_FINAL / DECAY_A)),
             ylabel = "N/N_0A",
             labels = ["Count A" "Count B"] #= "dt = " .* string.(transpose(time_step_vector)) =#,
             title = "τ_A / τ_B = " * string(DECAY_A / DECAY_B) # * "    N_0A / N_0B = " * string(INITIAL_COUNT_A / INITIAL_COUNT_B)
             # marker = "o",
             # size   = (1920, 1080)
         )
    savefig(p, "decay.png")
end

const DECAY_A :: Float64 = 1.0
const DECAY_B :: Float64 = 1.0 * DECAY_A
const TIME_FINAL :: Float64 = 10.0 * DECAY_A
# declare constants
const INITIAL_COUNT_A :: Float64 = 100
const INITIAL_COUNT_B :: Float64 = 0.0 * INITIAL_COUNT_A
decay_rate_a(count_a :: Float64) = - count_a / DECAY_A
decay_rate_b(count_a :: Float64, count_b :: Float64) = count_a / DECAY_A - count_b / DECAY_B

function main()
    # declare parameters
    time_step_vector :: Vector{Float64} = [0.01] #= Float64[1,2,5] =# #= logrange(0.01, 0.1, 2) =#

    # allocate
    time_matrix, count_a, count_b = allocate(time_step_vector, TIME_FINAL)

    # initialize
    initialize!(count_a, INITIAL_COUNT_A, count_b, INITIAL_COUNT_B)

    # evolve
    for j in 1:length(time_step_vector)
        time_step = time_step_vector[j]
        for i in 2:length(time_matrix[j])
            count_a_old :: Float64 = count_a[j][i-1]
            count_b_old :: Float64 = count_b[j][i-1]
            count_a[j][i], count_b[j][i] = evolve(count_a_old, count_b_old, time_step)
        end
    end

    plot_decay(time_matrix, count_a, count_b, time_step_vector)
end

main()
