using Plots

function euler(yn :: Number, dydt :: Number, dt :: Number)
    yn + dydt * dt
end

function allocate(nratios :: Int, nsteps :: Int)
    time_matrix :: Matrix{Float64} = Matrix{Float64}(undef, nsteps + 1, nratios)
    count_a :: Matrix{Float64} = similar(time_matrix)
    count_b :: Matrix{Float64} = similar(time_matrix)
    return (time_matrix, count_a, count_b)
end

function initialize!(
                     time_matrix :: Matrix{Float64},
                     count_a :: Matrix{Float64}, initial_count_a :: Float64,
                     count_b :: Matrix{Float64}, initial_count_b :: Float64
                    )
    # set the initial value in each count sequence
    eachcol(time_matrix) .= fill(collect(0:TIME_STEP:TIME_FINAL), size(time_matrix, 2))
    count_a[1, :] .= initial_count_a
    count_b[1, :] .= initial_count_b
end

function evolve(count_a_old :: Float64, count_b_old :: Float64, time_step :: Float64)
    count_a_new = euler(count_a_old, decay_rate_a(count_a_old, count_b_old), time_step)
    count_b_new = euler(count_b_old, decay_rate_b(count_a_old, count_b_old), time_step)
    return (count_a_new, count_b_new)
end

function plot_decay(time_matrix, count_a, count_b)
    p = plot(
             time_matrix / DECAY,
             [count_a, count_b] / INITIAL_COUNT_A,
             xlabel = "t/τ",
             ylabel = "N/N_A0",
             labels = ["Count A" "Count B"]
         )
    savefig(p, "counts.png")
end

const DECAY :: Float64 = 1.0
const TIME_FINAL :: Float64 = 10.0 * DECAY
const TIME_STEP :: Float64 = 0.01 * DECAY
# declare constants
const INITIAL_COUNT_A :: Float64 = 100
const INITIAL_COUNT_B :: Float64 = 0.0 * INITIAL_COUNT_A
decay_rate_a(count_a :: Float64, count_b :: Float64) = - count_a / DECAY + count_b / DECAY
decay_rate_b(count_a :: Float64, count_b :: Float64) = count_a / DECAY - count_b / DECAY

function main()
    # allocate
    time_matrix, count_a, count_b = allocate(1, Int(TIME_FINAL / TIME_STEP))

    # initialize
    initialize!(time_matrix, count_a, INITIAL_COUNT_A, count_b, INITIAL_COUNT_B)

    # evolve
    j = 1
    for i in 2:length(time_matrix[:, j])
        count_a_old :: Float64 = count_a[i-1, j]
        count_b_old :: Float64 = count_b[i-1, j]
        count_a[i, j], count_b[i, j] = evolve(count_a_old, count_b_old, TIME_STEP)
    end

    plot_decay(time_matrix, count_a[:,1], count_b)
end

main()
