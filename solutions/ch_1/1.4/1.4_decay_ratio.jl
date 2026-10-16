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

function evolve(count_a_old :: Float64, count_b_old :: Float64, decay_ratio :: Float64, time_step :: Float64)
    count_a_new = euler(count_a_old, decay_rate_a(count_a_old), time_step)
    count_b_new = euler(count_b_old, decay_rate_b(count_a_old, count_b_old, decay_ratio), time_step)
    return (count_a_new, count_b_new)
end

function plot_decay(time_matrix, count_a, count_b, decay_ratio_vector)
    p = plot(
             time_matrix / DECAY_A,
             [count_a, count_b] / INITIAL_COUNT_A,
             xlabel = "t/τ_A",
             ylabel = "N/N_A0",
             labels = ["Count A" ((ratio -> "Count B - τ_B / τ_A = $ratio").(decay_ratio_vector))...],
             # title = "τ_A / τ_B = " * string(DECAY_A / DECAY_B)
         )
    savefig(p, "decay.png")
end

const DECAY_A :: Float64 = 1.0
# const DECAY_B :: Float64 = DECAY_A
const TIME_FINAL :: Float64 = 10.0 * DECAY_A
const TIME_STEP :: Float64 = 0.1
# declare constants
const INITIAL_COUNT_A :: Float64 = 100
const INITIAL_COUNT_B :: Float64 = 0.0 * INITIAL_COUNT_A
decay_rate_a(count_a :: Float64) = - count_a / DECAY_A
decay_rate_b(count_a :: Float64, count_b :: Float64, decay_ratio :: Float64) = count_a / DECAY_A - count_b / (decay_ratio * DECAY_A)

function main()
    # declare parameters
    decay_ratio_vector :: Vector{Float64} = [0.5, 1.0, 2.0]
    # @show decay_ratio_vector

    # allocate
    time_matrix, count_a, count_b = allocate(length(decay_ratio_vector), Int(TIME_FINAL / TIME_STEP))

    # initialize
    initialize!(time_matrix, count_a, INITIAL_COUNT_A, count_b, INITIAL_COUNT_B)
    # @show time_matrix count_a count_b

    # evolve
    Threads.@threads for j in 1:length(decay_ratio_vector)
        decrat = decay_ratio_vector[j]
        @info "Processing decay ratio" decrat
        for i in 2:length(time_matrix[:, j])
            count_a_old :: Float64 = count_a[i-1, j]
            count_b_old :: Float64 = count_b[i-1, j]
            count_a[i, j], count_b[i, j] = evolve(count_a_old, count_b_old, decay_ratio_vector[j], TIME_STEP)
        end
    end

    plot_decay(time_matrix, count_a[:,1], count_b, decay_ratio_vector)
end

main()
