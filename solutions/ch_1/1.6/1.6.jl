using Plots

function euler(yn :: Number, dydt :: Number, dt :: Number)
    yn + dydt * dt
end

function allocate(time_step_vector :: Vector{Float64}, time_final :: Float64)
    time_matrix :: Vector{Vector{Float64}} = (step -> 0.0:step:time_final).(time_step_vector) .|> collect
    # there should be a velocity for every moment in time, so copy the time and overwrite it
    population_matrix :: Vector{Vector{Float64}} = deepcopy(time_matrix)
    return (time_matrix, population_matrix)
end

function initialize!(population :: Vector{Vector{Float64}}, init_pop :: Float64)
    # set the initial value in each count sequence
    for i in eachindex(population)
        population[i][1] = init_pop
    end
end

function evolve(pop_old :: Float64, time_step :: Float64)
    pop_new = euler(pop_old, pop_growth(pop_old), time_step)
    return pop_new
end

function plot_pop(time_matrix, pop_matrix, time_step_vector)
    p = plot(
             time_matrix * BIRTH_RATE,
             pop_matrix / INIT_POP,
             ylims = [0,1],
             xlabel = "t*f_B",
             # xticks = (0:(TIME_FINAL / DECAY_A)),
             ylabel = "N/N_0",
             labels = "dt = " .* string.(transpose(time_step_vector ./ BIRTH_RATE)) .* "τ_B",
             # title = "τ_A / τ_B = " * string(DECAY_A / DECAY_B) # * "    N_0A / N_0B = " * string(INITIAL_COUNT_A / INITIAL_COUNT_B)
             # marker = "o",
             # size   = (1920, 1080)
         )
    savefig(p, "population_growth.png")
end

const DEATH_RATE :: Float64 = 0.01
const BIRTH_RATE :: Float64 = 3.0 # 1.0 * DEATH_RATE
const TIME_FINAL :: Float64 = 5.0 / BIRTH_RATE
const INIT_POP :: Float64 = 1000
pop_growth(pop :: Float64) = BIRTH_RATE * pop - DEATH_RATE * pop^2

function main()
    time_step_vector :: Vector{Float64} = (0.01:0.01:0.04) .* BIRTH_RATE
    time_matrix, pop_matrix = allocate(time_step_vector, TIME_FINAL)
    initialize!(pop_matrix, INIT_POP)

    # evolve
    for j in eachindex(time_step_vector)
        time_step = time_step_vector[j]
        for i in 2:length(time_matrix[j])
            pop_old :: Float64 = pop_matrix[j][i-1]
            pop_matrix[j][i] = evolve(pop_old, time_step)
        end
    end

    plot_pop(time_matrix, pop_matrix, time_step_vector)
end

main()
