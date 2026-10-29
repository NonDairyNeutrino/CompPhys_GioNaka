using Base: init_regex
using Plots
const VF = Vector{Float64}

function euler(y0 :: VF, dydt :: VF, dt :: Float64)
    y = y0 .+ dydt .* dt
    return y
end

function acceleration()
    ax :: Float64 = 0.0
    ay :: Float64 = -GRAVITY
    a :: VF = Float64[ax, ay]
    return a :: VF
end

function evolve!(position_seq :: Vector{VF}, velocity_seq :: Vector{VF}, time_step :: Float64)
    step = 1
    # find velocity from acceleration
    while all(ri -> ri >= 0, position_seq[step])
        acc = acceleration()
        velocity_old = velocity_seq[step]
        position_old = position_seq[step]

        velocity_new = euler(velocity_old, acc, time_step)
        position_new = euler(position_old, velocity_old, time_step)

        push!(velocity_seq, velocity_new)
        push!(position_seq, position_new)

        step += 1
    end
end

const ANGLE = 30:5:55
const GRAVITY :: Float64 = 9.8
const INIT_POS :: VF = Float64[0.0, 0.0]
const INIT_VEL :: Vector{VF} = (angle -> 700 * [sincosd(angle)...]).(ANGLE) .|> reverse
const TIME_STEP :: Float64 = 0.01

function @main(args :: Vector{String})
    for (i, angle) in enumerate(ANGLE)
        position_seq :: Vector{VF} = VF[INIT_POS]
        velocity_seq :: Vector{VF} = VF[INIT_VEL[i]]

        evolve!(position_seq, velocity_seq, TIME_STEP)

        plot!(
              stack(position_seq)[1, :] ./ 1000, stack(position_seq)[2, :] ./ 1000,
              xlabel = "x (km)",
              ylabel = "y (km)",
              labels = angle
          )
        savefig("trajectories.png")
    end
end
