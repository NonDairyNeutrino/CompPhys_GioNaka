using Plots

function velocity_no_air(v0 :: Number, power :: Number, mass :: Number, time :: Number)
    return sqrt(v0^2 + 2 * power * time / mass)
end

function allocate(time_step :: Number, time_final :: Number)
    time :: Vector{Float64} = Vector{Float64}(undef, Int(time_final / time_step) + 1)
    velocity :: Vector{Float64} = similar(time)
    return (time, velocity)
end

function initialize!(time, time_step, time_final, velocity, initial_velocity)
    time .= 0.0:time_step:time_final
    velocity[1] = initial_velocity
end

function euler(y0 :: Number, dydt :: Number, dt :: Number)
    return y0 + dydt * dt
end

struct Biker
    mass :: Number
    power :: Number
    area :: Number
    C :: Number
end

struct Environment
    air_density :: Number
    slope :: Number
end

function power(biker :: Biker, velocity :: Number)
    # acceleration for low velocity
    v_star = 7
    f0 = biker.power / v_star
    lv = f0 / biker.mass
    # acceleration for high velocity
    hv = biker.power / (biker.mass * velocity)
    # biker power is independent of velocity when biker is slow
    p = f0 * velocity < biker.power ? lv : hv
    return p
end

function acceleration(biker :: Biker, env :: Environment, velocity :: Number)
    GRAVITY = 9.8
    THETA = atan(env.slope)
    p = power(biker, velocity)
    grav = GRAVITY * sin(THETA)
    b2 = biker.C * env.air_density * biker.area * velocity^2 / (2 * biker.mass)
    return p + grav - b2
end

function acceleration(biker :: Biker, env :: Environment)
    return (v -> acceleration(biker, env, v))
end

function evolve!(acc :: Function, velocity :: Vector{Float64}, time_step :: Number)
    threshold :: Float64 = 10^(-10)
    push!(velocity, euler(velocity[1], acc(velocity[1]), time_step))
    i :: UInt = 2
    while abs(velocity[i] - velocity[i - 1]) > threshold
        i += 1
        v_old = velocity[i - 1]
        v_new = euler(v_old, acc(v_old), time_step)
        push!(velocity, v_new)
    end
    # println("Threshold met after $i steps; Terminal velocity: ", velocity[end])
end

function mps2mph(mps :: Number)
    MILES_PER_METER = 1 / 1609
    SECONDS_PER_HOUR = 3600
    return mps * MILES_PER_METER * SECONDS_PER_HOUR
end

function @main(args :: Vector{String})
    MASS = 70
    POWER = 400
    AREA = 0.25:0.01:0.5
    C = 1
    AIR_DENSITY = 1
    SLOPE = 0.0:0.01:0.25
    INITIAL_VELOCITY = 0
    TIME_STEP = 0.01
    TIME_FINAL = 200
    terminal_velocity = Matrix{Float64}(undef, length(AREA), length(SLOPE))
    v = Float64[]

    for s in eachindex(SLOPE)
        for a in eachindex(AREA)
            # makes a new biker object every iteration when one could be used in-place
            biker :: Biker = Biker(MASS, POWER, AREA[a], C)
            env :: Environment = Environment(AIR_DENSITY, SLOPE[s])
            acc = acceleration(biker, env)
            velocity :: Vector{Float64} = Float64[INITIAL_VELOCITY]

            println("Evolving for (slope, area) = ($(SLOPE[s]), $(AREA[a]))")
            evolve!(acc, velocity, TIME_STEP)

            if SLOPE[s] == 0.0 && AREA[a] == 0.33
                v = velocity
            end
            terminal_velocity[a, s] = velocity[end]
        end
    end

    time = range(start = 0.0, step = 0.1, length = length(v))
    p = plot(
        time, v
    )
    savefig(p, "velocity.png")

    tvhm = heatmap(
         SLOPE, AREA, terminal_velocity .|> mps2mph,
         xlabel = "tan(θ) (rad)",
         ylabel = "Area (m^2)"
    )
    savefig(tvhm, "terminal_velocity_heatmap.png")
    return;
end
