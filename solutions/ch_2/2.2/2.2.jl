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

function acceleration(biker :: Biker, air_density :: Number, velocity :: Number)
    return biker.power / (biker.mass * velocity) - biker.C * air_density * biker.area * velocity^2 / (2 * biker.mass)
end

function acceleration(biker :: Biker, air_density :: Number)
    return (v -> acceleration(biker, air_density, v))
end

function evolve!(acc :: Function, velocity :: Vector{Float64}, time_step :: Number)
    threshold :: Float64 = 10^(-10)
    i :: UInt = 2
    push!(velocity, euler(velocity[1], acc(velocity[1]), time_step))
    while abs(velocity[i] - velocity[i - 1]) > threshold
        i += 1
        v_old = velocity[i - 1]
        v_new = euler(v_old, acc(v_old), time_step)
        push!(velocity, v_new)
    end
end

function @main(args :: Vector{String})
    MASS = 70
    POWER = 400
    AREA = 0.33
    C = 1
    AIR_DENSITY = 1
    INITIAL_VELOCITY = 4
    TIME_STEP = 0.01
    TIME_FINAL = 200
    biker = Biker(MASS, POWER, AREA, C)

    # time :: Vector{Float64}, velocity :: Vector{Float64} = allocate(TIME_STEP, TIME_FINAL)
    # initialize!(time, TIME_STEP, TIME_FINAL, velocity, INITIAL_VELOCITY)
    velocity = Float64[INITIAL_VELOCITY]

    acc = acceleration(biker, AIR_DENSITY)
    evolve!(acc, velocity, TIME_STEP)
    time = range(; start = 0.0, length = length(velocity), step = TIME_STEP)

    velocity_airless = velocity_no_air.(INITIAL_VELOCITY, biker.power, biker.mass, time)

    plot(
         time, [velocity velocity_airless],
         xlabel = "time (s)",
         ylabel = "velocity (m/s)",
         labels = ["with air" "no air"]
     )
    savefig("velocity.png")
    return;
end
