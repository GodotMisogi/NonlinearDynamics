## Packages
using DifferentialEquations
using DataFrames
using Plots, LaTeXStrings
using StaticArrays

## ODE definition
function dubby_pendy!(dx, x, ps, t)
    θ₁, θ₂, p₁, p₂     = x
    l₁, l₂, m₁, m₂, g  = ps

    C₀  = l₁ * l₂ * (m₁ + m₂ * sin(θ₁ - θ₂)^2)
    C₁  = (p₁ * p₂ * sin(θ₁ - θ₂)) / C₀
    C₂  = (m₂ * ((l₂ * p₁)^2 + (m₁ + m₂) * (l₁ * p₂)^2) - 2 * l₁ * l₂ * m₂ * p₁ * p₂ * cos(θ₁ - θ₂)) * sin(2(θ₁ - θ₂)) / (2C₀^2)

    dx[1] = (l₂ * p₁ - l₁ * p₂ * cos(θ₁ - θ₂)) / (l₁ * C₀)
    dx[2] = (l₁ * (m₁ + m₂) * p₂ - l₂ * m₂ * p₁ * cos(θ₁ - θ₂)) / (l₂ * m₂ * C₀)
    dx[3] = -(m₁ + m₂) * g * l₁ * sin(θ₁) - C₁ + C₂
    dx[4] = -m₂ * g * l₂ * sin(θ₂) + C₁ - C₂

    nothing
end 

## Inputs and parameters
m1, m2 = 6.0, 8.0  # Masses
l1, l2 = 3.0, 2.0  # Lengths
g      = 9.81      # Gravitational acceleration

ps    = [ l1, l2, m1, m2, g ]
x0    = [ π/2, π/4, 0.1, -0.1 ]
tspan = (0.0, 30.0)
tstep = 60
dt    = maximum(tspan) / tstep

## ODEProblem and solution
run   = ODEFunction(dubby_pendy!, syms = [:θ₁, :θ₂, :p₁, :p₂])
prob  = ODEProblem(run, x0, tspan, ps)
sol   = solve(prob, dt = dt)

## DataFrame
df = DataFrame(sol)

## Plotting
unicodeplots()

## Time-evolution
plt_time_θ = plot(sol, vars = [:θ₁,:θ₂], xlabel = "Time (t)", linewidth = 1.5)
plt_time_p = plot(sol, vars = [:p₁,:p₂], xlabel = "Time (t)", linewidth = 1.5);

## Phase plots
plt_phase_θ = plot(sol, vars = (:θ₁,:θ₂), xlabel = L"$\theta_1$", ylabel = L"$\theta_2$", label = :none)
plt_phase_p = plot(sol, vars = (:p₁,:p₂), xlabel = L"$p_{\theta_1}$", ylabel = L"$p_{\theta_2}$", label = :none);

## Time layout
plt_time  = plot(plt_time_θ, plt_time_p, layout = (1,2), size = (900, 400))

## Phase layout
plt_phase = plot(plt_phase_θ, plt_phase_p, layout = (1,2), size = (900, 400))


## Animation
gr()

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 1, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

## Plotting variables
polar_to_cartesian(r, θ) = SVector(r * cos(θ), r * sin(θ))

ts       = df[!,1]
θ1s, θ2s = df[!,2], df[!,3]
p1s, p2s = df[!,4], df[!,5]

pen_1 = polar_to_cartesian.(l1, θ1s)
pen_2 = pen_1 + polar_to_cartesian.(l2, θ2s)

# Circles
x1, y1 = first.(pen_1), last.(pen_1)
x2, y2 = first.(pen_2), last.(pen_2)

# Lines
m        = 2
lx1, ly1 = range.(0, x1, length = m), range.(0, y1, length = m)
lx2, ly2 = range.(x1, x2, length = m), range.(y1, y2, length = m);

## Loop
plt_phase_θ = plot((x0[1], x0[2]), xlabel = "θ₁", ylabel = "θ₂", color = :darkred, label = :none)
plt_phase_p = plot((x0[3], x0[4]), xlabel = "p₁", ylabel = "p₂", color = :darkblue, label = :none)
plt_time_θ  = plot([0.], [x0[1:2]'], xlabel = "t", label = ["θ₁" "θ₂"]);
plt_time_p  = plot([0.], [x0[3:4]'], xlabel = "t", label = ["p₁" "p₂"]);

n    = length(ts)
anim = @animate for i ∈ 1:n
    # Coordinates
    plt_coords = plot()
    plot!(ly1[i], -lx1[i], c = :black, label = :none)
    plot!(ly2[i], -lx2[i], c = :black, label = :none)
    scatter!((y1[i], -x1[i]), c = :blues, markersize = m1, label = "Bob 1")
    scatter!((y2[i], -x2[i]), c = :reds, markersize = m2, label = "Bob 2")
    circleplot!(y2, -x2, i, line_z = 1:n, cbar = false, c = :reds)

    # # Phases
    push!(plt_phase_θ, θ1s[i], θ2s[i])
    push!(plt_phase_p, p1s[i], p2s[i])

    # # Time
    push!(plt_time_θ, ts[i], [θ1s[i], θ2s[i]])
    push!(plt_time_p, ts[i], [p1s[i], p2s[i]])

    # Layouts
    plt_states = plot(plt_time_θ, plt_time_p, plt_phase_θ, plt_phase_p, layout = (2, 2), size = (800, 800))
    plt = plot(plt_coords, plt_states, layout = (2,1))
end

##
gif(anim, "plots/dubbypendy.gif", fps = 60)