##
f1(x, y, α, η) =  α * x - η * x * y
f2(x, y, β, λ) = -β * y + λ * x * y


function volterra!(R, xs, ps, t)
    R[1] = f1(xs[1], xs[2], ps[1,1], ps[1,2])
    R[2] = f2(xs[1], xs[2], ps[2,1], ps[2,2])

    nothing
end

##
α, η, β, λ = 0.1, 0.1, 0.2, 0.3
ps = [ α η ; 
       β λ ]
x0 = [ 2.0, 2.0 ]

##
using DifferentialEquations

tspan = (0.0, 500.0)
run   = ODEFunction(volterra!, syms = [:Wabbits, :Foxes])
prob  = ODEProblem(run, x0, tspan, ps)
sol   = solve(prob, maxiters = 500)

##
using DataFrames

df = DataFrame(sol)

##
using Plots
plotly(dpi = 300, size = (800, 600))

plt = plot(sol,linewidth=2)
savefig(plt, "wabbitsfoxes.json")