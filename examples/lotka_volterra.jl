## Packages
using DifferentialEquations
using DataFrames
using Plots, LaTeXStrings

## ODE definition
volterra = @ode_def Volterra begin
    dx =  α * x - η * x * y    # Wabbit growth/decay
    dy = -β * y + λ * x * y    # Fox growth/decay
end α β η λ

## Inputs and parameters
α, β  = 2, 1     # Growth/decay parameters
η, λ  = 1, 1     # Coupling parameters
ps    = [ α η ;  
          β λ ]
# x0    = [5, 1]
x0    = [β / λ + 0.01, α / η + 0.01] # Near stable critical point
tspan = (0.0, 180.0)

## ODEProblem and solution
run   = ODEFunction(volterra, syms = [:Wabbits, :Foxes])
prob  = ODEProblem(run, x0, tspan, ps)
sol   = solve(prob)

## DataFrame generation
df = DataFrame(sol)

## Plotting

# PyPlot with LaTeX
# pyplot(dpi = 300, size = (600, 200))
# fntsm = Plots.font("Latin Modern Math", 9)
# fntlg = Plots.font("Latin Modern Math", 10)
# default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)

# PlotlyJS
plotlyjs()

## History plot
plt_time = plot(sol, xlabel = "Time (t)", ylabel = "Population", linewidth = 1.5);
plt_phase = plot(sol, xlabel = "Wabbits", ylabel = "Foxes", vars = (1,2), aspect_ratio = :auto, color = :green, linewidth = 1.5, label = :none);

# savefig(plt_time, "plots/wabbitsfoxes.svg")
# savefig(plt_phase, "plots/wabbitsfoxes_phase.json")

## Combined plot
plt = plot(plt_phase, plt_time, layout = (1,2), size = (900, 400))
# savefig(plt, "plots/wabbitsfoxes-critical.json")