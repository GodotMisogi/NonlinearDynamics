## Packages
using DifferentialEquations
using DataFrames
using Plots, LaTeXStrings

## ODE definition
"""
Falkner-Skan ODE:

    F = ψ/m
    U = u/uₑ
    S = τξ/(ρuₑ²δ) = τ/(mu/ξ)
    
    where m = ρₑuₑδ
"""
falkner_skan = @ode_def FalknerSkan begin
    dF  = U                                     # F' = U
    dU  = S                                     # U' = F'' = S
    dS  = -(1 + a) / 2 * F * S - a * (1 - U^2)  # S' = F''' = -(1 + a)/2 * F * F'' - a(1 - (F')²)
end a

function falkner_skan_BC!(bcs, x, p, η)
    bcs[1] = x[1][1]        # η = 0,  F  = 0
    bcs[2] = x[1][2]        # η = 0,  F' = 0
    bcs[3] = x[end][2] - 1  # η = ηₑ, F' = 1

    nothing
end

##
ηs    = (0.0, 12.0)
x0    = [5.0, 2.0, 0.0] 
as    = sort([ range(-0.11, 0.2, length = 10); -0.092 ])

##
prob  = BVProblem.(Ref(falkner_skan), Ref(falkner_skan_BC!), Ref(x0), Ref(ηs), as, syms = Ref([:y, :δ]))
sol   = solve.(prob, Ref(GeneralMIRK4()), dt = 5e-1)

##
dfs  = DataFrame.(sol)

##
unicodeplots()

##
plt = plot()
map((solution, a) -> plot!(solution, label = "a = $a", xaxis = "U = F'", yaxis = "η", linewidth = 1.5, vars = [(2, 0)]), sol, as)
plot!()
# savefig(plt, "plots/falknerskan.json")