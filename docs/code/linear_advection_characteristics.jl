# Linear Advection Equation: u_t + a*u_x = 0
using Plots

a = 1.5
u0(x) = exp(-10.0 * (x - 0.5)^2)
x = range(-1.0, 5.0, length=200)

for t in 0.0:0.04:2.0
    u = [u0(xi - a * t) for xi in x]
    display(plot(x, u, ylims=(0, 1.2), xlabel="x", ylabel="u(x,t)",
                 title="Linear Advection: t = $(round(t, digits=2))", legend=false))
    sleep(0.02)
end
