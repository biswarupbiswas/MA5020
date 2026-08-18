# Inviscid Burgers Equation: u_t + u*u_x = 0
using Plots

u0(x) = exp(-x^2)
u0_prime(x) = -2.0 * x * exp(-x^2)
t_s = -1.0 / minimum([u0_prime(xi) for xi in range(-3.0, 3.0, length=1000)])

function find_x0(x, t)
    a, b = x - 3.0*(t + 1), x + 3.0*(t + 1)
    for _ in 1:60
        c = (a + b) / 2.0
        (a + u0(a)*t - x) * (c + u0(c)*t - x) < 0 ? b = c : a = c
    end
    return (a + b) / 2.0
end

x = range(-3.0, 4.0, length=200)

for t in 0.0:0.02:0.99*t_s
    u = [u0(find_x0(xi, t)) for xi in x]
    display(plot(x, u, ylims=(0, 1.2), xlabel="x", ylabel="u(x,t)",
                 title="Burgers Wave Steepening: t = $(round(t, digits=2))", legend=false))
    sleep(0.02)
end
