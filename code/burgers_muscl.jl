# ==============================================================================
# Inviscid Burgers' Equation: u_t + (u^2/2)_x = 0
# Method: Second-order MUSCL finite volume scheme (limited slopes + Godunov flux
#         + SSP-RK2), compared with the first-order Godunov scheme
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 11)
# ==============================================================================

using Plots

f(u) = 0.5 * u^2                      # Burgers' flux

# Godunov flux: exact Riemann solution at the interface (Lecture 10)
function godunov_flux(uL, uR)
    if uL >= uR                       # shock
        return (uL + uR) / 2 > 0 ? f(uL) : f(uR)
    elseif uL >= 0                    # rarefaction moving right
        return f(uL)
    elseif uR <= 0                    # rarefaction moving left
        return f(uR)
    else                              # transonic rarefaction: u* = 0
        return 0.0
    end
end

# Slope limiters: dx * sigma_i from backward jump a and forward jump b
minmod(a, b) = a * b <= 0 ? 0.0 : sign(a) * min(abs(a), abs(b))
zero_slope(a, b) = 0.0                                  # first-order Godunov
mc(a, b) = a * b <= 0 ? 0.0 : sign(a) * min(2abs(a), abs(a + b) / 2, 2abs(b))
centered(a, b) = (a + b) / 2                            # unlimited (try it!)

# Semi-discrete operator L(u) = -(F_{i+1/2} - F_{i-1/2}) / dx
function L(u, dx, limiter)
    N = length(u)
    ug = [u[1]; u[1]; u; u[end]; u[end]]    # 2 ghost cells each side (zero gradient)
    s = [limiter(ug[i] - ug[i-1], ug[i+1] - ug[i]) for i in 2:N+3]   # slopes in cells 0..N+1
    # interface j = 1..N+1 sits between cells j-1 and j
    F = [godunov_flux(ug[j+1] + s[j] / 2, ug[j+2] - s[j+1] / 2) for j in 1:N+1]
    return -(F[2:end] .- F[1:end-1]) ./ dx
end

# One SSP-RK2 step (Heun): average of u^n and two forward Euler steps
function ssprk2_step(u, dt, dx, limiter)
    u1 = u .+ dt .* L(u, dx, limiter)
    return 0.5 .* u .+ 0.5 .* (u1 .+ dt .* L(u1, dx, limiter))
end

# ------------------------------------------------------------------------------
# Riemann problem: change u_L, u_R to try other cases
#   u_L = 1.0,  u_R = 0.0 : moving shock
#   u_L = -0.5, u_R = 1.0 : transonic rarefaction
# ------------------------------------------------------------------------------
u_L, u_R = 1.0, 0.0

N = 100                                  # number of cells
dx = 2 / N
x = range(-1 + dx/2, 1 - dx/2, length=N) # cell centres on [-1, 1]
dt = 0.4 * dx / max(abs(u_L), abs(u_R))  # CFL = 0.4 (TVD needs <= 0.5)
t_final = 0.6

# Exact solution: shock or rarefaction fan
function u_exact(x, t)
    if u_L > u_R
        return x < (u_L + u_R) / 2 * t ? u_L : u_R
    else
        return clamp(x / t, u_L, u_R)
    end
end

u_god = [xi < 0 ? u_L : u_R for xi in x]
u_mm = copy(u_god)
u_mc = copy(u_god)
t = 0.0

while t < t_final
    global u_god = ssprk2_step(u_god, dt, dx, zero_slope)
    global u_mm = ssprk2_step(u_mm, dt, dx, minmod)
    global u_mc = ssprk2_step(u_mc, dt, dx, mc)
    global t += dt

    plot(x, u_exact.(x, t), label="Exact", color=:black, lw=2, linestyle=:dash)
    plot!(x, u_god, label="Godunov (1st order)", color=:gray, lw=2)
    plot!(x, u_mm, label="MUSCL-minmod", color=:blue, lw=2)
    plot!(x, u_mc, label="MUSCL-MC", color=:red, lw=2)
    display(plot!(xlabel="x", ylabel="u(x,t)", title="Burgers' equation: t = $(round(t, digits=2))",
                  ylims=(min(u_L, u_R) - 0.2, max(u_L, u_R) + 0.2)))
    sleep(0.03)
end

# Keep the plot window open when run as `julia burgers_muscl.jl`
isinteractive() || readline()
