# ==============================================================================
# Inviscid Burgers' Equation: u_t + (u^2/2)_x = 0
# Method: Godunov Finite Volume Scheme (compared with the Rusanov flux)
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 10)
# ==============================================================================

using Plots

f(u) = 0.5 * u^2                      # Burgers' flux

# Godunov flux: exact Riemann solution at the interface (five-branch rule)
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

# Rusanov (Local Lax-Friedrichs) flux
rusanov_flux(uL, uR) = 0.5 * (f(uL) + f(uR)) - 0.5 * max(abs(uL), abs(uR)) * (uR - uL)

# One finite volume step: u_i^{n+1} = u_i^n - (dt/dx) (F_{i+1/2} - F_{i-1/2})
function fv_step(u, lam, flux)
    ug = [u[1]; u; u[end]]            # ghost cells (zero-gradient boundaries)
    F = [flux(ug[i], ug[i+1]) for i in 1:length(u)+1]
    return u .- lam .* (F[2:end] .- F[1:end-1])
end

# ------------------------------------------------------------------------------
# Riemann problem: change u_L, u_R to try other cases
#   u_L = 1.0,  u_R = 0.0 : moving shock
#   u_L = -0.5, u_R = 1.0 : transonic rarefaction
# ------------------------------------------------------------------------------
u_L, u_R = 1.0, 0.0

N = 200                                  # number of cells
dx = 2 / N
x = range(-1 + dx/2, 1 - dx/2, length=N) # cell centres on [-1, 1]
dt = 0.9 * dx / max(abs(u_L), abs(u_R))  # CFL = 0.9
lam = dt / dx
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
u_rus = copy(u_god)
t = 0.0

while t < t_final
    global u_god = fv_step(u_god, lam, godunov_flux)
    global u_rus = fv_step(u_rus, lam, rusanov_flux)
    global t += dt

    plot(x, u_exact.(x, t), label="Exact", color=:black, lw=2, linestyle=:dash)
    plot!(x, u_god, label="Godunov", color=:red, lw=2)
    plot!(x, u_rus, label="Rusanov", color=:blue, lw=2)
    display(plot!(xlabel="x", ylabel="u(x,t)", title="Burgers' equation: t = $(round(t, digits=2))",
                  ylims=(min(u_L, u_R) - 0.2, max(u_L, u_R) + 0.2)))
    sleep(0.03)
end

# Keep the plot window open when run as `julia burgers_godunov.jl`
isinteractive() || readline()
