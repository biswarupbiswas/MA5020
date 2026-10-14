# ==============================================================================
# Linear Advection Equation: u_t + a u_x = 0  (a > 0, periodic domain [0, 1])
# Method: finite volume scheme with ENO3 / WENO5 reconstruction, upwind flux
#         F_{i+1/2} = a u^L_{i+1/2}, and SSP-RK3 time stepping
#         (compared with second-order MUSCL-minmod)
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 13)
# ==============================================================================

using Plots

const eps_weno = 1e-6                 # try 1e-40 and 1e-2 (Exercise E7)

# Third-order values of u at x_{i+1/2} from the stencils S0, S1, S2.
# Arguments: cell averages a = u_{i-2}, b = u_{i-1}, c = u_i, d = u_{i+1}, e = u_{i+2}
function candidates(a, b, c, d, e)
    p0 = (2a - 7b + 11c) / 6
    p1 = (-b + 5c + 2d) / 6
    p2 = (2c + 5d - e) / 6
    return p0, p1, p2
end

# WENO5 (Jiang-Shu): weighted average of the three candidates
function weno5(a, b, c, d, e)
    p0, p1, p2 = candidates(a, b, c, d, e)
    b0 = 13/12 * (a - 2b + c)^2 + 1/4 * (a - 4b + 3c)^2      # smoothness indicators
    b1 = 13/12 * (b - 2c + d)^2 + 1/4 * (b - d)^2
    b2 = 13/12 * (c - 2d + e)^2 + 1/4 * (3c - 4d + e)^2
    a0 = 0.1 / (eps_weno + b0)^2                              # linear weights 1/10, 6/10, 3/10
    a1 = 0.6 / (eps_weno + b1)^2
    a2 = 0.3 / (eps_weno + b2)^2
    return (a0 * p0 + a1 * p1 + a2 * p2) / (a0 + a1 + a2)
end

# ENO3: grow the stencil towards the smaller difference
function eno3(a, b, c, d, e)
    p0, p1, p2 = candidates(a, b, c, d, e)
    if abs(c - b) < abs(d - c)                                # stencil {i-1, i}
        return abs(a - 2b + c) < abs(b - 2c + d) ? p0 : p1
    else                                                      # stencil {i, i+1}
        return abs(b - 2c + d) < abs(c - 2d + e) ? p1 : p2
    end
end

# MUSCL with the minmod slope (Lecture 11), for comparison
minmod(p, q) = p * q <= 0 ? 0.0 : sign(p) * min(abs(p), abs(q))
muscl_minmod(a, b, c, d, e) = c + minmod(c - b, d - c) / 2

# ------------------------------------------------------------------------------
# Choose the reconstruction here: weno5, eno3, muscl_minmod
# ------------------------------------------------------------------------------
recon = weno5

a = 1.0
N = 200                                  # number of cells
dx = 1 / N
x = range(dx/2, 1 - dx/2, length=N)      # cell centres
nu = 0.4                                 # CFL number
dt = nu * dx / a
periods = 1                              # try 10 (Exercise E7)
nsteps = round(Int, periods / (a * dt))

# Initial condition: smooth Gaussian + square wave
u0(x) = exp(-300 * (x - 0.25)^2) + (0.55 < x < 0.8 ? 1.0 : 0.0)

# Semi-discrete operator L(u) = -(F_{i+1/2} - F_{i-1/2}) / dx, periodic boundaries
function L(u, recon)
    N = length(u)
    ug = [u[end-2:end]; u; u[1:3]]       # 3 ghost cells each side: ug[j+3] = u_j
    uL = [recon(ug[i+1], ug[i+2], ug[i+3], ug[i+4], ug[i+5]) for i in 1:N]   # u^L_{i+1/2}
    F = a .* uL
    return -(F .- circshift(F, 1)) ./ dx
end

# One SSP-RK3 step (Shu-Osher)
function ssprk3_step(u, recon)
    u1 = u .+ dt .* L(u, recon)
    u2 = 0.75 .* u .+ 0.25 .* (u1 .+ dt .* L(u1, recon))
    return u ./ 3 .+ (2/3) .* (u2 .+ dt .* L(u2, recon))
end

u = u0.(x)
u_mm = u0.(x)                            # MUSCL-minmod for comparison
for n in 1:nsteps
    global u = ssprk3_step(u, recon)
    global u_mm = ssprk3_step(u_mm, muscl_minmod)

    if n % 5 == 0 || n == nsteps
        xe = mod.(x .- a * n * dt, 1.0)  # exact solution: shifted profile
        plot(x, u0.(xe), label="Exact", color=:black, lw=2, linestyle=:dash)
        plot!(x, u_mm, label="MUSCL-minmod", color=:blue, lw=1.5)
        plot!(x, u, label="$(recon)", color=:red, lw=2)
        display(plot!(xlabel="x", ylabel="u(x,t)", ylims=(-0.2, 1.2),
                      title="$(recon) + SSP-RK3,  t = $(round(n * dt, digits=2))"))
        sleep(0.02)
    end
end

println("Gaussian peak: $(recon) = ", round(maximum(u[x .< 0.45]), digits=4),
        ",  MUSCL-minmod = ", round(maximum(u_mm[x .< 0.45]), digits=4),
        ",  exact = ", round(maximum(u0.(x[x .< 0.45])), digits=4))
println("min(u) = ", minimum(u), "   (negative values: small undershoots, not TVD)")

# Keep the plot window open when run as `julia advection_weno5.jl`
isinteractive() || readline()
