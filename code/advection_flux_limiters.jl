# ==============================================================================
# Linear Advection Equation: u_t + a u_x = 0  (a > 0, periodic domain [0, 1])
# Method: TVD flux-limiter schemes
#   F_{i+1/2} = a u_i + (a/2)(1 - nu) phi(r_i) (u_{i+1} - u_i),
#   r_i = (u_i - u_{i-1}) / (u_{i+1} - u_i)
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 12)
# ==============================================================================

using Plots

# Flux limiters phi(r)  (Sweby TVD region: phi = 0 for r <= 0, phi <= min(2r, 2))
minmod(r)   = max(0.0, min(1.0, r))
superbee(r) = max(0.0, min(2r, 1.0), min(r, 2.0))
vanleer(r)  = (r + abs(r)) / (1 + abs(r))
mc(r)       = max(0.0, min(2r, (1 + r) / 2, 2.0))
laxwendroff(r) = 1.0                  # not TVD: try it!

# ------------------------------------------------------------------------------
# Choose the limiter here: minmod, superbee, vanleer, mc, laxwendroff
# ------------------------------------------------------------------------------
phi = mc

a = 1.0
N = 200                                  # number of cells
dx = 1 / N
x = range(dx/2, 1 - dx/2, length=N)      # cell centres
nu = 0.5                                 # CFL number
dt = nu * dx / a
nsteps = round(Int, 1 / (a * dt))        # one full period

# Initial condition: smooth Gaussian + square wave
u0(x) = exp(-300 * (x - 0.25)^2) + (0.55 < x < 0.8 ? 1.0 : 0.0)

# One step of the flux-limiter scheme (periodic boundaries via circshift)
function step(u, phi)
    du = circshift(u, -1) .- u                      # u_{i+1} - u_i
    dm = u .- circshift(u, 1)                       # u_i - u_{i-1}
    r = [du[i] == 0 ? 0.0 : dm[i] / du[i] for i in eachindex(u)]
    F = a .* u .+ 0.5 * a * (1 - nu) .* phi.(r) .* du
    return u .- (dt / dx) .* (F .- circshift(F, 1))
end

u = u0.(x)
u_up = u0.(x)                            # first-order upwind for comparison
for n in 1:nsteps
    global u = step(u, phi)
    global u_up = step(u_up, r -> 0.0)

    if n % 4 == 0 || n == nsteps
        xe = mod.(x .- a * n * dt, 1.0)  # exact solution: shifted profile
        plot(x, u0.(xe), label="Exact", color=:black, lw=2, linestyle=:dash)
        plot!(x, u_up, label="Upwind", color=:gray, lw=1.5)
        plot!(x, u, label="$(phi)", color=:red, lw=2)
        display(plot!(xlabel="x", ylabel="u(x,t)", ylims=(-0.2, 1.2),
                      title="Flux limiter: $(phi),  t = $(round(n * dt, digits=2))"))
        sleep(0.02)
    end
end

println("TV(u) after one period = ", sum(abs.(circshift(u, -1) .- u)),
        "   (initial TV = ", sum(abs.(circshift(u0.(x), -1) .- u0.(x))), ")")

# Keep the plot window open when run as `julia advection_flux_limiters.jl`
isinteractive() || readline()
