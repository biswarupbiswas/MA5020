# ==============================================================================
# Inviscid Burgers' Equation: u_t + (u^2/2)_x = 0
# Demonstration: Non-Conservative Schemes Compute the WRONG Shock Position
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 9)
# ==============================================================================
#
# Riemann problem:  u(x,0) = u_L  for x < x0,   u(x,0) = u_R  for x > x0
# with u_L > u_R > 0, so a single shock forms and moves to the right.
#
#   Exact (Rankine-Hugoniot) shock speed :  s_exact = (u_L + u_R) / 2
#   Non-conservative scheme              :  converges to a shock with a DIFFERENT
#                                           speed (measured and printed below)
#
# Both schemes below use the SAME first-order upwind stencil (valid since u > 0):
#
#   Conservative     : u_j^{n+1} = u_j^n - (dt/dx) * ( (u_j^n)^2/2 - (u_{j-1}^n)^2/2 )
#   Non-conservative : u_j^{n+1} = u_j^n - (dt/dx) * u_j^n * ( u_j^n - u_{j-1}^n )
#
# They agree for smooth data but differ by the hidden source term
#   S_j = -(1/(2 dx)) * (u_j - u_{j-1})^2,
# which removes a finite amount of mass at the shock, independent of dx.
#
# Part 1: live animation of both schemes against the exact solution.
# Part 2: grid refinement study showing the error does NOT go away as dx -> 0.
# ==============================================================================

using Plots

# ------------------------------------------------------------------------------
# Problem data
# ------------------------------------------------------------------------------
const u_L = 1.2           # Left state (faster)
const u_R = 0.2           # Right state (slower); try u_R = 0.0 for a frozen shock
const x0 = 0.5            # Initial shock location
const x_min = 0.0
const x_max = 2.0
const t_final = 1.5
const CFL = 0.8           # Courant number based on max|u| = u_L

const s_exact = 0.5 * (u_L + u_R)   # Rankine-Hugoniot speed

u0(x) = x < x0 ? u_L : u_R
u_exact(x, t) = x < x0 + s_exact * t ? u_L : u_R

# ------------------------------------------------------------------------------
# One time step of each scheme
# Boundary conditions: inflow u = u_L on the left (ghost cell), outflow on the right
# (upwind stencil needs only the left neighbour, so no right ghost cell is needed).
# ------------------------------------------------------------------------------
function step_conservative!(u_new, u, lam)
    N = length(u)
    for j in 1:N
        um = (j == 1) ? u_L : u[j-1]
        u_new[j] = u[j] - lam * (0.5 * u[j]^2 - 0.5 * um^2)
    end
    return u_new
end

function step_nonconservative!(u_new, u, lam)
    N = length(u)
    for j in 1:N
        um = (j == 1) ? u_L : u[j-1]
        u_new[j] = u[j] - lam * u[j] * (u[j] - um)
    end
    return u_new
end

# Shock location = point where u crosses the mid value (u_L + u_R)/2,
# found by linear interpolation between neighbouring nodes.
function shock_position(x, u)
    um = 0.5 * (u_L + u_R)
    for j in 1:length(u)-1
        if u[j] >= um && u[j+1] < um
            return x[j] + (um - u[j]) / (u[j+1] - u[j]) * (x[j+1] - x[j])
        end
    end
    return NaN
end

# Run both schemes on an N-point grid up to t_final
function run_both(N)
    dx = (x_max - x_min) / N
    x = collect(range(x_min + dx / 2, x_max - dx / 2, length=N))
    dt = CFL * dx / u_L
    nsteps = round(Int, t_final / dt)
    lam = dt / dx

    uc = u0.(x);  un = u0.(x)
    tmp = similar(uc)
    for n in 1:nsteps
        step_conservative!(tmp, uc, lam);    uc .= tmp
        step_nonconservative!(tmp, un, lam); un .= tmp
    end
    return x, uc, un, nsteps * dt
end

# ==============================================================================
# Part 1: Live animation
# ==============================================================================
function live_demo(; N=200)
    dx = (x_max - x_min) / N
    x = collect(range(x_min + dx / 2, x_max - dx / 2, length=N))
    dt = CFL * dx / u_L
    nsteps = round(Int, t_final / dt)
    lam = dt / dx

    println("----------------------------------------------------------------")
    println("Burgers' Riemann problem: u_L = $u_L, u_R = $u_R, x0 = $x0")
    println("Exact shock speed  s = (u_L + u_R)/2 = $s_exact")
    println("Grid N = $N, dx = $(round(dx, digits=4)), CFL = $CFL, steps = $nsteps")
    println("----------------------------------------------------------------")

    uc = u0.(x);  un = u0.(x)
    tmp = similar(uc)
    xf = range(x_min, x_max, length=2000)   # fine grid for the exact solution
    every = max(1, round(Int, nsteps / 120))

    for n in 1:nsteps
        step_conservative!(tmp, uc, lam);    uc .= tmp
        step_nonconservative!(tmp, un, lam); un .= tmp

        if n % every == 0 || n == nsteps
            t = n * dt
            p = plot(xf, u_exact.(xf, t), label="Exact (s = $s_exact)",
                     color=:black, lw=2.5, linestyle=:dash, legend=:topright)
            plot!(p, x, uc, label="Conservative", color=:blue, lw=2)
            plot!(p, x, un, label="Non-conservative", color=:red, lw=2)
            xlabel!(p, "x"); ylabel!(p, "u(x, t)")
            title!(p, "Burgers' shock: t = $(round(t, digits=2)),  N = $N")
            ylims!(p, 0.0, 1.5); xlims!(p, x_min, x_max)
            display(p)
            sleep(0.02)
        end
    end

    t_end = nsteps * dt
    mass_exact = sum(u_exact.(x, t_end)) * dx
    println("\nAt t = $(round(t_end, digits=3)):")
    println("  Shock position  exact            = $(round(x0 + s_exact * t_end, digits=4))")
    println("  Shock position  conservative     = $(round(shock_position(x, uc), digits=4))")
    println("  Shock position  non-conservative = $(round(shock_position(x, un), digits=4))")
    println("  Shock speed     exact            = $s_exact")
    println("  Shock speed     conservative     = $(round((shock_position(x, uc) - x0) / t_end, digits=4))")
    println("  Shock speed     non-conservative = $(round((shock_position(x, un) - x0) / t_end, digits=4))")
    println("  Total mass      exact            = $(round(mass_exact, digits=4))")
    println("  Total mass      conservative     = $(round(sum(uc) * dx, digits=4))")
    println("  Total mass      non-conservative = $(round(sum(un) * dx, digits=4))  (mass lost!)")
end

# ==============================================================================
# Part 2: Grid refinement: the wrong shock position does NOT converge away
# ==============================================================================
function refinement_study(; Ns=[50, 100, 200, 400, 800, 1600])
    println("\n----------------------------------------------------------------")
    println("Grid refinement study at t = $t_final")
    println("Exact shock position = $(x0 + s_exact * t_final)")
    println("----------------------------------------------------------------")
    println(rpad("N", 8), rpad("x_s (cons)", 16), rpad("x_s (non-cons)", 18), "error (non-cons)")

    p = plot(xlabel="x", ylabel="u(x, t_final)", legend=:topright,
             title="Non-conservative: grid refinement, t = $t_final")
    xf = range(x_min, x_max, length=2000)
    plot!(p, xf, u_exact.(xf, t_final), label="Exact", color=:black, lw=2.5, linestyle=:dash)

    colors = cgrad(:reds, length(Ns), categorical=true)
    for (k, N) in enumerate(Ns)
        x, uc, un, t_end = run_both(N)
        xs_c = shock_position(x, uc)
        xs_n = shock_position(x, un)
        err = abs(xs_n - (x0 + s_exact * t_end))
        println(rpad(N, 8), rpad(round(xs_c, digits=4), 16),
                rpad(round(xs_n, digits=4), 18), round(err, digits=4))
        plot!(p, x, un, label="Non-cons, N = $N", color=colors[k], lw=1.8)
    end
    ylims!(p, 0.0, 1.5); xlims!(p, x_min, x_max)
    display(p)
    println("\nRefining the grid sharpens the non-conservative shock but it stays at the WRONG place.")
end

# ------------------------------------------------------------------------------
# Execute
# ------------------------------------------------------------------------------
live_demo(N=200)
refinement_study()

# Keep the plot window open when run as `julia burgers_nonconservative_shock.jl`
if !isinteractive()
    println("\nPress Enter to exit.")
    readline()
end
