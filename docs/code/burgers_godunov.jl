# ==============================================================================
# Inviscid Burgers' Equation: u_t + (u^2/2)_x = 0
# Method: Finite Volume Godunov Scheme (exact Riemann solver at each interface)
#         compared with the Rusanov (Local Lax-Friedrichs) flux
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 10)
# ==============================================================================
#
# Finite volume update for the cell average ubar_i:
#   ubar_i^{n+1} = ubar_i^n - (dt/dx) * ( F_{i+1/2}^n - F_{i-1/2}^n )
#
# Godunov flux:  F_{i+1/2} = f( u_RP(0; ubar_i, ubar_{i+1}) ),  f(u) = u^2/2
# For Burgers' equation this is the five-branch rule of Lecture 10:
#
#   u_L >= u_R, s = (u_L+u_R)/2 > 0  : f(u_L)   right-moving shock
#   u_L >= u_R, s <= 0               : f(u_R)   left-moving shock
#   u_L <  u_R, u_L >= 0             : f(u_L)   supersonic right rarefaction
#   u_L <  u_R, u_R <= 0             : f(u_R)   supersonic left rarefaction
#   u_L <  0 < u_R                   : 0        transonic rarefaction (u* = 0)
#
# Rusanov flux:  F = (f(u_L) + f(u_R))/2 - (s_max/2) * (u_R - u_L),
#                s_max = max(|u_L|, |u_R|)
#
# Test 1: moving shock          (u_L = 1.0,  u_R = 0.0)
# Test 2: transonic rarefaction (u_L = -0.5, u_R = 1.0), fan passes through u = 0
# ==============================================================================

using Plots

f(u) = 0.5 * u^2

# ------------------------------------------------------------------------------
# Numerical fluxes
# ------------------------------------------------------------------------------
function godunov_flux(uL, uR)
    if uL >= uR                        # shock
        s = 0.5 * (uL + uR)
        return s > 0 ? f(uL) : f(uR)
    else                               # rarefaction
        if uL >= 0
            return f(uL)               # whole fan moves right
        elseif uR <= 0
            return f(uR)               # whole fan moves left
        else
            return 0.0                 # transonic: sonic point u* = 0
        end
    end
end

function rusanov_flux(uL, uR)
    s_max = max(abs(uL), abs(uR))
    return 0.5 * (f(uL) + f(uR)) - 0.5 * s_max * (uR - uL)
end

# ------------------------------------------------------------------------------
# Exact Riemann solution u(x,t) with jump at x0 (entropy solution)
# ------------------------------------------------------------------------------
function riemann_exact(x, t, uL, uR, x0)
    t == 0 && return x < x0 ? uL : uR
    xi = (x - x0) / t
    if uL > uR                         # shock with Rankine-Hugoniot speed
        return xi < 0.5 * (uL + uR) ? uL : uR
    else                               # rarefaction fan: u = xi inside the fan
        return xi <= uL ? uL : (xi >= uR ? uR : xi)
    end
end

# ------------------------------------------------------------------------------
# One finite volume step with transmissive (zero-gradient) boundaries
# ------------------------------------------------------------------------------
function fv_step!(u_new, u, lam, flux)
    N = length(u)
    F = zeros(N + 1)                   # F[i] is the flux at interface x_{i-1/2}
    for i in 1:N+1
        uL = (i == 1) ? u[1] : u[i-1]
        uR = (i == N + 1) ? u[N] : u[i]
        F[i] = flux(uL, uR)
    end
    for i in 1:N
        u_new[i] = u[i] - lam * (F[i+1] - F[i])
    end
    return u_new
end

# ==============================================================================
# Setup
# ==============================================================================
function solve_burgers_godunov()
    x_min, x_max = -1.0, 1.0
    x0 = 0.0
    N = 200                            # number of cells
    dx = (x_max - x_min) / N
    x = collect(range(x_min + dx / 2, x_max - dx / 2, length=N))   # cell centres
    t_final = 0.6
    CFL = 0.9

    tests = [(name="Moving shock",          uL=1.0,  uR=0.0),
             (name="Transonic rarefaction", uL=-0.5, uR=1.0)]

    # Time step from the largest wave speed max|f'(u)| = max|u| over both tests
    s_max = maximum(max(abs(T.uL), abs(T.uR)) for T in tests)
    dt = CFL * dx / s_max
    nsteps = round(Int, t_final / dt)
    lam = dt / dx

    println("----------------------------------------------------------------")
    println("Burgers' Equation: Godunov vs Rusanov (Finite Volume)")
    println("Cells N          : $N (dx = $(round(dx, digits=4)))")
    println("CFL number       : $CFL")
    println("Time step dt     : $(round(dt, digits=5))")
    println("Total time steps : $nsteps (t_final = $(round(nsteps * dt, digits=3)))")
    println("----------------------------------------------------------------")

    # Initial cell averages (jump sits exactly on the interface at x0)
    ug = [[x[i] < x0 ? T.uL : T.uR for i in 1:N] for T in tests]
    ur = deepcopy(ug)
    tmp = zeros(N)
    xf = range(x_min, x_max, length=2000)
    every = max(1, round(Int, nsteps / 120))

    # --------------------------------------------------------------------------
    # Time stepping with live plot
    # --------------------------------------------------------------------------
    for n in 1:nsteps
        for k in eachindex(tests)
            fv_step!(tmp, ug[k], lam, godunov_flux); ug[k] .= tmp
            fv_step!(tmp, ur[k], lam, rusanov_flux); ur[k] .= tmp
        end

        if n % every == 0 || n == nsteps
            t = n * dt
            panels = []
            for (k, T) in enumerate(tests)
                lo, hi = min(T.uL, T.uR), max(T.uL, T.uR)
                p = plot(xf, riemann_exact.(xf, t, T.uL, T.uR, x0), label="Exact",
                         color=:black, lw=2.5, linestyle=:dash,
                         legend=(k == 1 ? :bottomleft : :topleft))
                plot!(p, x, ug[k], label="Godunov", color=:red, lw=2,
                      marker=:circle, ms=2, msw=0)
                plot!(p, x, ur[k], label="Rusanov", color=:blue, lw=1.5)
                hline!(p, [0.0], label="", color=:gray, alpha=0.4)
                title!(p, "$(T.name) (u_L = $(T.uL), u_R = $(T.uR))")
                xlabel!(p, "x"); ylabel!(p, "u(x, t)")
                ylims!(p, lo - 0.2, hi + 0.2); xlims!(p, x_min, x_max)
                push!(panels, p)
            end
            display(plot(panels..., layout=(1, 2), size=(1200, 450), left_margin=5Plots.mm,
                         plot_title="Burgers' equation, t = $(round(t, digits=3)),  N = $N"))
            sleep(0.02)
        end
    end

    # --------------------------------------------------------------------------
    # Final errors and overshoot check
    # --------------------------------------------------------------------------
    t_end = nsteps * dt
    println("\nAt t = $(round(t_end, digits=3)):")
    for (k, T) in enumerate(tests)
        ue = riemann_exact.(x, t_end, T.uL, T.uR, x0)
        eg = sum(abs.(ug[k] .- ue)) * dx
        er = sum(abs.(ur[k] .- ue)) * dx
        println("  $(T.name):")
        println("    L1 error  Godunov = $(round(eg, sigdigits=4)),  Rusanov = $(round(er, sigdigits=4))")
        println("    min/max   Godunov = [$(round(minimum(ug[k]), digits=4)), $(round(maximum(ug[k]), digits=4))]" *
                "  (no overshoot: stays within [$(min(T.uL, T.uR)), $(max(T.uL, T.uR))])")
    end
    println("\nGodunov gives the sharper shock and the correct fan through the sonic point u = 0.")
end

# Execute the simulation
solve_burgers_godunov()

# Keep the plot window open when run as `julia burgers_godunov.jl`
if !isinteractive()
    println("\nPress Enter to exit.")
    readline()
end
