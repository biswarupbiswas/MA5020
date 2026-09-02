# ==============================================================================
# Numerical Solution of 1D Linear Advection Equation: u_t + a * u_x = 0
# Method: Beam-Warming Scheme (2nd-Order Upwind Differencing)
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow
# ==============================================================================

using Plots

function solve_advection_beam_warming()
    # -------------------------------------------------------------------------
    # 1. Physical Parameters and Domain Setup
    # -------------------------------------------------------------------------
    a = 1.0               # Wave advection speed (a > 0: right-traveling wave)
    x_min = 0.0           # Left boundary of spatial domain
    x_max = 2.0           # Right boundary of spatial domain
    L = x_max - x_min     # Domain length
    t_final = 2.0         # Final simulation time (one complete period)

    # -------------------------------------------------------------------------
    # 2. Numerical Discretization Parameters
    # -------------------------------------------------------------------------
    N = 200               # Number of spatial grid intervals
    dx = L / N            # Grid spacing
    x = range(x_min, x_max - dx, length=N)  # Node locations: x_j

    nu = 0.8              # Courant number (CFL condition: 0 <= nu <= 2 for stability)
    dt = nu * dx / abs(a) # Time step
    nsteps = round(Int, t_final / dt)
    actual_t_final = nsteps * dt

    println("----------------------------------------------------------------")
    println("Linear Advection Solver: Beam-Warming Scheme (2nd-Order Upwind)")
    println("Grid points N    : $N (dx = $(round(dx, digits=4)))")
    println("Wave speed a     : $a")
    println("CFL Number nu    : $nu")
    println("Time step dt     : $(round(dt, digits=5))")
    println("Total time steps : $nsteps (t_final = $(round(actual_t_final, digits=3)))")
    println("----------------------------------------------------------------")

    # -------------------------------------------------------------------------
    # 3. Initial Condition: u(x, 0) = u0(x)
    # -------------------------------------------------------------------------
    # Combination of a smooth Gaussian wave and a sharp square pulse
    u0(xi) = exp(-100.0 * (xi - 0.5)^2) + (1.1 <= xi <= 1.5 ? 1.0 : 0.0)

    u = [u0(xi) for xi in x]          # Current numerical solution: u^n
    u_new = copy(u)                   # Next time level solution: u^{n+1}

    # -------------------------------------------------------------------------
    # 4. Time Stepping Loop with Live Animation
    # -------------------------------------------------------------------------
    # Beam-Warming update formula for a > 0 (one-sided backward differences):
    #   u_j^{n+1} = u_j^n - 0.5 * nu * (3*u_j^n - 4*u_{j-1}^n + u_{j-2}^n)
    #                     + 0.5 * nu^2 * (u_j^n - 2*u_{j-1}^n + u_{j-2}^n)
    # Note: Leading error is positive dispersion -> leading oscillations!

    anim = @animate for n in 1:nsteps
        t_current = n * dt

        # Spatial update using periodic boundary conditions
        for j in 1:N
            if a > 0
                jm1 = (j == 1) ? N : j - 1      # Periodic index for j-1
                jm2 = (j == 1) ? N - 1 : ((j == 2) ? N : j - 2) # Periodic index for j-2

                u_new[j] = u[j] - 0.5 * nu * (3.0 * u[j] - 4.0 * u[jm1] + u[jm2]) +
                                  0.5 * nu^2 * (u[j] - 2.0 * u[jm1] + u[jm2])
            else
                jp1 = (j == N) ? 1 : j + 1      # Periodic index for j+1
                jp2 = (j == N) ? 2 : ((j == N - 1) ? 1 : j + 2) # Periodic index for j+2

                u_new[j] = u[j] + 0.5 * nu * (3.0 * u[j] - 4.0 * u[jp1] + u[jp2]) +
                                  0.5 * nu^2 * (u[j] - 2.0 * u[jp1] + u[jp2])
            end
        end

        u .= u_new

        # Analytical exact solution: u(x, t) = u0( (x - a*t) mod L )
        u_exact = [u0(mod(xi - a * t_current, L) + x_min) for xi in x]

        # Plot comparison
        plot(x, u_exact, label="Exact Solution", color=:black, lw=2, linestyle=:dash)
        plot!(x, u, label="Beam-Warming (2nd Upwind)", color=:magenta, lw=2)
        xlabel!("Spatial coordinate x")
        ylabel!("u(x, t)")
        title!("1D Advection (Beam-Warming): t = $(round(t_current, digits=2)) (CFL = $nu)")
        ylims!(-0.4, 1.5)
        xlims!(x_min, x_max)
    end every max(1, round(Int, nsteps / 100))

    gif(anim, "advection_beam_warming_simulation.gif", fps=25)
    println("Simulation complete! Animation saved as 'advection_beam_warming_simulation.gif'.")

    # -------------------------------------------------------------------------
    # 5. Final State Comparison Plot
    # -------------------------------------------------------------------------
    u_exact_final = [u0(mod(xi - a * actual_t_final, L) + x_min) for xi in x]
    p_final = plot(x, u_exact_final, label="Exact Solution", color=:black, lw=2, linestyle=:dash)
    plot!(p_final, x, u, label="Beam-Warming Scheme", color=:magenta, lw=2)
    xlabel!(p_final, "Spatial coordinate x")
    ylabel!(p_final, "u(x, t)")
    title!(p_final, "Beam-Warming at t = $(round(actual_t_final, digits=2)) (CFL = $nu, N = $N)")
    ylims!(p_final, -0.4, 1.5)
    savefig(p_final, "advection_beam_warming_final.png")
    println("Final comparison plot saved as 'advection_beam_warming_final.png'.")
end

solve_advection_beam_warming()
