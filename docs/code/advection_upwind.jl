# ==============================================================================
# Numerical Solution of 1D Linear Advection Equation: u_t + a * u_x = 0
# Method: First-Order Upwind (FOU) Finite Difference Scheme
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow
# ==============================================================================

using Plots

function solve_advection_upwind()
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
    N = 200               # Number of spatial grid intervals (cells)
    dx = L / N            # Grid spacing
    x = range(x_min, x_max - dx, length=N)  # Node locations: x_j

    nu = 0.8              # Courant number (CFL condition: 0 < nu <= 1 for stability)
    dt = nu * dx / abs(a) # Time step determined by CFL criterion
    nsteps = round(Int, t_final / dt)
    actual_t_final = nsteps * dt

    println("----------------------------------------------------------------")
    println("Linear Advection Solver: First-Order Upwind Scheme")
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
    # Upwind formula for a > 0 (backward difference):
    #   u_j^{n+1} = u_j^n - nu * (u_j^n - u_{j-1}^n)
    #
    # Upwind formula for a < 0 (forward difference):
    #   u_j^{n+1} = u_j^n - nu * (u_{j+1}^n - u_j^n)

    anim = @animate for n in 1:nsteps
        t_current = n * dt

        # Spatial update using periodic boundary conditions
        for j in 1:N
            if a > 0
                jm = (j == 1) ? N : j - 1  # Periodic wrap-around at left boundary
                u_new[j] = u[j] - nu * (u[j] - u[jm])
            else
                jp = (j == N) ? 1 : j + 1  # Periodic wrap-around at right boundary
                u_new[j] = u[j] + nu * (u[jp] - u[j])
            end
        end

        u .= u_new  # Update state for next step

        # Analytical exact solution: u(x, t) = u0( (x - a*t) mod L )
        u_exact = [u0(mod(xi - a * t_current, L) + x_min) for xi in x]

        # Plot comparison
        plot(x, u_exact, label="Exact Solution", color=:black, lw=2, linestyle=:dash)
        plot!(x, u, label="Upwind Scheme", color=:blue, lw=2)
        xlabel!("Spatial coordinate x")
        ylabel!("u(x, t)")
        title!("1D Advection (Upwind Scheme): t = $(round(t_current, digits=2)) (CFL = $nu)")
        ylims!(-0.2, 1.4)
        xlims!(x_min, x_max)
    end every max(1, round(Int, nsteps / 100))

    # Save animation as GIF
    gif(anim, "advection_upwind_simulation.gif", fps=25)
    println("Simulation complete! Animation saved as 'advection_upwind_simulation.gif'.")

    # -------------------------------------------------------------------------
    # 5. Final State Comparison Plot
    # -------------------------------------------------------------------------
    u_exact_final = [u0(mod(xi - a * actual_t_final, L) + x_min) for xi in x]
    p_final = plot(x, u_exact_final, label="Exact Solution", color=:black, lw=2, linestyle=:dash)
    plot!(p_final, x, u, label="Upwind Scheme (FOU)", color=:blue, lw=2)
    xlabel!(p_final, "Spatial coordinate x")
    ylabel!(p_final, "u(x, t)")
    title!(p_final, "Advection at t = $(round(actual_t_final, digits=2)) (CFL = $nu, N = $N)")
    ylims!(p_final, -0.2, 1.4)
    savefig(p_final, "advection_upwind_final.png")
    println("Final comparison plot saved as 'advection_upwind_final.png'.")
end

# Execute the simulation
solve_advection_upwind()
