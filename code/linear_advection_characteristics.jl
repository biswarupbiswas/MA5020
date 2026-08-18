# ====================================================================
# Solving Linear Advection Equation: u_t + a*u_x = 0
# Method of Characteristics with Live Running Plot
# ====================================================================

using Plots

# 1. Problem parameters
const a = 1.5                                      # Constant wave speed
u0(x) = exp(-10.0 * (x - 0.5)^2)                   # Gaussian initial profile

# 2. Exact solution along characteristics
function solve_linear_advection(x_grid, t; speed=a)
    u_sol = zeros(length(x_grid))
    for (i, xi) in enumerate(x_grid)
        x0 = xi - speed * t                        # Trace back characteristic origin
        u_sol[i] = u0(x0)                          # Value is constant along ray
    end
    return u_sol
end

# 3. Live Running Plot Simulation
function run_simulation()
    x_grid = range(-1.0, 5.0, length=300)
    t_final = 2.0
    time_steps = range(0.0, t_final, length=60)

    println("Running live animation: Linear Advection (a = $a)...")

    for t in time_steps
        u_current = solve_linear_advection(x_grid, t)
        
        p = plot(x_grid, u_current,
                 lw=2.5,
                 color=:crimson,
                 xlims=(-1.0, 5.0),
                 ylims=(-0.1, 1.2),
                 xlabel="Spatial Coordinate x",
                 ylabel="Solution u(x,t)",
                 title="Linear Advection: u_t + 1.5 u_x = 0 (t = $(round(t, digits=2)))",
                 legend=false,
                 grid=true)
        
        display(p)
        sleep(0.03)  # Pause for animation effect
    end
    
    println("Simulation finished.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_simulation()
end
