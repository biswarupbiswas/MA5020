# ====================================================================
# Solving Burgers' Equation with u0(x) = exp(-x^2) up to t <= t_s
# Method of Characteristics with Live Running Plot
# ====================================================================

using Plots

# 1. Define initial condition and its derivative
u0(x) = exp(-x^2)
u0_prime(x) = -2.0 * x * exp(-x^2)

# 2. Compute the exact breaking time: t_s = -1 / min(u0'(x))
function compute_breaking_time()
    x_test = range(-3.0, 3.0, length=1000)
    min_slope = minimum([u0_prime(xi) for xi in x_test])
    t_s = -1.0 / min_slope
    return t_s
end

# 3. For given (x, t), find initial coordinate x0: x0 + u0(x0)*t = x
function find_x0(x, t; tol=1e-7, max_iter=100)
    a = x - 3.0 * (t + 1.0)
    b = x + 3.0 * (t + 1.0)
    g(x0) = x0 + u0(x0) * t - x
    
    for _ in 1:max_iter
        c = (a + b) / 2.0
        if abs(g(c)) < tol || (b - a) / 2.0 < tol
            return c
        end
        if g(a) * g(c) < 0.0
            b = c
        else
            a = c
        end
    end
    return (a + b) / 2.0
end

# 4. Evaluate solution u(x,t) = u0(x0) on a spatial grid
function solve_burgers(x_grid, t)
    u_sol = zeros(length(x_grid))
    for (i, xi) in enumerate(x_grid)
        x0 = find_x0(xi, t)
        u_sol[i] = u0(x0)
    end
    return u_sol
end

# 5. Live Running Plot Simulation
function run_simulation()
    t_s = compute_breaking_time()
    println("Analytical Breaking Time t_s = ", round(t_s, digits=4))
    println("Running live animation of wave steepening up to t <= t_s...")

    x_grid = range(-3.0, 4.0, length=300)
    time_steps = range(0.0, 0.99 * t_s, length=60)

    for t in time_steps
        u_current = solve_burgers(x_grid, t)
        
        p = plot(x_grid, u_current,
                 lw=2.5,
                 color=:darkblue,
                 xlims=(-3.0, 4.0),
                 ylims=(-0.1, 1.2),
                 xlabel="Spatial Coordinate x",
                 ylabel="Solution u(x,t)",
                 title="Burgers Equation Steepening (t = $(round(t, digits=3)), t_s = $(round(t_s, digits=3)))",
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
