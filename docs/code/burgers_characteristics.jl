# ====================================================================
# Solving Burgers' Equation with u0(x) = exp(-x^2) up to t <= t_s
# Method of Characteristics Solver
# ====================================================================

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

# Main execution:
function main()
    t_s = compute_breaking_time()
    println("Analytical Breaking Time t_s = ", round(t_s, digits=4))

    # Evaluate profile at t = 0, t = 0.5*t_s, and t = 0.99*t_s
    for t_val in [0.0, 0.5 * t_s, 0.99 * t_s]
        u_vals = solve_burgers([0.0, 1.0, 2.0], t_val)
        println("t = $(round(t_val, digits=3)): u(0.0) = $(round(u_vals[1], digits=4)), u(1.0) = $(round(u_vals[2], digits=4)), u(2.0) = $(round(u_vals[3], digits=4))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
