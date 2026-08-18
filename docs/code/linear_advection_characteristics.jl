# ====================================================================
# Solving Linear Advection Equation: u_t + a*u_x = 0
# Method of Characteristics Solver: u(x, t) = u0(x - a*t)
# ====================================================================

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

# Main execution:
function main()
    println("Linear Advection Wave Speed a = ", a)
    
    # Test evaluation at x = [0.5, 1.25, 2.0] at different times
    x_test = [0.5, 1.25, 2.0]
    
    for t_val in [0.0, 0.5, 1.0]
        u_vals = solve_linear_advection(x_test, t_val)
        println("t = $(round(t_val, digits=2)): u(0.5) = $(round(u_vals[1], digits=4)), u(1.25) = $(round(u_vals[2], digits=4)), u(2.0) = $(round(u_vals[3], digits=4))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
