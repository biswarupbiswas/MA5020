# ⚡ Julia Quickstart Guide for Computational Methods

This quick primer covers the essential Julia syntax and patterns needed for implementing numerical methods in this course. It is designed to get you up and running in **under 10 minutes**.

---

## 1. Why Julia for Numerical Methods?
* **Fast like C/Fortran:** Just-In-Time (JIT) compiled to native machine code via LLVM.
* **Readable like Python/MATLAB:** Clean, mathematical syntax.
* **1-Based Indexing:** Array indices start at `1` (matching standard mathematical and textbook notation: $x_1, x_2, \dots, x_N$).
* **Native Linear Algebra & Vectorization:** First-class support for vectors, matrices, and broadcasting.

---

## 2. Installation & Running Code

### Installation
Download and install Julia from the official site: **[https://julialang.org/downloads/](https://julialang.org/downloads/)**

### Running Scripts
Save your code in a file, e.g. `solver.jl`, and run it from your terminal:
```bash
julia solver.jl
```

### Interactive REPL
Type `julia` in your terminal to open the interactive Read-Eval-Print-Loop:
```julia
julia> 2 + 2
4
```
*(To exit the REPL, type `exit()` or press `Ctrl+D`)*.

---

## 3. Core Syntax & Numerical Essentials

### Variables & Math Operations
```julia
a = 1.5           # Float64
N = 100           # Int64
pi_val = pi       # Built-in mathematical constants (pi, ℯ)

# Standard operators: +, -, *, /, ^ (power)
c = 3.0 * 10^8
```

### Defining Functions
```julia
# Compact one-line mathematical functions:
f(x) = x^2 - 4.0
u0(x) = exp(-x^2)

# Multi-line functions:
function wave_speed(u)
    if u > 0.0
        return u
    else
        return 0.0
    end
end
```

---

## 4. Grids, Arrays & Vectorization

### Spatial Discretization & Ranges
```julia
# Create an evenly spaced grid of N points from x_min to x_max
x = range(0.0, 1.0, length=100)
dx = x[2] - x[1]                # Grid spacing Δx

# Preallocate vectors (zeros, ones)
u = zeros(100)                  # Vector of 100 zeros
F = ones(100)                   # Vector of 100 ones
```

### 1-Based Indexing & Slicing
```julia
u[1]            # First element
u[end]          # Last element
u[end-1]        # Second to last element
u[2:end-1]      # Interior slice (from index 2 to N-1)
```

### The Dot Operator (Broadcasting / Element-wise Math)
In Julia, adding a dot `.` before an operator or function applies it **element-by-element** to an array without writing a loop:
```julia
x = range(0.0, 2*pi, length=50)

# Evaluate function at every grid point:
u = sin.(x)                     # Element-wise sine
u_sq = u.^2                     # Element-wise square
g = exp.(-x.^2)                 # Element-wise exponential
```

### Array Comprehensions
```julia
u_init = [u0(xi) for xi in x]
```

---

## 5. Control Flow (Loops & Conditionals)

### `for` Loops
```julia
N = 100
for i in 1:N
    # do something with index i
end

# Stepping in reverse or with stride:
for i in 1:2:N        # Stride of 2
end
```

### `while` Loops (Time Marching)
```julia
t = 0.0
t_final = 1.5
dt = 0.01

while t < t_final
    # Time integration update here...
    t += dt
end
```

### Conditional Statements
```julia
if x > 0.0
    u = 1.0
else
    u = 0.0
end
```

---

## 6. Common Built-in Functions for Numerics

| Task | Julia Syntax | Example |
| :--- | :--- | :--- |
| Minimum value | `minimum(u)` | `min_val = minimum(u)` |
| Maximum value | `maximum(u)` | `max_val = maximum(abs.(u))` |
| Array length | `length(u)` | `N = length(x)` |
| Round numbers | `round(val, digits=4)` | `round(pi, digits=3) # 3.142` |
| Absolute value | `abs(x)` or `abs.(u)` | `abs.([-1.0, 2.0])` |

---

## 7. Numerical Methods Examples

### Example 1: Root Finding via the Bisection Method
Finding the root of a non-linear equation $g(x) = 0$ on an interval $[a, b]$:

```julia
# Define function: g(x) = x^3 - x - 2 = 0
g(x) = x^3 - x - 2.0

function bisection(g, a, b; tol=1e-6, max_iter=50)
    for iter in 1:max_iter
        c = (a + b) / 2.0
        if abs(g(c)) < tol || (b - a) / 2.0 < tol
            println("Converged: root = ", round(c, digits=5), " in $iter iterations")
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

root = bisection(g, 1.0, 2.0)
```

### Example 2: Solving an Initial Value Problem (Forward Euler Method)
Solving a first-order decay equation $\frac{dy}{dt} = -2y$ with $y(0) = 1.0$:

```julia
function solve_decay_ode()
    dt = 0.05
    t_grid = 0.0:dt:2.0
    N = length(t_grid)
    
    y = zeros(N)
    y[1] = 1.0  # Initial condition y(0) = 1.0
    
    for n in 1:N-1
        # Forward Euler update: y_{n+1} = y_n + dt * f(t_n, y_n)
        y[n+1] = y[n] + dt * (-2.0 * y[n])
    end
    
    println("At t = 2.0: Numerical y = $(round(y[end], digits=4)), Exact y = $(round(exp(-4.0), digits=4))")
    return t_grid, y
end

t_vals, y_vals = solve_decay_ode()
```

---

## 💡 Quick Tips
* **No `import numpy` needed:** Vectors, matrices, and linear algebra are built directly into base Julia.
* **Copying Arrays:** Use `u_new = copy(u)` rather than `u_new = u` (assignment alone creates a reference, not a copy).
* **Performance:** Avoid global variables inside performance-critical loops; place your simulation code inside a `function`.
