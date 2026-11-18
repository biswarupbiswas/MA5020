# ==============================================================================
# Shallow Water Equations: h_t + (hu)_x = 0,  (hu)_t + (hu^2 + g h^2/2)_x = 0
# Problem: Dam break (Riemann problem with water at rest)
# Method: Finite volume scheme with three approximate Riemann solvers
#         (Rusanov, HLL, Roe), compared with the exact solution
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 15)
# ==============================================================================

using Plots

g = 1.0
h_l, u_l = 3.0, 0.0                       # left state  (deep water)
h_r, u_r = 1.0, 0.0                       # right state (shallow water)

# ------------------------------------------------------------------------------
# Approximate Riemann solvers: flux at an interface with states (hL, mL), (hR, mR)
# m = hu is the discharge
# ------------------------------------------------------------------------------
flux(h, m) = [m, m^2 / h + 0.5 * g * h^2]

function rusanov_flux(hL, mL, hR, mR)
    uL, uR = mL / hL, mR / hR
    a = max(abs(uL) + sqrt(g * hL), abs(uR) + sqrt(g * hR))
    return 0.5 * (flux(hL, mL) + flux(hR, mR)) - 0.5 * a * [hR - hL, mR - mL]
end

function hll_flux(hL, mL, hR, mR)
    uL, uR = mL / hL, mR / hR
    cL, cR = sqrt(g * hL), sqrt(g * hR)
    SL = min(uL - cL, uR - cR)            # slowest wave speed
    SR = max(uL + cL, uR + cR)            # fastest wave speed
    FL, FR = flux(hL, mL), flux(hR, mR)
    SL >= 0 && return FL
    SR <= 0 && return FR
    return (SR * FL - SL * FR + SL * SR * [hR - hL, mR - mL]) / (SR - SL)
end

function roe_flux(hL, mL, hR, mR)
    uL, uR = mL / hL, mR / hR
    uh = (sqrt(hL) * uL + sqrt(hR) * uR) / (sqrt(hL) + sqrt(hR))   # Roe averages
    ch = sqrt(0.5 * g * (hL + hR))
    dh, dm = hR - hL, mR - mL
    a1 = 0.5 * (dh - (dm - uh * dh) / ch)                          # wave strengths
    a2 = 0.5 * (dh + (dm - uh * dh) / ch)
    l1, l2 = uh - ch, uh + ch                                      # wave speeds
    diss = abs(l1) * a1 * [1, l1] + abs(l2) * a2 * [1, l2]         # |A_hat| (q_R - q_L)
    return 0.5 * (flux(hL, mL) + flux(hR, mR)) - 0.5 * diss
end

# ------------------------------------------------------------------------------
# Finite volume solver
# ------------------------------------------------------------------------------
function solve(numflux, N, t_final; L=5.0, cfl=0.9)
    dx = 2L / N
    x = range(-L + dx/2, L - dx/2, length=N)       # cell centres on [-L, L]
    h = [xi < 0 ? h_l : h_r for xi in x]
    m = [xi < 0 ? h_l * u_l : h_r * u_r for xi in x]
    t = 0.0
    while t < t_final - 1e-12
        smax = maximum(abs.(m ./ h) .+ sqrt.(g .* h))
        dt = min(cfl * dx / smax, t_final - t)      # CFL condition, every step
        hg = [h[1]; h; h[end]]                      # ghost cells (zero gradient)
        mg = [m[1]; m; m[end]]
        F = [numflux(hg[i], mg[i], hg[i+1], mg[i+1]) for i in 1:N+1]
        h = h .- (dt / dx) .* ([Fi[1] for Fi in F[2:end]] .- [Fi[1] for Fi in F[1:end-1]])
        m = m .- (dt / dx) .* ([Fi[2] for Fi in F[2:end]] .- [Fi[2] for Fi in F[1:end-1]])
        t += dt
    end
    return x, h, m ./ h
end

# ------------------------------------------------------------------------------
# Exact solution (reference only; the exact Riemann solver is NOT covered in
# the lecture). f_K gives the change in u across a shock or rarefaction.
# ------------------------------------------------------------------------------
fK(h, hK)  = h <= hK ? 2 * (sqrt(g * h) - sqrt(g * hK)) : (h - hK) * sqrt(0.5 * g * (1/h + 1/hK))
dfK(h, hK) = h <= hK ? sqrt(g / h) :
             sqrt(0.5 * g * (1/h + 1/hK)) - (h - hK) * 0.25 * g / (h^2 * sqrt(0.5 * g * (1/h + 1/hK)))

function middle_state(hl, ul, hr, ur)
    h = 0.5 * (hl + hr)
    for k in 1:50                                   # Newton iteration
        dh = (fK(h, hl) + fK(h, hr) + ur - ul) / (dfK(h, hl) + dfK(h, hr))
        h = max(h - dh, 1e-12)
        abs(dh) < 1e-14 && break
    end
    return h, 0.5 * (ul + ur) + 0.5 * (fK(h, hr) - fK(h, hl))
end

h_m, u_m = middle_state(h_l, u_l, h_r, u_r)
c_l, c_r, c_m = sqrt(g * h_l), sqrt(g * h_r), sqrt(g * h_m)
println("Exact middle state: h_m = ", h_m, ",  u_m = ", u_m)

function exact(s)                                   # solution at x/t = s
    if s <= u_m                                     # 1-wave side
        if h_m > h_l                                # 1-shock
            return s < (h_m * u_m - h_l * u_l) / (h_m - h_l) ? (h_l, u_l) : (h_m, u_m)
        elseif s < u_l - c_l
            return (h_l, u_l)
        elseif s > u_m - c_m
            return (h_m, u_m)
        else                                        # inside 1-rarefaction
            c = (u_l + 2c_l - s) / 3
            return (c^2 / g, s + c)
        end
    else                                            # 2-wave side
        if h_m > h_r                                # 2-shock
            return s < (h_m * u_m - h_r * u_r) / (h_m - h_r) ? (h_m, u_m) : (h_r, u_r)
        elseif s > u_r + c_r
            return (h_r, u_r)
        elseif s < u_m + c_m
            return (h_m, u_m)
        else                                        # inside 2-rarefaction
            c = (s - u_r + 2c_r) / 3
            return (c^2 / g, s - c)
        end
    end
end

# ------------------------------------------------------------------------------
# Run the three solvers and compare
# ------------------------------------------------------------------------------
N, t_final = 100, 1.0
solvers = [("Rusanov", rusanov_flux, :blue), ("HLL", hll_flux, :green), ("Roe", roe_flux, :red)]

xe = range(-5, 5, length=2001)
qe = [exact(xi / t_final) for xi in xe]
p1 = plot(xe, first.(qe), color=:black, lw=2, label="Exact", xlabel="x", ylabel="h",
          title="Depth, t = $t_final, N = $N", ylims=(0.8, 3.2))
p2 = plot(xe, last.(qe), color=:black, lw=2, label="Exact", xlabel="x", ylabel="u",
          title="Velocity, t = $t_final, N = $N", ylims=(-0.05, 0.85), legend=:topleft)
for (name, F, col) in solvers
    x, h, u = solve(F, N, t_final)
    scatter!(p1, x, h, ms=2.5, color=col, label=name)
    scatter!(p2, x, u, ms=2.5, color=col, label=name)
end
display(plot(p1, p2, layout=(1, 2), size=(1000, 380)))

# L1 error in h for each solver (exact cell averages on a fine sub-grid)
println("\nL1 error in h at t = $t_final:")
for (name, F, col) in solvers
    errs = Float64[]
    for N in (100, 200, 400)
        x, h, u = solve(F, N, t_final)
        dx = 10 / N
        hbar = [sum(first(exact((xi - dx/2 + (k - 0.5) * dx / 40) / t_final)) for k in 1:40) / 40 for xi in x]
        push!(errs, dx * sum(abs.(h .- hbar)))
    end
    println(rpad(name, 8), ": ", join(round.(errs, sigdigits=3), ",  "))
end
