# ==============================================================================
# 2D Euler equations: q = (rho, rho u, rho v, E),  p = (gamma-1)(E - rho (u^2+v^2)/2)
# Problems: (1) a 2D Riemann problem on [0,1]^2, (2) the isentropic vortex (smooth, exact solution)
# Method: finite volume, Rusanov and HLLC fluxes, dimensional splitting (x- and y-sweeps)
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 17)
#
# The file is split into cells (lines starting with ##). In VS Code run one
# cell at a time with Shift+Enter, or run the whole file with  julia euler2d_solvers.jl
# ==============================================================================

## 0. Packages and the gas constant
using Plots

const γ = 1.4

## 1. Conversions between primitive (rho, u, v, p) and conserved variables
cons(ρ, u, v, p) = [ρ, ρ * u, ρ * v, p / (γ - 1) + 0.5 * ρ * (u^2 + v^2)]
function prim(q)
    ρ, u, v = q[1], q[2] / q[1], q[3] / q[1]
    return ρ, u, v, (γ - 1) * (q[4] - 0.5 * ρ * (u^2 + v^2))
end
function flux_x(q)
    ρ, u, v, p = prim(q)
    return [ρ * u, ρ * u^2 + p, ρ * u * v, u * (q[4] + p)]
end
soundspeed(ρ, p) = sqrt(γ * p / ρ)

## 2a. Rusanov flux in the x-direction
function rusanov_x(qL, qR)
    ρL, uL, vL, pL = prim(qL); ρR, uR, vR, pR = prim(qR)
    a = max(abs(uL) + soundspeed(ρL, pL), abs(uR) + soundspeed(ρR, pR))
    return 0.5 * (flux_x(qL) + flux_x(qR)) - 0.5 * a * (qR - qL)
end

## 2b. HLLC flux in the x-direction (v is carried across the contact unchanged)
function hllc_x(qL, qR)
    ρL, uL, vL, pL = prim(qL); ρR, uR, vR, pR = prim(qR)
    cL, cR = soundspeed(ρL, pL), soundspeed(ρR, pR)
    SL, SR = min(uL - cL, uR - cR), max(uL + cL, uR + cR)
    FL, FR = flux_x(qL), flux_x(qR)
    SL >= 0 && return FL
    SR <= 0 && return FR
    Ss = (pR - pL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)) / (ρL * (SL - uL) - ρR * (SR - uR))
    ps = pL + ρL * (SL - uL) * (Ss - uL)
    function qstar(ρ, u, v, p, E, S)
        ρs = ρ * (S - u) / (S - Ss)
        Es = (E * (S - u) - p * u + ps * Ss) / (S - Ss)
        return [ρs, ρs * Ss, ρs * v, Es]          # tangential velocity v unchanged
    end
    Ss >= 0 && return FL + SL * (qstar(ρL, uL, vL, pL, qL[4], SL) - qL)
    return FR + SR * (qstar(ρR, uR, vR, pR, qR[4], SR) - qR)
end

## 3. y-flux from the x-flux: swap the two momentum components
swap(q) = [q[1], q[3], q[2], q[4]]
flux_y(numflux, qB, qT) = swap(numflux(swap(qB), swap(qT)))

# Quick check: Sod in the y-direction gives the 1D fluxes of Lecture 16
for (name, F) in [("Rusanov", rusanov_x), ("HLLC", hllc_x)]
    println(rpad(name, 8), " G at a Sod interface in y: ",
            round.(flux_y(F, cons(1.0, 0, 0, 1.0), cons(0.125, 0, 0, 0.1)), digits=5))
end

## 4. Ghost cells: transmissive (copy), wall (mirror the normal velocity); periodic is handled in the sweeps
function ghost(q, bc, normal)          # normal = 2 for x-walls (u), 3 for y-walls (v)
    g = copy(q)
    bc == :wall && (g[normal] = -g[normal])
    return g
end

## 5. One x-sweep and one y-sweep over the whole grid
#    Q is an Nx × Ny matrix of 4-vectors
function x_sweep(numflux, Q, dt, dx, bc)
    Nx, Ny = size(Q)
    Qn = similar(Q)
    for j in 1:Ny
        row = bc == :periodic ? [[Q[Nx, j]]; Q[:, j]; [Q[1, j]]] :
                               [[ghost(Q[1, j], bc, 2)]; Q[:, j]; [ghost(Q[Nx, j], bc, 2)]]
        F = [numflux(row[i], row[i+1]) for i in 1:Nx+1]
        for i in 1:Nx
            Qn[i, j] = Q[i, j] - (dt / dx) * (F[i+1] - F[i])
        end
    end
    return Qn
end

function y_sweep(numflux, Q, dt, dy, bc)
    Nx, Ny = size(Q)
    Qn = similar(Q)
    for i in 1:Nx
        col = bc == :periodic ? [[Q[i, Ny]]; Q[i, :]; [Q[i, 1]]] :
                               [[ghost(Q[i, 1], bc, 3)]; Q[i, :]; [ghost(Q[i, Ny], bc, 3)]]
        G = [flux_y(numflux, col[j], col[j+1]) for j in 1:Ny+1]
        for j in 1:Ny
            Qn[i, j] = Q[i, j] - (dt / dy) * (G[j+1] - G[j])
        end
    end
    return Qn
end

## 6. Time loop with dimensional splitting (order of sweeps alternates)
function solve2d(numflux, Q, dx, dy, t_final; cfl=0.9, bc=:transmissive)
    t, n = 0.0, 0
    while t < t_final - 1e-12
        sx = maximum(q -> (W = prim(q); abs(W[2]) + soundspeed(W[1], W[4])), Q)
        sy = maximum(q -> (W = prim(q); abs(W[3]) + soundspeed(W[1], W[4])), Q)
        dt = min(cfl * min(dx / sx, dy / sy), t_final - t)
        if iseven(n)
            Q = y_sweep(numflux, x_sweep(numflux, Q, dt, dx, bc), dt, dy, bc)
        else
            Q = x_sweep(numflux, y_sweep(numflux, Q, dt, dy, bc), dt, dx, bc)
        end
        t += dt; n += 1
    end
    println(rpad(string(numflux), 10), " reached t = $t in $n steps")
    return Q
end

## 7. Test 1: a 2D Riemann problem on [0,1]^2 (four constant quadrants), t = 0.3
#    (rho, u, v, p) in each quadrant
TR = (1.5,    0.0,   0.0,   1.5)      # top right
TL = (0.5323, 1.206, 0.0,   0.3)      # top left
BL = (0.138,  1.206, 1.206, 0.029)    # bottom left
BR = (0.5323, 0.0,   1.206, 0.3)      # bottom right

N = 200                               # 400 gives sharper pictures but takes longer
dx = 1 / N
x = range(dx/2, 1 - dx/2, length=N)
function quadrant(xi, yj)
    xi > 0.5 && yj > 0.5 && return cons(TR...)
    xi <= 0.5 && yj > 0.5 && return cons(TL...)
    xi <= 0.5 && return cons(BL...)
    return cons(BR...)
end
Q0 = [quadrant(xi, yj) for xi in x, yj in x]

## 8. Run both solvers and plot the density
plots_rp = []
for (name, F) in [("Rusanov", rusanov_x), ("HLLC", hllc_x)]
    Q = solve2d(F, Q0, dx, dx, 0.3)
    push!(plots_rp, heatmap(x, x, [q[1] for q in Q]', aspect_ratio=1, c=:viridis,
                            xlabel="x", ylabel="y", title="$name, density, t = 0.3"))
end
display(plot(plots_rp..., layout=(1, 2), size=(1100, 450)))

## 9. Test 2: isentropic vortex on the periodic square [0,10]^2, mean flow (1,1)
const ε = 5.0
function vortex(xi, yj, t; L=10.0)
    xb = mod(xi - 5 - t + L/2, L) - L/2          # distance from the moving centre
    yb = mod(yj - 5 - t + L/2, L) - L/2
    r2 = xb^2 + yb^2
    T = 1 - (γ - 1) * ε^2 / (8γ * π^2) * exp(1 - r2)
    ρ = T^(1 / (γ - 1))
    u = 1 - ε / (2π) * exp((1 - r2) / 2) * yb
    v = 1 + ε / (2π) * exp((1 - r2) / 2) * xb
    return ρ, u, v, ρ^γ
end

function run_vortex(numflux, N; t_final=2.0)
    dx = 10 / N
    x = range(dx/2, 10 - dx/2, length=N)
    Q = solve2d(numflux, [cons(vortex(xi, yj, 0.0)...) for xi in x, yj in x], dx, dx, t_final;
                bc=:periodic)
    ρ = [q[1] for q in Q]
    ρex = [vortex(xi, yj, t_final)[1] for xi in x, yj in x]
    println("  N = $N:  L1 error in density = ", round(sum(abs.(ρ - ρex)) * dx^2, digits=4))
    return x, ρ, ρex
end

## 10. Vortex: density along a line through the centre at t = 2 (centre at (7,7))
xv, ρ_rus, ρex = run_vortex(rusanov_x, 100)
xv, ρ_hllc, _  = run_vortex(hllc_x, 100)
j = argmin(abs.(xv .- 7))
pv = plot(xv, ρex[:, j], color=:black, lw=2, label="exact", xlabel="x", ylabel="ρ",
          title="vortex: density along y = $(round(xv[j], digits=2)), t = 2")
scatter!(pv, xv, ρ_rus[:, j], ms=3, color=:blue, label="Rusanov")
scatter!(pv, xv, ρ_hllc[:, j], ms=3, color=:red, label="HLLC")
display(pv)

## 11. Vortex: first order, so halving dx roughly halves the error (HLLC)
for N in (50, 100, 200)
    run_vortex(hllc_x, N)
end
