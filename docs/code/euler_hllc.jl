# ==============================================================================
# 1D Euler equations: q = (rho, rho u, E),  p = (gamma - 1)(E - rho u^2 / 2)
# Problem: Sod shock tube; optional solid wall at x = 1
# Method: Finite volume scheme with Rusanov, HLL and HLLC fluxes
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 16)
#
# The file is split into cells (lines starting with ##). In VS Code run one
# cell at a time with Shift+Enter, or run the whole file with  julia euler_hllc.jl
# ==============================================================================

## 0. Packages and the gas constant
using Plots

const γ = 1.4

## 1. Conversions between primitive (rho, u, p) and conserved (rho, rho u, E) variables
cons(ρ, u, p) = [ρ, ρ * u, p / (γ - 1) + 0.5 * ρ * u^2]
function prim(q)
    ρ, u = q[1], q[2] / q[1]
    return ρ, u, (γ - 1) * (q[3] - 0.5 * ρ * u^2)
end
function flux(q)
    ρ, u, p = prim(q)
    return [ρ * u, ρ * u^2 + p, u * (q[3] + p)]
end
soundspeed(ρ, p) = sqrt(γ * p / ρ)

## 2a. Rusanov flux
function rusanov_flux(qL, qR)
    ρL, uL, pL = prim(qL); ρR, uR, pR = prim(qR)
    a = max(abs(uL) + soundspeed(ρL, pL), abs(uR) + soundspeed(ρR, pR))
    return 0.5 * (flux(qL) + flux(qR)) - 0.5 * a * (qR - qL)
end

## 2b. HLL flux (wave speed estimates S_L, S_R)
function wave_speeds(qL, qR)
    ρL, uL, pL = prim(qL); ρR, uR, pR = prim(qR)
    cL, cR = soundspeed(ρL, pL), soundspeed(ρR, pR)
    return min(uL - cL, uR - cR), max(uL + cL, uR + cR)
end

function hll_flux(qL, qR)
    SL, SR = wave_speeds(qL, qR)
    FL, FR = flux(qL), flux(qR)
    SL >= 0 && return FL
    SR <= 0 && return FR
    return (SR * FL - SL * FR + SL * SR * (qR - qL)) / (SR - SL)
end

## 2c. HLLC flux (restores the contact wave)
function hllc_flux(qL, qR)
    ρL, uL, pL = prim(qL); ρR, uR, pR = prim(qR)
    SL, SR = wave_speeds(qL, qR)
    FL, FR = flux(qL), flux(qR)
    SL >= 0 && return FL
    SR <= 0 && return FR
    # speed of the contact and pressure in the star region
    Ss = (pR - pL + ρL * uL * (SL - uL) - ρR * uR * (SR - uR)) / (ρL * (SL - uL) - ρR * (SR - uR))
    ps = pL + ρL * (SL - uL) * (Ss - uL)
    # star state on side K (K = L or R)
    function qstar(ρ, u, p, E, S)
        ρs = ρ * (S - u) / (S - Ss)
        Es = (E * (S - u) - p * u + ps * Ss) / (S - Ss)
        return [ρs, ρs * Ss, Es]
    end
    if Ss >= 0
        return FL + SL * (qstar(ρL, uL, pL, qL[3], SL) - qL)
    else
        return FR + SR * (qstar(ρR, uR, pR, qR[3], SR) - qR)
    end
end

## 3. Finite volume solver on [0, 1] with ghost cells
#   right_bc = :transmissive (copy the last cell) or :wall (mirror the velocity)
function solve(numflux, WL, WR, N, t_final; x0=0.5, cfl=0.9, right_bc=:transmissive)
    dx = 1 / N
    x = range(dx/2, 1 - dx/2, length=N)
    Q = [xi < x0 ? cons(WL...) : cons(WR...) for xi in x]
    t = 0.0
    while t < t_final - 1e-12
        smax = maximum(q -> (W = prim(q); abs(W[2]) + soundspeed(W[1], W[3])), Q)
        dt = min(cfl * dx / smax, t_final - t)
        gL = copy(Q[1])                            # left: transmissive
        gR = copy(Q[end])
        right_bc == :wall && (gR[2] = -gR[2])      # wall: u -> -u
        Qg = [[gL]; Q; [gR]]
        F = [numflux(Qg[i], Qg[i+1]) for i in 1:N+1]
        Q = [Q[i] - (dt / dx) * (F[i+1] - F[i]) for i in 1:N]
        t += dt
    end
    W = prim.(Q)
    return x, getindex.(W, 1), getindex.(W, 2), getindex.(W, 3)
end

## 4. Sod shock tube: initial data
WL, WR = (1.0, 0.0, 1.0), (0.125, 0.0, 0.1)     # (rho, u, p)
N, t_final = 100, 0.2

## 5. Flux at the diaphragm at t = 0 (the hand computation in the notes)
for (name, F) in [("Rusanov", rusanov_flux), ("HLL", hll_flux), ("HLLC", hllc_flux)]
    println(rpad(name, 8), " flux at the diaphragm: ", round.(F(cons(WL...), cons(WR...)), digits=5))
end

## 6. Sod shock tube: compare the three solvers
p1 = plot(xlabel="x", ylabel="ρ", title="Density, t = $t_final, N = $N")
p2 = plot(xlabel="x", ylabel="u", title="Velocity", legend=:topleft)
p3 = plot(xlabel="x", ylabel="p", title="Pressure")
for (name, F, col) in [("Rusanov", rusanov_flux, :blue), ("HLL", hll_flux, :green), ("HLLC", hllc_flux, :red)]
    x, ρ, u, p = solve(F, WL, WR, N, t_final)
    scatter!(p1, x, ρ, ms=2.5, color=col, label=name)
    scatter!(p2, x, u, ms=2.5, color=col, label=name)
    scatter!(p3, x, p, ms=2.5, color=col, label=name)
end
display(plot(p1, p2, p3, layout=(1, 3), size=(1200, 360)))

## 7. Solid wall at x = 1: the shock reflects
anim = @animate for tf in 0.02:0.02:0.4
    x, ρ, u, p = solve(hllc_flux, WL, WR, 200, tf; right_bc=:wall)
    plot(x, ρ, color=:red, lw=2, label="HLLC, N = 200", ylims=(0, 1.1),
         xlabel="x", ylabel="ρ", title="Shock reflecting from a wall at x = 1, t = $tf")
    vline!([1.0], color=:black, lw=4, label="wall")
end
gif(anim, "euler_wall.gif", fps=6)
