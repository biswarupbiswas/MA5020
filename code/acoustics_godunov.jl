# ==============================================================================
# Linear Acoustics: q_t + A q_x = 0,  q = (p, u)
#     A = [u0  K0; 1/rho0  u0],  eigenvalues u0 - c0, u0 + c0
# Method: Godunov (characteristic upwind) finite volume scheme,
#         compared with the Rusanov flux and the exact Riemann solution
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 14)
# ==============================================================================

using Plots, LinearAlgebra

# Background state: change u0 to see the difference between the two fluxes
rho0, c0, u0 = 1.0, 1.0, 0.5
K0, Z0 = rho0 * c0^2, rho0 * c0                 # bulk modulus, impedance

A   = [u0 K0; 1/rho0 u0]
lam = [u0 - c0, u0 + c0]                        # eigenvalues (wave speeds)
R   = [-Z0 Z0; 1.0 1.0]                         # eigenvectors as columns
Rinv = inv(R)
Aplus  = R * Diagonal(max.(lam, 0)) * Rinv      # right-going part of A
Aminus = R * Diagonal(min.(lam, 0)) * Rinv      # left-going part of A
smax = maximum(abs.(lam))                       # spectral radius of A

# Numerical fluxes at an interface with states qL (left) and qR (right)
godunov_flux(qL, qR) = Aplus * qL + Aminus * qR
rusanov_flux(qL, qR) = 0.5 * (A * qL + A * qR) - 0.5 * smax * (qR - qL)

# One finite volume step for all cells; q is a 2 x N matrix
function fv_step(q, dt, dx, flux)
    qg = hcat(q[:, 1], q, q[:, end])            # ghost cells (zero gradient)
    N = size(q, 2)
    F = hcat([flux(qg[:, i], qg[:, i+1]) for i in 1:N+1]...)
    return q .- (dt / dx) .* (F[:, 2:end] .- F[:, 1:end-1])
end

# ------------------------------------------------------------------------------
# Riemann problem: a pressure jump, fluid initially at rest
# ------------------------------------------------------------------------------
qL_data, qR_data = [1.0, 0.0], [0.0, 0.0]

# Exact solution: decompose the jump into the two eigenvectors
alpha = Rinv * (qR_data - qL_data)
qm = qL_data + alpha[1] * R[:, 1]                 # middle state
q_exact(x, t) = x < lam[1] * t ? qL_data : (x < lam[2] * t ? qm : qR_data)
println("Middle state: p_m = ", qm[1], ",  u_m = ", qm[2])

N = 100
dx = 2 / N
x = range(-1 + dx/2, 1 - dx/2, length=N)        # cell centres on [-1, 1]
dt = 0.9 * dx / smax                            # CFL = 0.9
t_final = 0.4

q0 = hcat([xi < 0 ? qL_data : qR_data for xi in x]...)
qG, qRu = copy(q0), copy(q0)
t = 0.0
anim = @animate while t < t_final - 1e-12
    global qG, qRu, t
    dtn = min(dt, t_final - t)
    qG  = fv_step(qG,  dtn, dx, godunov_flux)
    qRu = fv_step(qRu, dtn, dx, rusanov_flux)
    t += dtn

    xe = range(-1, 1, length=1001)
    qe = hcat([q_exact(xi, t) for xi in xe]...)
    p1 = plot(xe, qe[1, :], color=:black, label="Exact", ylims=(-0.1, 1.1),
              xlabel="x", ylabel="p", title="Pressure, t = $(round(t, digits=3))")
    scatter!(p1, x, qG[1, :],  ms=2.5, color=:red,  label="Godunov")
    scatter!(p1, x, qRu[1, :], ms=2.5, color=:blue, label="Rusanov", markershape=:square)
    p2 = plot(xe, qe[2, :], color=:black, label="Exact", ylims=(-0.1, 1.1),
              xlabel="x", ylabel="u", title="Velocity, t = $(round(t, digits=3))")
    scatter!(p2, x, qG[2, :],  ms=2.5, color=:red,  label="Godunov")
    scatter!(p2, x, qRu[2, :], ms=2.5, color=:blue, label="Rusanov", markershape=:square)
    plot(p1, p2, layout=(1, 2), size=(1000, 380))
end
gif(anim, "acoustics_godunov.gif", fps=15)
