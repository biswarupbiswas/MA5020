# ==============================================================================
# 1D Euler equations: symbolic eigenvalues and eigenvectors of the flux Jacobian
#   q = (rho, rho u, E),  f(q) = (rho u, rho u^2 + p, u (E + p)),
#   p = (gamma - 1) (E - rho u^2 / 2),  c^2 = gamma p / rho,  H = (E + p) / rho
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 16)
#
# The file is split into cells (lines starting with ##). In VS Code run one
# cell at a time with Shift+Enter, or run the whole file with  julia euler_eigen_symbolic.jl
# ==============================================================================

## 0. Packages (SymPy.jl is installed automatically the first time)
import Pkg
if Base.find_package("SymPy") === nothing
    Pkg.add("SymPy")
end
using SymPy, LinearAlgebra

@syms q1::positive q2::real q3::positive          # conserved variables
@syms ρ::positive u::real p::positive c::positive H::positive γ::positive λ

## 1. Flux in conserved variables and its Jacobian A = df/dq
P = (γ - 1) * (q3 - q2^2 / (2q1))                   # pressure in terms of q
f = [q2, q2^2 / q1 + P, q2 / q1 * (q3 + P)]
q = [q1, q2, q3]
A = [diff(f[i], q[j]) for i in 1:3, j in 1:3]

println("Flux Jacobian A in conserved variables:")
display(simplify.(A))

## 2. The same Jacobian written with u and the enthalpy H
#    q1 = rho, q2 = rho u, q3 = E = rho H - p,  p = rho (γ-1)/γ (H - u^2/2)
p_of_H = ρ * (γ - 1) / γ * (H - u^2 / 2)
to_uH = Dict(q1 => ρ, q2 => ρ * u, q3 => ρ * H - p_of_H)
A_H = simplify.(expand.(subs.(A, Ref(to_uH))))

println("Flux Jacobian A in terms of u and H:")
display(A_H)

## 3. Eigenvalues: factor the characteristic polynomial det(A - λI)
#    Use H = c^2/(γ - 1) + u^2/2 so that the sound speed appears.
A_c = simplify.(subs.(A_H, H => c^2 / (γ - 1) + u^2 / 2))
charpoly = simplify(det(A_c - λ * I(3)))

println("det(A - λI) = ", factor(charpoly))
println("eigenvalues  = ", solve(charpoly, λ))

## 4. Eigenvectors in primitive variables w = (rho, u, p)
#    w_t + B w_x = 0; B has the same eigenvalues as A.
B = [u ρ 0; 0 u 1/ρ; 0 ρ*c^2 u]

println("Primitive matrix B:")
display(B)
for (val, mult, vecs) in B.eigenvects()
    println("λ = ", val, "   r_w = ", simplify.(vec(vecs[1])))
end

## 5. Eigenvectors of A (conserved variables), checked by substitution
#    r1 = (1, u - c, H - u c),  r2 = (1, u, u^2/2),  r3 = (1, u + c, H + u c)
r    = [[1, u - c, H - u * c],          # 1-wave, λ = u - c
        [1, u,     u^2 / 2],            # 2-wave (contact), λ = u
        [1, u + c, H + u * c]]          # 3-wave, λ = u + c
lams = [u - c, u, u + c]

println("Check A r_k - λ_k r_k = 0:")
for k in 1:3
    res = simplify.(subs.(A_H * r[k] - lams[k] * r[k], H => c^2 / (γ - 1) + u^2 / 2))
    println("k = $k   λ = ", lams[k], "   r = ", r[k], "   residual = ", res)
end
