# ==============================================================================
# 2D linear advection: q_t + a q_x + b q_y = 0 on the periodic unit square
# Problem: a square moving diagonally with velocity (a, b) = (1, 1)
# Method: donor cell (upwind) finite volume, unsplit vs dimensional splitting
# Course: MA4110/MA5020 - Computational Methods for Fluid Flow (Lecture 17)
#
# The file is split into cells (lines starting with ##). In VS Code run one
# cell at a time with Shift+Enter, or run the whole file with  julia advection2d_splitting.jl
# ==============================================================================

## 0. Packages
using Plots

## 1. Grid and initial square
N = 100
dx = 1 / N
x = range(dx/2, 1 - dx/2, length=N)
square(cx, cy) = [abs(xi - cx) < 0.1 && abs(yj - cy) < 0.1 ? 1.0 : 0.0 for xi in x, yj in x]
q0 = square(0.25, 0.25)
t_final = 0.4
exact = square(0.25 + t_final, 0.25 + t_final)

## 2. Index of the upwind neighbour with periodic wrap-around
left(i) = i == 1 ? N : i - 1

## 3. Unsplit donor cell: Q_ij - νx (Q_ij - Q_i-1,j) - νy (Q_ij - Q_i,j-1); needs νx + νy ≤ 1
function unsplit_step(q, νx, νy)
    return [q[i, j] - νx * (q[i, j] - q[left(i), j]) - νy * (q[i, j] - q[i, left(j)])
            for i in 1:N, j in 1:N]
end

## 4. Dimensional splitting: an x-sweep, then a y-sweep; needs νx ≤ 1 and νy ≤ 1
function split_step(q, νx, νy)
    qs = [q[i, j] - νx * (q[i, j] - q[left(i), j]) for i in 1:N, j in 1:N]       # x-sweep
    return [qs[i, j] - νy * (qs[i, j] - qs[i, left(j)]) for i in 1:N, j in 1:N]  # y-sweep
end

## 5. Run both schemes (a = b = 1, so ν = dt/dx)
function run(step, ν)
    q = copy(q0)
    nsteps = round(Int, t_final / (ν * dx))
    for _ in 1:nsteps
        q = step(q, ν, ν)
    end
    L1 = sum(abs.(q - exact)) * dx^2
    println(rpad(string(step), 14), " ν = $ν, $nsteps steps, max = ", round(maximum(q), digits=3),
            ", L1 error = ", round(L1, digits=4), ", mass = ", round(sum(q) * dx^2, digits=6))
    return q
end
q_unsplit = run(unsplit_step, 0.4)
q_split   = run(split_step, 0.8)

## 6. Plot
p1 = heatmap(x, x, q_unsplit', aspect_ratio=1, c=:viridis, clims=(0, 1), title="unsplit, ν = 0.4")
p2 = heatmap(x, x, q_split', aspect_ratio=1, c=:viridis, clims=(0, 1), title="split, ν = 0.8")
display(plot(p1, p2, layout=(1, 2), size=(1000, 430)))

## 7. Try it: the unsplit scheme with νx + νy > 1 blows up
# q_bad = run(unsplit_step, 0.55)
