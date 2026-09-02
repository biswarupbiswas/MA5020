# 💻 Course Codes & Numerical Implementations

This page provides standalone code implementations for the computational methods and solution techniques discussed in the lectures. Reference codes are provided preferably in **Julia** (students are free to implement their solutions in Python, MATLAB, C++, or their preferred language). Students can download, run, and modify these scripts locally.

> ⚡ **New to Julia?** Check out our quick 10-minute [Julia Guide for Computational Methods](julia.md).

---

## 🚀 Code Repository

| Topic / Method | Language | Associated Lecture | Source File |
| :--- | :--- | :--- | :--- |
| **Linear Advection Equation: Method of Characteristics** | Julia (`.jl`) | [Lecture 4](pdfs/lecture04.pdf) | [`linear_advection_characteristics.jl`](code/linear_advection_characteristics.jl) |
| **Inviscid Burgers' Equation: Method of Characteristics** | Julia (`.jl`) | [Lecture 5](pdfs/lecture05.pdf) | [`burgers_characteristics.jl`](code/burgers_characteristics.jl) |
| **Linear Advection Equation: First-Order Upwind Scheme** | Julia (`.jl`) | [Lecture 6](pdfs/lecture06.pdf) / [Lecture 7](pdfs/lecture07.pdf) | [`advection_upwind.jl`](code/advection_upwind.jl) |
| **Linear Advection Equation: Lax-Wendroff Scheme (2nd-Order Central)** | Julia (`.jl`) | [Lecture 7](pdfs/lecture07.pdf) / [Lecture 8](pdfs/lecture08.pdf) | [`advection_lax_wendroff.jl`](code/advection_lax_wendroff.jl) |
| **Linear Advection Equation: Beam-Warming Scheme (2nd-Order Upwind)** | Julia (`.jl`) | [Lecture 8](pdfs/lecture08.pdf) | [`advection_beam_warming.jl`](code/advection_beam_warming.jl) |

---

## 🛠️ Instructions for Running Codes

### Julia Scripts
To execute any of the Julia programs:
1. Make sure [Julia](https://julialang.org/downloads/) is installed.
2. Install the `Plots` package (run once in Julia):
   ```julia
   using Pkg; Pkg.add("Plots")
   ```
3. Run the script directly from your terminal to view the live animated plot:
   ```bash
   julia linear_advection_characteristics.jl
   # or
   julia burgers_characteristics.jl
   ```
   A dynamic plot window will open and animate the wave propagation / steepening in real time.
