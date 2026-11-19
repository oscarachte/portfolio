# Chebyshev Spectral Collocation for a 1D Diffusion Problem

**EPFL · BSc Mechanical Engineering (BA6) · Méthodes de discrétisation en fluides · 3 ECTS**  
**Academic group project · Team of four · MATLAB**

## Overview

This project investigates a **Chebyshev spectral collocation method** for solving a stationary one-dimensional diffusion boundary-value problem. The numerical approximation is compared against a known analytical solution, with particular attention to how solution regularity affects convergence.

## Problem

We solve

\[u''(x)=f(x),\qquad x\in(0,1),\]

with Dirichlet conditions \(u(0)=0\) and \(u(1)=(1/2)^\alpha\). The reference solution is piecewise:

\[u(x)=\begin{cases}0,&x\leq 1/2,\\(x-1/2)^\alpha,&x>1/2,\end{cases}\qquad \alpha\in\{2,3,4,5\}.\]

Varying \(\alpha\) allows the study of the relationship between regularity and numerical convergence.

## Numerical approach

1. Construct Chebyshev–Gauss–Lobatto nodes on \([0,1]\).
2. Evaluate the Chebyshev basis and its second derivatives.
3. Assemble the collocation system and impose Dirichlet boundary conditions.
4. Solve the linear system in MATLAB.
5. Compare the approximation with the exact solution on a fine grid and estimate the \(L^2\) error.
6. Plot convergence curves and solution comparisons for several values of \(\alpha\).

## Findings

The report finds agreement between numerical and analytical solutions and observes convergence consistent with the supplied regularity-based estimate. For the tested cases, the spectral method becomes especially advantageous relative to a second-order finite-difference approach for smoother solutions (notably \(\alpha\geq4\)).

## Repository contents

- [`src/Discretisation_Groupe_2_Miniprojet_1.m`](src/Discretisation_Groupe_2_Miniprojet_1.m) — original project script, including the numerical study and plotting commands.

### Reproducibility note

**The original script is not standalone.** It calls `ChebyshevNodes` and `ChebyshevInterpolation`, two helper functions that were not included in the files available for this repository. To run the script as written, these functions must be supplied from the original course materials (or equivalent compatible implementations). This repository preserves the submitted script rather than inventing replacements or claiming it was executed here.

The script is designed to save `convergence_L2.png` and `solutions_comparison.png` when run with the required dependencies.

## Academic context and attribution

Developed as a **four-person group project** at EPFL. The original script credits Oscar Achte, Justin Rozzio, Jules Van Cauwenberghe, and Titouan Vaganay. Individual contributions are not assigned because the division of work is no longer documented.

The original report is not included in this public-ready package because it contains student identifiers and was prepared collaboratively. It can be added later after reviewing privacy and co-author sharing permissions.

**Topics:** numerical methods · spectral collocation · diffusion equations · convergence analysis · MATLAB
