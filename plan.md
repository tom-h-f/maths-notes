
## Overview
- Structure: Phase-based (no calendar dates). Move forward when you meet mastery criteria.
- Learning loop each session:
  1) Concept: read/watch, derive key ideas.
  2) Practice: targeted problems with solutions.
  3) Code: Julia exploration to test/visualize.
  4) Review: error log + spaced repetition.
  5) penis4

---

## Tooling Setup
- Julia packages (Pkg mode: press `]`):
  - Core: `IJulia`, `Pluto`, `Plots`, `StatsPlots`, `LinearAlgebra`, `SparseArrays`, `Random`
  - Math/Analysis: `SpecialFunctions`, `QuadGK`, `ForwardDiff`, `NLsolve`, `Optim`
  - Differential equations: `DifferentialEquations`, `ModelingToolkit`, `StaticArrays`
  - Probability/Stats: `Distributions`, `StatsBase`, `DataFrames`, `CSV`, `GLM`
  - Signals/FFT (optional): `FFTW`
- Suggested install:
  - `] add IJulia Pluto Plots StatsPlots SpecialFunctions QuadGK ForwardDiff NLsolve Optim DifferentialEquations ModelingToolkit StaticArrays Distributions StatsBase DataFrames CSV GLM FFTW`
- Notes workflow:
  - Use Pluto.jl or Jupyter for each phase. Keep an Obsidian vault with:
    - Cheatsheets (definitions/theorems)
    - Worked examples
    - Problem solutions
    - Experiment notes and plots
- Coding conventions:
  - Keep functions pure; document with docstrings; plot with consistent styles.
  - Save small utilities per phase in a `utils.jl`.

---

## Mastery Criteria
- 60–100 problems per phase (varied difficulty, including applications).
- 1–2 Julia mini-projects per phase with a short write-up (What/Why/How/Results).
- Timed checkpoint (self-quiz or subset of an OCW exam).
- One-page summary sheet: key definitions, theorems, techniques, pitfalls.
- Explain the “big ideas” aloud and solve mixed problems without aids.

---

## Resource Key
- Primary (as requested): Paul’s Online Math Notes (algebra → calculus → DE).
- Free textbooks: OpenStax (Precalculus, Calculus, Linear Algebra, Statistics).
- MIT OCW: 18.01, 18.02, 18.03, 18.06, 6.042J (lectures, problem sets, exams).
- Intuition boosts: 3Blue1Brown (calculus, linear algebra), StatQuest (stats).
- Alternatives: Stewart (Calculus), Strang (Linear Algebra), Boyce & DiPrima (ODE),
  Blitzstein & Hwang (Probability), Ross (Probability), Rosen (Discrete).

---

## Phase 0: Algebra/Precalculus Refresh

**NOTE: We are using NotebookLM as our source of knowledge.**
### Targets
- Functions/graphs; exponentials/logs; trig identities; inequalities; complex numbers; basic vectors.

### Primary
- Paul’s Notes: Algebra + Precalculus

### Alternatives
- OpenStax Precalculus; Khan Academy (Algebra II/Precalculus)

### Practice (select 50–80 problems)
- Mixed Paul’s problem sets on functions, logs/exp, trig identities, inequalities, complex numbers.

### Julia Exercises
- Function transforms: compose/shift/scale and visualize changes.
- Root-finding (bisection) for polynomial and exp/log equations; compare brackets.
- Complex plane:
  - Plot roots of `z^n − 1`; verify De Moivre numerically for random `θ`.

### Mini-Project
- Function Explorer (Pluto):
  - Sliders for `a, b, c` in `f(x) = a*exp(bx) + c*sin(x)`.
  - Linearize exponential data with `log`; semilog/lin-log comparisons.

---

## Phase 1: Single-Variable Calculus
### Targets
- Limits/continuity; derivatives/applications; definite integrals; FTC; methods of integration; sequences/series; Taylor/Maclaurin; convergence tests.

### Primary
- Paul’s Calculus I–II notes + problems

### Alternatives
- OpenStax Calculus Vol. 1–2; MIT OCW 18.01; Stewart Calculus

### Practice (80–120 problems)
- Differentiation, optimization, related rates, Riemann sums, substitution, parts,
  partial fractions, improper integrals, series tests, Taylor approximations.
- Checkpoint: one MIT 18.01 exam.

### Julia Exercises
- Numerical derivative (central difference) with step-size sweep; plot truncation vs. round-off error.
- Newton’s method with safeguards; visualize iterations on 1D and basins for 2D extensions.
- Quadrature: trapezoid/Simpson vs. `QuadGK`; error vs. N plots.
- Series approximations: truncated `exp`, `sin`, `cos`; compare to `SpecialFunctions`.

### Mini-Projects
- Optimization Sandbox:
  - Compare Newton, secant (`NLsolve`, `Optim`) on nonconvex functions; analyze convergence.
- Quadrature Explorer:
  - Oscillatory integrals; adaptive vs. fixed-step behavior; performance notes.

---

## Phase 2: Introductory Proofs (Parallel, Light)
### Targets
- Logic, sets, functions, induction; direct/contradiction; epsilon–delta intuition.

### Primary
- Hammack: Book of Proof (free)

### Alternatives
- Velleman: How to Prove It

### Practice (20–40 short proofs)
- Basic set/logic proofs, induction (sum, product, divisibility), simple analysis proofs.

### Julia Exercises
- Empirical identity tester: randomized checks of algebraic/trig identities to form intuition.
- Counterexample search: heuristics to find edge cases violating a conjecture.

### Mini-Project
- Continuity vs. differentiability gallery:
  - Construct piecewise/absolute value/fractal-like functions; visualize and annotate properties.

---

## Phase 3: Linear Algebra
### Targets
- Systems and elimination; LU/QR; vector spaces, linear maps; bases/dimension; rank-nullity; eigenvalues/eigenvectors; orthogonality; least squares; SVD (basics).

### Primary
- MIT OCW 18.06 (Strang) + problem sets or Strang’s text

### Alternatives
- Axler: Linear Algebra Done Right (conceptual), OpenStax Linear Algebra, 3Blue1Brown videos

### Practice
- 18.06 sets/exams; rank-nullity proofs; subspace/basis; conditioning and projections.

### Julia Exercises
- Implement Gaussian elimination with partial pivoting; compare to `LinearAlgebra.lu`.
- QR via Modified Gram–Schmidt and Householders; test on ill-conditioned matrices.
- Power method (dominant eigenvalue); simple deflation for next eigenpairs.
- Least squares via normal equations and QR; compare residuals/conditioning.

### Mini-Projects
- Data Fitting:
  - Polynomial vs. spline fits on noisy data; cross-validation; bias–variance discussion.
- Image Compression with SVD:
  - Reconstruct images with k singular values; PSNR vs. k; storage trade-offs.

---

## Phase 4: Multivariable Calculus
### Targets
- Partial derivatives, gradients/Jacobians; chain rule; multiple integrals; change of variables; vector fields; line/surface integrals; Green/Stokes/Divergence.

### Primary
- Paul’s Calculus III notes

### Alternatives
- OpenStax Calculus Vol. 3; MIT OCW 18.02

### Practice (60–100 problems)
- Tangent planes, gradient-based optimization, Lagrange multipliers, double/triple integrals, vector calculus identities, flux/circulation.

### Julia Exercises
- Visualize scalar fields and gradient directions; contour + quiver plots.
- Automatic differentiation with `ForwardDiff` for gradients/Jacobians; compare to finite differences.
- Monte Carlo integration in 2–5D; error vs. N; importance sampling (light).
- Line integrals over parameterized curves; verify path independence by checking curl.

### Mini-Projects
- Change of Variables:
  - Integrate over an irregular region via a mapping; compute Jacobian; validate numerically.
- Vector Field Clinic:
  - Sample/plot random fields; estimate curl/divergence numerically; streamline visualizations.

---

## Phase 5: Ordinary Differential Equations
### Targets
- First-order methods, linear ODEs with constant coefficients, Laplace transforms, systems of ODEs, stability and phase planes, convolution, numerical methods (Euler/RK4).

### Primary
- Paul’s Differential Equations notes

### Alternatives
- MIT OCW 18.03; Boyce & DiPrima

### Practice (60–100 problems)
- Analytical solution methods; Laplace applications; system stability classification.

### Julia Exercises
- Implement Euler, Heun, RK4 from scratch; stability on `y' = λy`.
- `DifferentialEquations.jl`: nonstiff vs. stiff solvers (`Tsit5`, `Rodas5`); event handling.
- Phase portraits for 2D linear systems; classify using eigenvalues.

### Mini-Projects
- Mass–Spring–Damper:
  - Analytic vs. numeric solutions; step and frequency responses; parameter sweeps.
- RC/RLC Circuits:
  - Laplace-transform sketch vs. numeric simulation; Bode-style plots using `FFTW` if needed.

---

## Phase 6: Probability and Statistics
### Targets
- Counting/combinatorics; random variables/distributions; expectation/variance; LLN/CLT; joint/conditional distributions; Bayes; MLE; confidence intervals; hypothesis tests; linear regression.

### Primary
- Blitzstein & Hwang: Introduction to Probability

### Alternatives
- Ross: A First Course in Probability; OpenIntro Statistics; Montgomery (engineering stats)

### Practice (80–120 problems)
- Mix puzzles, distribution calculations, estimators/tests, and regression analysis.

### Julia Exercises
- `Distributions.jl`: simulate PMF/PDF/CDF/quantiles; inverse transform sampling.
- Monte Carlo estimators: estimate π; integrate functions; variance reduction (antithetic/control variates).
- Bootstrap: empirical CIs; compare to normal-theory CIs.
- `GLM.jl`: fit linear regression; residual diagnostics; CI/PI.

### Mini-Projects
- CLT Lab:
  - Sums of varied distributions; visualize convergence; Berry–Esseen style numerics.
- MLE Sandbox:
  - Implement MLE for Poisson/Exponential/Normal; profile likelihood; simulation study of estimator bias/variance.

---

## Phase 7: Discrete Mathematics (CS-Focused)
### Targets
- Logic/sets; induction/strong induction; recurrences; graphs/trees; counting; light generating functions.

### Primary
- MIT OCW 6.042J Mathematics for Computer Science

### Alternatives
- Rosen: Discrete Mathematics and Its Applications (selected chapters)

### Practice
- Proof-based and computational problems: graph properties, recurrences, counting arguments.

### Julia Exercises
- Graphs using `Graphs.jl`:
  - Implement BFS/DFS/topological sort; shortest paths; connectivity checks.
- Recurrences:
  - Solve numerically and verify closed forms; use generating functions for simple cases.
- Combinatorics:
  - Randomized checks of identities; Monte Carlo estimates for large counts.

### Mini-Projects
- Scheduler:
  - Task dependency DAG; topo orders; cycle detection; critical path.
- Random Graphs:
  - Explore connectivity thresholds in `G(n, p)`; plot probability of connectedness vs. `p`.

---

## Phase 8: Numerical Methods
### Targets
- Floating-point/round-off/conditioning; root-finding; linear solves & stability; interpolation; quadrature; ODE solvers; intro to PDE discretization (finite differences).

### Primary
- Burden & Faires (selected) or Heath: Scientific Computing

### Alternatives
- MIT 18.330; UW AMATH 301 materials

### Practice
- Error analysis; stability vs. conditioning; convergence studies on test problems.

### Julia Exercises
- Floating-point:
  - Estimate machine epsilon; demonstrate catastrophic cancellation; Kahan summation.
- Linear solves:
  - Compare `\` (LU), QR, and CG on SPD systems; condition number effects.
- Interpolation:
  - Polynomial interpolation (Runge phenomenon) vs. cubic splines; error plots.
- PDE intro:
  - 1D heat equation with explicit vs. implicit schemes; CFL stability exploration.

### Mini-Projects
- Ill-Conditioned Systems:
  - Hilbert matrices; error amplification; Tikhonov regularization; L-curve style plot.
- Adaptive Quadrature:
  - Implement recursive Simpson; compare to `QuadGK` on challenging integrands.

---

## Optional Depth Tracks
- PDEs
  - Haberman: Applied PDEs; MIT 18.303/18.085/18.086 resources
  - Julia: Method of lines with `DifferentialEquations.jl`; 2D Poisson via finite differences (`SparseArrays`) + iterative solvers.
- Optimization/Convexity
  - Boyd & Vandenberghe (free PDF)
  - Julia: `Optim.jl`, `Convex.jl` + SCS; gradient checks via `ForwardDiff`.
- Signals/Systems
  - Oppenheim/Willsky; MIT 6.003
  - Julia: `FFTW.jl`; filter design; convolution in time vs. multiplication in frequency.
- Control
  - Ogata or Nise; MIT 6.302/6.241
  - Julia: `ControlSystems.jl`; Bode, Nyquist, root locus; state-space design/observers.
- Complex Analysis
  - Brown & Churchill
  - Julia: Visualize complex mappings; numeric residues/integrals.

---

## Assessment and Progression
- For each phase:
  - Pass a timed self-quiz or OCW subset.
  - Complete 1–2 mini-projects with clear write-ups.
  - Produce a one-page cheat sheet and error-log summary.
- Move forward when:
  - You can articulate core ideas and solve mixed problems without aids.
  - Julia experiments align with analytic expectations within known error bounds.

---

## Flexible Weekly Structure
- Session A (concepts, 60–120 min): Read/derive; 1–2 proofs/derivations.
- Session B (practice, 60–120 min): 8–15 problems; check with solutions.
- Session C (coding, 30–90 min): Implement/visualize; record findings.
- Optional Session D (90–150 min): Exam-style set or project progress.
- Review (10–15 min each session): spaced repetition + error log.

---

## How to Choose Resources
- Prefer Paul’s Notes for concise, example-driven learning.
- Pair with MIT OCW for structured problem sets and exams.
- Add 3Blue1Brown/StatQuest when intuition lags.
- Use OpenStax for free, textbook-style continuity.

---

## Getting Started: Next Steps
1) Install Julia packages (see Tooling Setup).
2) Create a Phase 0 Pluto notebook with a cheatsheet (logs, trig, inequalities).
3) Solve 15–20 Paul’s Algebra/Precalc problems; implement bisection and complex-root plots.
4) Start Paul’s Calc I (limits/derivatives); code Newton’s method + derivative step-size study.
5) Begin your error log and flashcards for key formulas and theorems.

---

## Notes
- Phases are independent of calendar time—advance on mastery.
- Interleave light proofs (Phase 2) alongside Phases 1–4.
- Tie math to Julia experiments frequently for deeper intuition and retention.
