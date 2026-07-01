---
id: HYP-3782
title: THE COVERING-MIN IS THE LORENTZIAN/MODULAR POLE, the residual is the NON-Lorentzian cusp form; + the lazy-cut prover status. PART A (LORENTZIAN merge): 'Lorentzian' (the (1,-,..,-) Minkowski/hyperbolic signature of Branden-Huh, a subclass of hyperbolic/stable polynomials) unifies the covering-min's whole regularizable structure -- (i) HYPERBOLIC: the self-concordant barrier 1/M=(n-1)+1/n (S71 HYP-3780) is a hyperbolicity-cone object and the apex is hyperbolic (2,3,7) geometry (S65); (ii) LOG-CONCAVE: the Dedekind margin |s(n,Phi6)|=T/(12T+6) is log-concave = univariate Lorentzian, rising to sup 1/12=|zeta(-1)| (S71); (iii) MODULAR: the construction margin = the classical (modular) Dedekind sum = the E2/eta -1/12 anomaly (opus-S6/S7, my S67-S70). These are one pole. The RESIDUAL is the NON-Lorentzian / NON-modular part: opus-S7 showed the covering-min BEATERS are non-modular (their margin is NOT a Zagier cotangent sum), and the covering-min irregularity a(n) = that non-modularity = the genus-1 cusp form f14 (born at the n=14 genus jump). So covering-min = Lorentzian/modular/hyperbolic bulk; residual = non-Lorentzian/non-modular cusp form (unifying S64-S71 + opus-S6/S7 under the Lorentzian banner). PART B (LAZY-CUT status): the cutting-plane covering-min prover -- search a strict beater (M<n/Phi6) via ILP feasibility + lazy danger-arc cuts (sum_{v: ||v t*||<r} x_v>=1 at each non-beater's M-witness t*); INFEASIBLE => covering-min = n/Phi6 at speeds<=n(n-1), RIGOROUS (cuts are Positivstellensatz-type, valid for any beater). n=12 = 12/133 RIGOROUS (collaborator, 208 cuts); my scipy.milp reimplementation reproduces the trajectory (215 cuts, M-witnesses trending 0.15-0.22, candidates cut off) but scipy (no warm-start, solve-time explodes past ~200 cuts) cannot close n=12/13/14 in budget -- needs a warm-starting backend (OR-tools) or multi-cut. PROOF-DIRECTION: the lazy-cut infeasibility certificate (a nonneg combination of danger constraints) is a Positivstellensatz; realizing it as a Lorentzian/hyperbolic-programming certificate would give the covering-min for all n (the hyperbolic-programming face, S71).
status: PART A synthesis (verified pieces: |s| log-concave [S71], margin=Dedekind sum [S64/S67], opus-S7 modular/non-modular; the Lorentzian unification is a framing). PART B: method correct + honest status -- n=12 rigorous stands via collaborator; my scipy prover reproduces the trajectory but is too slow to close (215 cuts still-feasible at 520s; needs warm-start/OR-tools); n=13,14 pending faster solver (task spawned). NOT a new proof of LRC14.
source: mac-mini-2026-06-30-S72
related:
  - HYP-3780   # S71 hyperconcavity: self-concordant barrier 1/M + log-concave |s| (the hyperbolic + concave)
  - HYP-3774   # S67 zeta-regularization carrier (margin = Dedekind -> -1/12)
  - HYP-3779   # S70 eta-invariant / p-adic avatars
  - HYP-3778   # klein-S60 set-cover ILP up to 4n (the lazy-cut extends this to n(n-1))
  - HYP-3777   # S69 annealing stress-test (independent confirmation covering-min=14/183)
  - HYP-3771   # S65 the (2,3,p) spine; apex-7 hyperbolic; genus/cusp-form residual
results:
  - 04-computation/lazy_cut_covering_min_prover_macmini_20260630.py
  - 05-knowledge/results/lazy_cut_covering_min_macmini_20260630.out
  - 04-computation/hyperconcavity_covering_min_macmini_20260630.py
---

# HYP-3782 -- the covering-min is the Lorentzian/modular pole; + lazy-cut status

## Part A -- "Lorentzian" unifies the covering-min's regularizable structure
**Lorentzian polynomials** (Brandén-Huh) are the polynomials whose Hessian has the **Minkowski signature
`(1, -, ..., -)`** -- a subclass of hyperbolic/stable polynomials, and exactly the multivariate home of
**log-concavity**. Reading the owner's "think Lorentzian" against the covering-min: the construction sits at a
single **Lorentzian/modular pole**, three faces of one object:
- **hyperbolic**: `1/M = (n-1)+1/n` is a self-concordant barrier (S71 HYP-3780) -- a hyperbolicity-cone object
  -- and the apex is **hyperbolic `(2,3,7)`** geometry (S65);
- **log-concave (univariate Lorentzian)**: the Dedekind margin `|s(n,Phi_6)| = T/(12T+6)` is log-concave
  (S71), rising to its supremum `1/12 = |zeta(-1)|`;
- **modular**: the construction margin **is** the classical (modular) **Dedekind sum** = the `E_2`/`eta`
  `-1/12` anomaly (opus-S6/S7; my S67-S70).

The **residual** is the **non-Lorentzian / non-modular** part. opus-S7 proved the covering-min **beaters**
(the spread family) are **non-modular** -- their margin is *not* a Zagier cotangent sum -- and the covering-min
irregularity `a(n)` *is* that non-modularity, `=` the genus-1 **cusp form `f_14`** (born at the `n=14` genus
jump). So:

> **covering-min = the Lorentzian / modular / hyperbolic bulk; residual = the non-Lorentzian / non-modular cusp
> form `f_14`.**

This is the S64-S71 + opus-S6/S7 thread under one banner: everything regularizable (hyperbolic, log-concave,
modular, `zeta`-values) is "Lorentzian"; the hard core is what is not.

## Part B -- the lazy-cut covering-min prover (method + honest status)
The cutting-plane prover of `covering-min = n/Phi_6` at the construction scale `V = n(n-1)` (where the full ILP
times out): search for a **strict beater** (a primitive covering `(n-1)`-set with all speeds `<= V` and
`M(S) < r := n/Phi_6`) by ILP feasibility [cardinality + covering + primitivity + accumulated cuts]; for each
non-beater candidate `S` (`M(S) >= r`), take its `M`-witness `t*` (where `min_v ||v t*|| = M(S) >= r`) and add
the **valid cut** `sum_{v: ||v t*|| < r} x_v >= 1` (any strict beater's open danger arcs cover `t*`).
**INFEASIBLE => no strict beater => covering-min `= n/Phi_6` for speeds `<= n(n-1)`, RIGOROUS.**

**Status (honest):**
- `n=12`: `12/133` **RIGOROUS** (collaborator's run, infeasible after 208 cuts). This closes the `(4n, n(n-1)]`
  residual for `n=12` (HYP-3778 had only `<= 4n`).
- My `scipy.milp` reimplementation **reproduces the trajectory** (215 cuts, `M`-witnesses trending `0.15-0.22`
  vs target `0.090`, candidates repeatedly cut off) but `scipy` (no warm-start, from-scratch re-solve whose
  time explodes past `~200` cuts) **cannot close** `n=12/13/14` within a `520s` budget. A warm-starting backend
  (OR-tools CP-SAT) or a multi-cut-per-iteration strategy is needed. A best-effort `n=12` run (budget `1500s`)
  is in progress; `n=13, 14` are **pending a faster solver** (task spawned).

**Proof-direction (the merge):** the lazy-cut infeasibility certificate is a **nonneg combination of danger
constraints** -- a Positivstellensatz. Realizing it as a **Lorentzian / hyperbolic-programming** certificate
(the hyperbolic-programming face of the covering-min, S71) would upgrade the finite `n(n-1)`-bounded result to
all `n` -- turning the bounded lazy-cut into a uniform proof.

## Honest scope
Part A is a synthesis (its pieces -- `|s|` log-concave, margin = Dedekind sum, opus-S7 modular/non-modular --
are verified/established; the Lorentzian unification is a framing). Part B is a correct method with an honest
status: `n=12` rigorous stands via the collaborator; my scipy prover corroborates the trajectory but is
solver-limited; `n=13, 14` are not yet closed here (need a warm-starting solver). No new LRC14 proof.
