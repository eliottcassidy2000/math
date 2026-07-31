---
id: THM-2832
title: "Degree-three Boolean functions have maximal point sensitivity exactly six"
status: >
  FINITE-EXACT (two independent engines, exhaustive) + PROVED reduction.
  A Boolean function of real multilinear degree at most 3, on any number of
  variables, has at most 6 sensitive coordinates at every input; Kushilevitz's
  6-variable cubic attains 6, so m(3) = 6 is sharp.  Also m(2) = 3 (NAE sharp),
  re-proved by hand and by both engines.  This closes the degree-3 route to
  improving the sensitivity-vs-degree separation exponent log 6 / log 3.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Degree vs Sensitivity for Boolean Functions")
depends_on: []
related:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
script: 04-computation/sens_degree_fullsens_search_macmini_S171.c
output: 05-knowledge/results/sens_degree_fullsens_macmini_S171.out
script_sha256: 1b9757c276e35329fe7a43e2d5a2493c2d52005939435b91f55fa88d78e805b9
output_sha256: b6c5dfb74a7b02c3739a736ac6d8beef5995258086490653318df3ee9a83ecad
hash_basis: LF-normalized bytes
---

# THM-2832 — max 0-sensitivity of cubic Boolean functions is 6

## Statement

Let `f : {0,1}^n -> {0,1}` be any function whose unique multilinear real
polynomial representation has degree `<= 3`.  Then at every input `x`, the
number of coordinates `i` with `f(x ⊕ e_i) != f(x)` is at most **6**.
The bound is attained (Kushilevitz's function, n = 6, degree 3, fully
sensitive at 0).  For degree `<= 2` the sharp bound is **3** (NAE).

## Reduction (proved)

Sensitivity at `x` is invariant under the substitution `x_i -> 1 - x_i`
(which preserves multilinearity and degree), so WLOG `x = 0`; replacing `f`
by `1 - f` fixes `f(0) = 0`.  Setting all non-sensitive coordinates to `0`
preserves the sensitive values and cannot raise the degree.  Hence a
degree-`d` function with `s` sensitive coordinates at a point restricts to a
**fully sensitive** function on exactly `n = s` variables with `f(0) = 0`,
`f(e_i) = 1` for all `i`, degree `<= d`.  The existence question for the pair
`(d, s)` is therefore a finite decision problem on `2^s` values.

## Computation (exhaustive, two engines)

Encoding: the Moebius coefficient of every subset `S` with `|S| > d` vanishes;
processing subsets in numeric mask order makes each such value a forced
integer that must lie in `{0,1}` — an exhaustive DFS with forced-level
propagation (engine 1, C).  Engine 2 is an independent OR-Tools CP-SAT model
of the same constraint system with degree-sequence symmetry breaking.

| (d, s) | verdict | witness / count |
|--------|---------|-----------------|
| (2, 3) | SAT     | NAE (unique) — control |
| (2, 4) | UNSAT   | 33 nodes — matches the hand cascade proof |
| (3, 6) | SAT     | 12 solutions, all genuinely degree 3 — control |
| (3, 7) | **UNSAT** | engine 1: 818,943 nodes exhaustive; engine 2: INFEASIBLE |

Hand proof for (2,4), recorded as a control: pair coefficients satisfy
`c_ij in {-2,-1}`; each triple value `3 + c_ij + c_ik + c_jk in {0,1}` forces
all `c = -1` (a `-2` would drive the sum below `-3`); then
`f(1111) = 4 - 6 = -2`, contradiction.

## Consequences and boundary

* Kushilevitz's construction is **optimal at degree 3** — not merely the best
  known.  Any improvement of the separation `s0 >= deg^{log_3 6}` must use
  degree `>= 4`.
* Next open cells (monotone in `s` by the restriction argument):
  `(4,10)` — SAT would give exponent `log_4 10 ~ 1.661 > log_3 6 ~ 1.6309`;
  UNSAT gives `m(4) = 9` (NAE∘NAE attains 9).  `(5,14)`, `(6,19)` follow.
  Composition gives `m(d1 d2) >= m(d1) m(d2)`, so `m(9) >= 36`.
* Scope: this is 0-sensitivity = full sensitivity at a point; it implies the
  same bound for the usual (max) sensitivity of degree-3 functions at any
  single point, i.e. `s(f) <= 6` for every cubic Boolean `f`.
