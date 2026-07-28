---
id: THM-2835
title: "Degree-four Boolean functions have maximal point sensitivity nine; record exponent needs degree >= 5"
status: >
  FINITE-EXACT, SINGLE ENGINE WITH BIDIRECTIONAL CONTROLS — AWAITING
  INDEPENDENT ENGINE AUDIT.  CP-SAT decides (d,s)=(4,10) INFEASIBLE
  (n=10, 375 free bits, 97 s wall), so no Boolean function of real degree
  <= 4 has 10 sensitive coordinates at any input; NAE∘NAE attains 9, so
  m(4) = 9 sharp.  The same encoder returns SAT on (4,9) (control) and its
  INFEASIBLE verdict on (3,7) is independently confirmed by the exhaustive
  C engine of THM-2832 (control of the UNSAT direction).  PROVED lemma:
  XOR on disjoint variables gives m(d1+d2) >= m(d1)+m(d2), complementing
  composition multiplicativity m(d1 d2) >= m(d1) m(d2).  COROLLARY (on
  THM-2832 + this): any multilinear Boolean P with 0-sensitivity
  s0 > deg^{log_3 6} has degree >= 5; at degree 5 the only open cells are
  s in {14,...,25} (5^{log_3 6} = 13.78, Huang cap 25).
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Degree vs Sensitivity for Boolean Functions")
depends_on:
  - THM-2832-degree-three-boolean-max-point-sensitivity-six
related: []
script: 04-computation/sens_degree_fullsens_cpsat_macmini_S171.py
output: 05-knowledge/results/sens_degree_fullsens_macmini_S171.out
script_sha256: see THM-2832 companion files
output_sha256: see THM-2832 companion files
hash_basis: LF-normalized bytes
---

# THM-2835 — m(4) = 9; beating exponent log_3 6 requires degree >= 5

## Statement

Let m(d) be the maximum number of coordinates at which a Boolean function of
real multilinear degree <= d can be sensitive at a single input (WLOG the
input 0 with f(0)=0, all listed f(e_i)=1; THM-2832 reduction).  Then

    m(1) = 1,   m(2) = 3,   m(3) = 6   (THM-2832, two engines),
    m(4) = 9   (this theorem, CP-SAT INFEASIBLE at (4,10); NAE∘NAE = 9).

## Engine and controls

OR-Tools CP-SAT model: Boolean variable per point of {0,1}^10, one integer
equality per subset S with |S| >= 5 (vanishing Moebius coefficient
sum_{T subset S} (-1)^{|S|-|T|} f(T) = 0), degree-sequence symmetry breaking
on the 45 pair values.  Verdict INFEASIBLE in 97 s (8 workers).
Controls through the *same code path*:
  * (4,9)  -> FEASIBLE, returning an explicit degree-4 fully-sensitive
    truth table on 9 variables (SAT direction validated);
  * (3,7)  -> INFEASIBLE, matching THM-2832's independent exhaustive C DFS
    (UNSAT direction validated against a second engine).
Raw exhaustive C DFS on (4,10) was abandoned at 4.3e9 nodes (~1% of tree);
(4,7)/(4,8) solution libraries are too large for a layered second engine.
An independent engine (e.g. PB-SAT with DRAT proof, or orbit-canonical DFS)
is REQUESTED before promotion past audit; until then treat m(4)=9 as
verified-single-engine.

## Superadditivity lemma (PROVED, elementary)

If f, g are fully sensitive at 0 on disjoint variable sets (sizes s_f, s_g,
degrees d_f, d_g, f(0)=g(0)=0), then h = f XOR g = f + g - 2fg is Boolean,
multilinear, h(0)=0, fully sensitive at 0 on all s_f + s_g variables, and
deg h = d_f + d_g (the top product terms cannot cancel: the coefficient of
a top monomial of fg is the product of nonzero top coefficients... verified
on generators used below).  Hence  m(d1 + d2) >= m(d1) + m(d2).
Together with composition (m(d1 d2) >= m(d1) m(d2)) the best known lower
bounds are:  m(5) >= 10 (9+1), m(6) >= 18 (6*3), m(7) >= 19, m(8) >= 27
(3*9), m(9) >= 36 (6*6), m(12) >= 54.  All satisfy m(d) <= d^{log_3 6},
consistent with the conjecture that Kushilevitz tensorization is extremal.

## Consequence for the Epoch problem

Any multilinear integer-coefficient P, Boolean on the cube, with
sum_i |P(0)-P(e_i)| = deg(P)^a and a > log 6/log 3 must have deg(P) >= 5
and at least 14 sensitive coordinates.  The first undecided cells are
(5, 14..25); (5,14) is CP-SAT-expressible (16k variables) but exceeded this
host's 8 GB during search — flagged as compute-wanted for a larger machine.
