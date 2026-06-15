# Baby Hodge is the det/permanent wall read through the moment problem

**Source:** mac-mini-2026-06-15-S1. Dispatch: "tr(A^k) = Σ over cycle-configs — power
sums are the moments of the tournament spectrum, decomposing into cycle-counts +
overlap corrections; work on the baby Hodge problems." Canon: THM-509, HYP-2526..2528,
OPEN-Q-099, T820. Builds on THM-118/499/500/502/505/506/507.

## The picture in one sentence

The realizable-invariant-vector problem ("which (c₃, c₅, …, H) arise from a tournament")
is a **truncated moment problem**, and its **holes** — the non-algebraic-Hodge-class
analogs — are exactly the **integer lattice points interior to the continuous moment
body that the #P count skips**. The continuous body is `det`-side (spectral) and blind
to the holes; the skip is `permanent`-side (integrality + conflict). Baby Hodge is the
**Valiant det/permanent = P/#P wall** read through the moment problem.

## The three layers, and which one cuts

| layer | object | cuts a hole? |
|---|---|---|
| **Spectral / det-side** | skew Hankel `[m_{i+j}]`, `m_r = tr((SS^T)^r)` | **NO** — PSD for every tournament (`SS^T = −S²⪰0`, Stieltjes), and *constant on cospectral classes* (THM-507), so it takes the same value on a hole and a cospectral realized point |
| **Continuous overlap / flag** | Cauchy–Schwarz Gram of overlap densities (p33, α₂) | **NO** — the hole is a convex (tournament-limit) interior point, e.g. `(8,10) = ½(8,8)+½(8,12)`; any SOS/flag matrix is PSD there |
| **Integrality / #P** | the count is an integer; `H = I(Ω,2)` packing | **YES** — the only thing that skips an interior lattice point |

So the **det-side Hodge–Riemann inequalities are real but non-cutting** (proved), the
**continuous flag/positivity relaxation is also non-cutting** (verified: the holes are
convex-interior tournament-limit points), and the entire obstruction is **integrality on
the permanent side**. This is the sharp content of THM-509.

## Why this is genuinely Hodge-shaped

A non-algebraic Hodge class is a rational `(p,p)` cohomology class that satisfies every
Hodge-theoretic (positivity, type) condition yet is not the class of an algebraic cycle —
the failure is one of *integral/algebraic realizability*, not of any continuous
inequality. The baby-Hodge holes have exactly this shape: each one (verified, all 36 at
n = 6, 7) satisfies every continuous moment inequality — the spectral Hankel ones AND the
convex/flag ones — and fails only the integral realizability "is it `tr A⁵` of an actual
tournament." The flagship `(c₃,c₅) = (8,10)` at n = 6: `tr A⁵ = 50 = ⅓·40 + ⅔·55`, a
convex combination of two realized spectra, score-stratification-forbidden as an integer
point. A non-algebraic Hodge class with a one-line certificate.

## The moment–cumulant engine (and the free-probability mirage)

The `tr(A^k) = Σ_configs rot(C)·emb(C)` census (THM-502) is an exact **moment–cumulant
pair**: moments `p_k = tr(A^k)` (closed-walk power sums), cumulants `W_k =
(1/k)Σ_{d|k} μ(d) p_{k/d}` (Witt/necklace numbers), joined by the **Artin–Mazur zeta**
`exp(Σ p_k u^k/k) = Π (1−u^k)^{−W_k}` — Möbius inversion on the **divisor lattice**. The
"overlap corrections" the user named are precisely `W_k − c_k` (the reducible/aperiodic
configs): the cumulant `W_k` is spectral; its **split** into a simple cycle `c_k` and
overlaps is the non-spectral content — the mechanism of `c_k`-non-spectrality from k = 6.

It is tempting to call this "free probability." It is **not** (verified, HYP-2526): the
free cumulants (non-crossing-partition Möbius) differ numerically from `W_k`, and the
posets are different sizes (divisor lattice 1,2,2,3,2,4 vs Catalan 1,2,5,14,42,132). The
project's *genuine* free-probability content (THM-438: Catalan/semicircle in the Paley
spectrum) lives in the **spectral distribution's** moments — a different object from this
**cyclic necklace census**. Two distinct "cumulant" stories; do not conflate them. The
honest name here is the **necklace/zeta moment–cumulant pair**.

## How this unifies the week's threads

- It **explains `d ⊥ H`** (the determinant-lens orthogonality, THM-468/499): `d =
  det(I+S) = Π(1+μ_j²)` is a det-side spectral coordinate (cospectral-constant); `H` is
  the permanent-side packing count. They are orthogonal because they sit on **opposite
  sides of the Valiant wall** — the same wall that creates the baby-Hodge holes.
- It connects to the **BSD/Hodge "self-dual middle"** reflection (the alternating →
  Pfaffian-square spine): there the rigorous content was that a det-side (alternating)
  pairing forces a square; here the det-side (spectral) data is *everything the
  determinant can see* and is provably blind to the realizability holes. Both say: the
  determinant/spectral world is clean and closed; the hard, hole-creating content is on
  the permanent/#P side.

## Open (OPEN-Q-099)

Prove the all-n **positivstellensatz**: no polynomial moment inequality (any degree,
spectral or overlap-side) cuts a baby-Hodge hole — making "hole = integer lattice point
interior to the flag-feasible body, skipped by the #P count" a theorem. Evidence (n=6,7,
explicit Gram) points firmly that way: the holes are pure integrality gaps, the cleanest
possible non-algebraic-Hodge-class analog.
