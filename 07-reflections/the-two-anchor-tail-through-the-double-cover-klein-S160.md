---
source: klein-2026-07-07-S160 (HYP-4841)
status: EXACT identities + independent confirmation + a dual-structure discovery.
  The 2-anchor tail (boxeph's discharge-all-legs object) factors through the Z/2 double
  cover: gap@0(2E,x) = 2(min r + min l) pointwise-exact; all-odd families have identical
  anchor laws (PA_2 = meas(T ∪ T−½)); mixed parity drops the ½-anchor tail by ~0.3.
  boxeph's inf table independently confirmed (discharges all legs). The pairwise-LP duals
  are INTEGER graphic forms: exactly a Hunter spanning tree at the barrier shape, tree +
  integer correction layer at the AP — the closed-form weighted-cherry candidate.
tags:
  - lonely-runner
  - LRC14
  - two-anchor
  - double-cover
  - voltage-lift
  - weighted-cherry
  - duality
  - relativity
---

# The two-anchor tail through the double cover

**klein-2026-07-07-S160.** Owner: work the 2-anchor joint tail and a closed-form weighted
cherry; tournaments from creative pair-statistic cutoffs; the relativity of perspectives;
BB(27)/Goldbach compression as inspiration.

## 1. The double-cover dictionary for PA₂ (exact)

boxeph's reduction made `PA₂(E) = P(max(gap@0, gap@½) > 1/7)` the discharge-all-legs
object. Writing `r_a/l_a` for the one-sided min distances at anchor `a ∈ {0, ½}`:

- **(i) `gap@0(2E, x) = 2·(min(r₀,r_½) + min(l₀,l_½))`** — machine-verified to 1e-16
  pointwise; one-line proof from `frac(2t) = 2·dist(t, {0,½})`-bookkeeping. The doubled
  family reads the TWO anchors of E collapsed into one; hence
  `gap@0(2E) ≤ 2·min(gap@0, gap@½) ≤ 2·max(...)` pointwise. The ½-anchor is not an extra
  device — it is the second sheet of the 2-lift (the same voltage structure as the record
  families, S157).
- **(ii) the parity dichotomy.** If E is ALL-ODD, the orbit at `x+½` is the orbit at `x`
  shifted by exactly ½, so `gap@½(E,x) = gap@0(E,x+½)` and the two anchor laws COINCIDE
  (verified: diffs 0.0000). Then **`PA₂ = meas(T ∪ (T−½))`** with `T = {gap@0 > 1/7}` —
  the 2-anchor tail is literally the 1-anchor tail set unioned with its half-rotation:
  a clean object for exact work (roof-computable on AP-like shapes). For MIXED-parity
  families `P(gap@½ > 1/7)` drops by 0.29–0.33 (the even part crowds the ½-anchor: its
  points are the 0-anchor points of the half-speed family) — the max still gains, but the
  mechanism differs by 2-adic type. Any proof of `inf PA₂ ≥ T_k` should case-split on
  parity type; the all-odd case reduces to a UNION-of-translates estimate on the single
  tail set.
- **(iii) independent adversarial confirmation of boxeph's table:** my jump-move infs
  `0.798 (k=8), 0.714 (k=9), 0.429 (k=13)` vs their `0.766/0.685/0.360` — all comfortably
  above `T_k = 0.6185/0.5057/0.0565`. The discharge stands from two independent engines.

## 2. The weighted-cherry duals are integer graphic forms (the discovery)

Extracting the optimal duals of the S159 pairwise LP (certificate form
`W ≥ y₀ + Σᵢyᵢθ + Σᵢⱼyᵢⱼmᵢⱼ`, feasibility = the quadratic form is dominated by `1[S=∅]`):

- **At the barrier shape `{1,3,5,6,7,9,25,38}`** (LP = 0.1233): the dual's support is
  singles (unit weights) plus SIX unit-weight pairs — `(13,32),(13,33),(13,35),(29,37),
  (31,32),(31,37)` — which form **exactly a spanning tree** on the seven differences
  (star-at-13 + path). **Hunter is dual-OPTIMAL at spread shapes**, not merely close: the
  S155 tree bound was the sharp pairwise certificate there, and the LP names the tree.
- **At the AP** (LP = W = 0.3347): unit singles plus an INTEGER pair layer — `+2` on
  `(3,7)`, `+1` on `(1,6)`, `−1` on `{(1,7),(6,7),(1,2),(3,6),(2,4),(3,5)}` (solver sign
  convention to be audited before asserting the certificate direction). The optimal dual
  is a tree PLUS an integer signed correction layer.

**Closed-form candidate (the weighted-cherry theorem to prove):** the optimal pairwise
dual is always an INTEGER graphic form — a spanning tree plus signed integer corrections
on few extra pairs — with the AP's layer determined by the law's G-table (its `+2` sits on
`(3,7)`: `r=(3,0)`, a resonant pair; its `+1` on `(1,6)`: the max-G pair `G=6/49`...).
Next: pin the sign convention, re-verify feasibility of the extracted forms by direct
`f(S) ≤ 1[S=∅]` enumeration (128 checks), and test transfer: a feasible integer form is a
VALID certificate at every shape with law-computable value.

## 3. The relativity frame (the owner's hint, recorded)

The observer formulation (13 moving + 1 fixed) is one ROW of the antisymmetric
relative-speed matrix `d_ij = e_j − e_i` on `E ∪ {0}` — a WEIGHTED TOURNAMENT, the
intersubjective binary relation itself. Every runner's perspective carries its own
leave-one-out difference family `D_i` and loneliness events; my `W_i`/`B` statistics
(S155) are the per-perspective row objects, `x ↦ 1−x` is global arc reversal (`T^op`),
and kps-S64's THM-640 Paley bridge is the QR-cutoff quotient of this same matrix. New
cutoffs proposed (leads, deliberately imaginative per the owner): the **valuation
tournament** (orient by `v₇(d_a)` comparisons — the 7-adic tree of S157), the **mass
tournament** (orient perspectives by Hunter-tree weight: who is most defendable), the
**anchor tournament** (orient anchors by stochastic dominance of their gap laws — the
parity dichotomy above says {0} beats {½} exactly on mixed-parity families). Untested
beyond definitions; the test is whether any invariant predicts LP/PA₂ minimizers.

## 4. The compression note (BB(27)/Goldbach)

Yedidia–Aaronson (arXiv:1505.02479) compress Goldbach/RH/ZFC-independence into tiny
explicit machines; the community's 27-state Goldbach machine makes BB(27) ≥ deciding
Goldbach. The fleet's week has the same shape: (A′) compressed into an 18-anchor
statistic (opus F₇), then a 2-anchor tail (boxeph), pair structure into a 49-class law
(THM-638), legs into finite rational ledgers (kps-S60) — each step trades an infinite
quantifier for a small exact object plus one honest residual. The Goldbach echo in the
mathematics itself: everything hard here lives on additive pair structure (midpoint
relations drive bisection alignment; the law is graded by residue pairs) — the density
floor is, like Goldbach, a statement that additive structure cannot conspire everywhere
at once. Guardrail from the same source: MISTAKE-116 (unbounded covering) is this
project's own lesson that "compressed into finite" must be proved, not pattern-matched.

## Honest status

(1i)/(1ii) exact + verified; (iii) numeric adversarial; the dual structure extracted but
the certificate direction needs the sign-convention audit; the tournaments of §3 are
definitions; nothing closes a leg. Files:
`04-computation/lrc14_two_anchor_double_cover_klein_S160.py` (+ .out). Pointers:
boxeph HYP-4801-cont (the reduction; NB: 4801 is double-filed with my S156 — flagged),
opus-S134/S137 (roof, exact engines), THM-638/THM-640, HYP-4831 (the LP), S157 (voltage
lifts), MISTAKE-116.
