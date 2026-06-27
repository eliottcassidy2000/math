# The Lee-Yang zeros literally are the Bonferroni failure

*kind-pasteur, 2026-06-27-S254. A reflection on the LRC(14) covering floor.*

## The convergence

Two completely different languages, the same boundary at the same place.

The LRC(14) few-apex covering problem asks: is the joint loneliness
`meas(R-safe ∩ Q-lonely)` bounded below? The documented obstruction (HYP-3121) was
phrased in **probabilistic-combinatorial** terms: the Bonferroni / union bound
`meas(A) + meas(B) - 1` goes **negative** for few-apex sets, so the intersection is not
forced positive by inclusion-exclusion — it survives "only via quasi-independence."

Reframing the same measure as a **partition function** — loneliness = `P(M=0)` for the
danger-count `M = Σ_s 1[‖s t‖ < 1/14]`, with single-fugacity transform
`Ξ(λ) = ∫ (1-λ)^{M(t)} dt = G_M(1-λ)` — turns "is it positive?" into a **Lee-Yang
zero-freeness** question: does `Ξ` have a zero in the unit λ-disk?

The two collide exactly:

- The **single-variable Lee-Yang property holds iff roughly the union-bound regime**
  (≤ 6 danger events at probability 1/7 each, so `Σ meas(D_s) = k/7 < 1`). The apex
  Q-block (r ≤ 6) is Lee-Yang with room to spare. The 14-free R-block (≥ 7 events) is
  not — and `7/7 = 1` is precisely where Bonferroni turns vacuous.
- The **bidisk is not zero-free** because the R-block factor `Ξ(λ,0) = G_R(1-λ)` already
  has a zero inside `|λ| < 1` (concretely near `λ = -0.43 ± 0.48 i`). That interior zero
  is the same object as the negative Bonferroni number. Naive Asano contraction cannot
  reach `(1,1)` because the zero is in the way.

So the "Bonferroni failure" is not a deficiency of a crude bound that a cleverer
inequality would fix. It is a **genuine zero of the partition function sitting inside the
unit disk**, and any method whose positivity certificate is a zero-free polydisk (Asano,
single-variable Lee-Yang, cluster expansion in the worst-case dependency graph) must fail
for the same structural reason. The honest negative result and the documented honest
negative result are the *same* result, seen from two sides.

## Why the apex comb escapes

The thing that makes the Q-block (apexes `14m`) Lee-Yang is not that `r/7 < 1` — the
maximally-correlated count `{0: 1-1/k, k: 1/k}` violates Lee-Yang already at `k = 4`
despite `kp = 1`. It is that the apex danger arcs are an **equidistributed comb**: their
joint statistics are close to *independent* (the Binomial extreme has zero-free clearance
growing like `k`, the opposite end from comonotone). The `14` inside every constraint
`‖14m t‖ ≥ 1/14` forces the comb onto the `14ℤ` Bohr lattice, and equidistribution does
the rest.

This is why the whole problem reduces, cleanly, to **equidistribution of the apex comb
against the coarse R-safe set** — which is exactly where the spectrum thread (THM-546,
HYP-2840, and the sibling HYP-3129 certifying `R' ≥ 0.642`) lives. The Lee-Yang frame
does not prove the floor; it *diagnoses* the floor, telling you with surgical precision
which factor is harmless (apex, Lee-Yang) and which factor carries the difficulty
(R-block, interior zero), and therefore which mathematics (equidistribution, not
zero-freeness) must finish the job.

## The transcending pattern

When a "crude" bound fails and a "sophisticated" reframing fails *at the same threshold
with the same witness*, that coincidence is information: the threshold is real, not an
artifact of either method. The constant `1/7` (half a sector, `1/14`, doubled) and the
count `7` (events before the union saturates) are the same `7` that runs through the whole
project — the seven sectors, the `‖·‖ ≥ 1/14`, the `A002854`-flavoured arithmetic of the
prime 7. The partition function did not give a new proof here; it gave a new *reason*, and
the reason is that 7 events of mass 1/7 is exactly one full circle. Past that, you are no
longer counting — you are equidistributing.

→ HYP-3128 (the dichotomy + obstruction), HYP-3129 (the elementary equidistribution
floor), HYP-3121 (the unified covering argument), HYP-2840 / THM-546 (the comb).
