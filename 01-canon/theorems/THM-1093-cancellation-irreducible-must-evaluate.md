---
id: THM-1093
title: THE RESIDUE-COSET SPLITTING IS EXACT, BUT ITS DIRECT TRIANGLE BOUND PRESERVES THE OBSERVED SPREAD — absolute convergence permits the exact mod-14 partition; for odd support the character sum and mean reciprocal coset sum both vanish by antipodal symmetry; on twelve k=3 families COSET-ABS is close to FULL-ABS and varies by 12.4x, so this specific refinement gains essentially nothing and points toward generalized Dedekind-sum evaluation
status: the coset decomposition is rigorous by THM-1092 absolute convergence (and numerically rebuilds the signed sum to 6 digits on six families); C ≡ 0 and Tbar = 0 for odd support are proved by antipodal symmetry; the three-granularity comparison is measured over twelve k=3 families and does not rule out every possible signed estimate; a closed-form higher-dimensional evaluation is a named target, not a result
source: opus-2026-07-17-S371 (owner: work the bound that keeps the signs)
depends_on: [THM-1092 (absolute convergence and legal regrouping), THM-1085 (the absolute bound whose variation prompted this), THM-1080 (m₇), THM-1075 (the resonance lattice), THM-1070 (the first instance of the same failure mode), THM-965 (the k=2 evaluation, which IS the Bernoulli/fold closed form)]
scripts: 04-computation/sign_cancellation_opus_S371.py, character_sum_opus_S371.py, variation_bound_opus_S371.py -> 05-knowledge/results/
---

# THM-1093 — the direct residue-coset triangle bound does not sharpen the sample

THM-1070 and THM-1085 both lost substantial information when signs were
discarded.  This file tests one precise repair: first retain cancellation
inside each residue coset, then apply a triangle inequality across cosets.
The exact symmetry simplifications are rigorous; the conclusion that the
resulting bound is still loose is a twelve-family `k=3` measurement.

## The exact splitting

ĥ(n) = sin(πn/7)/(πn), and sin(πn/7) depends only on n mod 14. Writing
c(r) = sin(πr/7):

> **δ(S) = (1/π^k) Σ_{r ∈ Λ mod 14} ∏ᵢ c(rᵢ) · T(r)**,
> T(r) = Σ over {n ∈ Λ : n ≡ r mod 14} of 1/∏nᵢ

The residues form a subgroup of (ℤ/14)^k.  THM-1092 proves the full-support
reciprocal lattice series converges absolutely, so this finite regrouping is
legal and every `T(r)` converges absolutely.  The numerical rebuild agrees
with the direct sum to 6 digits on all families tested. All signs sit in two
places: the finite table ∏c(rᵢ), and 1/∏nᵢ inside each coset.

## Where the cancellation lives — measured

| family | \|δ\| | COSET-ABS | FULL-ABS |
|---|---|---|---|
| (2,3,5) | 0.019412 | 0.031912 | 0.036197 |
| (31,37,41) | 0.000829 | 0.002401 | 0.002410 |
| (6,10,15) | 0.000334 | 0.004282 | 0.004321 |
| (101,103,107) | 0.003873 | 0.004239 | 0.004283 |

COSET-ABS (which keeps cancellation *inside* each coset) is within a few
percent of FULL-ABS (which keeps none). **In these sampled families,
essentially all observed cancellation is across cosets.**

## Two symmetry facts that close the route

Both follow from Λ being antipodally symmetric, with c odd and
1/∏(−nᵢ) = (−1)^k/∏nᵢ:

1. **C(A) = Σ_r ∏c(rᵢ) ≡ 0 for odd k.** Computed as 0.0000 for all
   twelve families tested. This is a triviality of the symmetry, *not* a
   discriminant between families — my framing of C as a leading
   coefficient whose vanishing would signal small δ was **wrong**, and
   the computation says so unambiguously: every family gives zero,
   including those with the largest and smallest deviations.

2. **T̄ = 0 for odd k**, by the same pairing. Hence the "variation bound"
   Σ_r |∏c(rᵢ)|·|T(r) − T̄| — motivated by the mean dropping out — is
   **exactly** the coset-absolute bound. Verified: the two columns agree
   to six digits on all twelve families. The reduction was vacuous.

Spreads: VAR-BOUND 1.03–12.83 (12.4×), FULL-ABS 1.10–12.95 (11.8×). The
refined bound is not better; it is the same bound.

## What this particular refinement rules out, and what remains

For these twelve families, the direct coset-wise triangle bound is stuck at
the measured ~12× spread because it discards the *between-coset*
cancellation.  Refining a coset and then summing absolute values cannot
improve its contribution: another triangle inequality can only increase it.
This rules out that specific class of refinements as a sharp certificate on
the sample.  It does **not** prove that every analytic signed estimate fails.

The evidence therefore points toward **evaluating δ(S), rather than applying
this triangle bound**. This is exactly what happens at k=2:

> δ(a,b) = (1/(π²a'b')) Σ_{s≠0} sin(πsa'/7)·sin(πsb'/7)/s²

which is the Fourier series of the second Bernoulli polynomial, summing
to the tent function fold_M(r) = r(M−r) — i.e. THM-965. It was never
bounded; it was summed.

The k ≥ 3 analogue, Σ over a rank-(k−1) lattice of ∏sin(πnᵢ/7)/∏nᵢ, is a
natural **generalized higher-dimensional Dedekind-sum** target. Reciprocity
laws are the mechanism that makes such sums computable rather than merely
bounded. Establishing the exact match and extracting a usable formula remain
open.

Containment/fragmentation (THM-1070), full absolute values (THM-1085), and
this coset decomposition all lose information in the tested regime.  The
consistent signal is useful but narrower than a no-go theorem: direct
triangle inequalities at these granularities are poor, while closed-form
evaluation and genuinely coupled signed estimates remain live.
