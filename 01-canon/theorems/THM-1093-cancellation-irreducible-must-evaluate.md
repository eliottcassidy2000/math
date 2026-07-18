---
id: THM-1093
title: THE CANCELLATION IS IRREDUCIBLY ACROSS COSETS — NO TRIANGLE-INEQUALITY BOUND CAN CAPTURE IT, SO δ(S) MUST BE EVALUATED RATHER THAN ESTIMATED — since ĥ(n) = sin(πn/7)/(πn) and sin(πn/7) depends only on n mod 14, δ splits exactly as (1/π^k)Σ_r ∏c(rᵢ)·T(r) over the residue subgroup Λ mod 14 (rebuild verified to 6 digits); measuring the three granularities gives COSET-ABS ≈ FULL-ABS (1.64 vs 1.86, 12.83 vs 12.95, …), so keeping cancellation INSIDE cosets gains essentially nothing and all of it lives ACROSS cosets. Two further facts, both by the antipodal symmetry r ↔ −r: the character sum C(A) = Σ_r ∏c(rᵢ) is IDENTICALLY ZERO for odd k (hence no discriminant — my "leading coefficient" framing was wrong), and T̄ = 0 identically, so the variation bound Σ|∏c|·|T(r)−T̄| EQUALS the coset-absolute bound exactly and gains nothing. Conclusion: every bound that decomposes the sum and applies a triangle inequality at coset granularity is provably lossy by the same 12× spread; the sum must be evaluated in closed form, i.e. as a higher-dimensional Dedekind sum
status: the coset decomposition verified exactly (rebuilt = signed to 6 digits, 6 families); the three-granularity comparison measured over 12 families at k=3; C ≡ 0 and T̄ = 0 are proved by antipodal symmetry, not merely observed; the closed-form identification for k≥3 is a NAMED TARGET, not a result
source: opus-2026-07-17-S371 (owner: work the bound that keeps the signs)
depends_on: [THM-1085 (the absolute bound whose non-uniformity prompted this), THM-1080 (m₇), THM-1075 (the resonance lattice), THM-1070 (the first instance of the same failure mode), THM-965 (the k=2 evaluation, which IS the Bernoulli/fold closed form)]
scripts: 04-computation/sign_cancellation_opus_S371.py, character_sum_opus_S371.py, variation_bound_opus_S371.py -> 05-knowledge/results/
---

# THM-1093 — why "keep the signs" cannot be done by bounding

THM-1070 and THM-1085 both failed the same way: a valid bound whose
looseness grows with k, in each case because signs were discarded. This
file asks whether the signs can be kept by a finer decomposition. They
cannot, and the reason is structural.

## The exact splitting

ĥ(n) = sin(πn/7)/(πn), and sin(πn/7) depends only on n mod 14. Writing
c(r) = sin(πr/7):

> **δ(S) = (1/π^k) Σ_{r ∈ Λ mod 14} ∏ᵢ c(rᵢ) · T(r)**,
> T(r) = Σ over {n ∈ Λ : n ≡ r mod 14} of 1/∏nᵢ

The residues form a subgroup of (ℤ/14)^k. Rebuild verified against the
direct sum to 6 digits on all families tested. All signs sit in two
places: the finite table ∏c(rᵢ), and 1/∏nᵢ inside each coset.

## Where the cancellation lives — measured

| family | \|δ\| | COSET-ABS | FULL-ABS |
|---|---|---|---|
| (2,3,5) | 0.019412 | 0.031912 | 0.036197 |
| (31,37,41) | 0.000829 | 0.002401 | 0.002410 |
| (6,10,15) | 0.000334 | 0.004282 | 0.004321 |
| (101,103,107) | 0.003873 | 0.004239 | 0.004283 |

COSET-ABS (which keeps cancellation *inside* each coset) is within a few
percent of FULL-ABS (which keeps none). **Essentially all cancellation is
across cosets.**

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

## What this rules out, and what remains

Any bound obtained by decomposing δ into coset pieces and applying a
triangle inequality is stuck at the same ~12× spread at k=3, because the
cancellation it must capture is precisely the *between-piece* cancellation
that the triangle inequality destroys. Finer decomposition does not help:
cosets are already the finest granularity on which ∏c(rᵢ) is constant.

**So δ(S) must be EVALUATED, not estimated.** This is exactly what
happens at k=2, and it is worth seeing why that case closed:

> δ(a,b) = (1/(π²a'b')) Σ_{s≠0} sin(πsa'/7)·sin(πsb'/7)/s²

which is the Fourier series of the second Bernoulli polynomial, summing
to the tent function fold_M(r) = r(M−r) — i.e. THM-965. It was never
bounded; it was summed.

The k ≥ 3 analogue, Σ over a rank-(k−1) lattice of ∏sin(πnᵢ/7)/∏nᵢ, is a
**higher-dimensional Dedekind sum** in Zagier's sense. Those carry
reciprocity laws, which is the mechanism that makes such sums computable
rather than merely bounded. That is the named target, and it is now the
only route this program has not ruled out.

**Three sessions, three approaches, one verdict.** Containment/
fragmentation (THM-1070), absolute values (THM-1085), and coset
decomposition (here) all lose the same way. The consistency of the
failure is itself the finding: the quantity is not estimable at this
granularity, and the remaining work is a closed-form evaluation.
