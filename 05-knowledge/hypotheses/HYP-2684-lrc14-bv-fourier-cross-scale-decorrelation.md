---
id: HYP-2684
title: LRC(14) BV-Fourier cross-scale decorrelation lemma
status: PARTIAL. THM-2054 proves the abstract nonresonant filter, the exact seven-sector signed-atom budget, and the two-scale alias cutoff |M|>H sum|b_i|. At H=2^20 its uniform factorwise error is <90816/1048577; rowwise, H=2^19 is already below each recorded pinned-base cap-Q margin. The remaining gates are correct lifted-model identification, a model-specific decorrelated plateau bound, and bounded-resonance classification; MISTAKE-080/082 forbid importing Q(k-1) from cardinality alone. LRC(14) is not proved.
source: codex-2026-06-20-S55
depends_on:
  - HYP-2675
  - THM-546
  - THM-548
  - HYP-2682
  - HYP-2681
  - HYP-2680
related:
  - THM-2054
  - HYP-2683
  - HYP-2679
  - HYP-2677
  - HYP-2676
  - HYP-2639
  - HYP-2637
  - OPEN-Q-108
external:
  - https://arxiv.org/pdf/2301.05561
  - https://www.math.bas.bg/infres/MathBalk/MB-19/MB-19-349-366.pdf
  - https://arxiv.org/abs/1106.4673
  - https://arxiv.org/pdf/0902.1717
---

# HYP-2684 - BV-Fourier Cross-Scale Decorrelation

## 2026-07-21 resolution of the abstract nonresonant step

THM-2054 proves a stronger factorwise substitute for the mixed-variation
target below. For characters `chi_i`, an integer line `lambda`, and degree-`H`
Fejer approximants, assume every height-`H` scalar relation along the line is
already a vector relation among the characters. Then

```text
|line product average-full-torus product average|
 <=2 sum_i ||f_i-Fejer_H*f_i||_1
              product_(j!=i)||f_j||_infinity.          (R)
```

The same holds after summing any finite signed Boolean atom expansion. For
two-scale speeds `e_i=b_i+M c_i`, the resonance implication is automatic once

```text
|M|>H sum_i |b_i|.
```

Thus the abstract nonresonant/Weyl lemma is no longer open, and the actual LRC
coverage algebra does not require a mixed-variation theorem: expand it into
finitely many products of interval predicates and apply (R). THM-2054 now
also performs that expansion for THM-534's six inner sectors. After deleting
the pinned zero factor, `sum_A |A|=192`, so for `r<=11` active offsets and
`H=2^20` the complete signed error is

```text
<90816/1048577.                                       (R1)
```

This is below the smallest recorded `cap-Q` numerical margin by

```text
46851988331/1028633065460>0.                          (R2)
```

Using the five margins row by row, the same theorem permits `H=2^19`; the
tightest residual is the `k=9` value
`2064067449/171439007740>0`. This halves the later two-scale alias cutoff.

However, (R2) is not yet an LRC comparison. MISTAKE-080/082 show that a
decorrelated limit and the pinned-base `Q(k-1)` model cannot be imported as a
majorant for an arbitrary cluster shape; the exact row `[0,19,...,25]`
already exceeds `Q(7)`. What remains is to identify the full-torus model
produced by the chosen character lift, prove its correct plateau bound, and
route the finite bounded-resonance branch. The older mixed-BV formulation is
retained below as a valid alternative and historical derivation.

Namespace note: an incoming concurrent session had already claimed
`HYP-2683/T922` for wide-branch address repair.  This BV-Fourier route was
renumbered to `HYP-2684/T923`.  The two notes are complementary: HYP-2684
describes the analytic decorrelation error, while HYP-2683 records address data
that may route the low-height finite-resonant branch.

## Claim Being Tested

The surviving HYP-2675 route should be written as a Fourier resonance filter,
not as another scalar discrepancy constant.

For a two-scale cluster decomposition, let `H(x,phi)` be the exact Boolean
coverage function of the bounded/slow variable `x` and the fast cluster anchor
phase `phi`.  The actual row with scale separation `M` is

```text
I_M(H) = integral_T H(x, Mx) dx,
```

while the fully decorrelated model is

```text
J(H) = integral_{T^2} H(x, phi) dx dphi.
```

With Fourier convention `e(t)=exp(2*pi*i*t)`,

```text
I_M(H) - J(H) = sum_{s != 0} Hhat(-M*s, s).
```

Thus decorrelation is not mystical randomness.  It is the statement that the
line `r+M*s=0` samples only high-frequency coefficients of `H`, unless a
low-height resonance survives.

The target lemma is:

```text
if |Hhat(r,s)| <= V_mix(H)/(4*pi^2*|r*s|) for r*s != 0,
then |I_M(H)-J(H)| <= V_mix(H)/(12*M).
```

The proof is immediate from the displayed resonance identity and
`sum_{s>=1} 1/s^2 = pi^2/6`. Historically the proposed next step was to prove
the required mixed-variation budget for the seven-sector coverage functions.
THM-2054 shows that this is optional for the actual finite product algebra;
its factorwise coefficient budget is now the more direct target.

## Why This Came From The Web Search

The lacunary-sequence literature is the right analogy.  Aistleitner, Berkes,
and Tichy describe sums `f(n_k x)` with Hadamard or super-lacunary gaps, with
`f` often of bounded variation, and emphasize that the dependence structure is
controlled jointly by the analytic regularity of `f` and the arithmetic of the
frequency sequence.  That is exactly the LRC wide branch: analytic regularity is
the BV/sector-wall budget; arithmetic is the small relation lattice among
cluster anchors.

The Erdos-Turan-Koksma inequality expresses discrepancy through trigonometric
sums, and Koksma-Hlawka variants for piecewise smooth indicators explain why
piecewise sector-coverage functions should be handled by variation/discrepancy,
not by smooth majorants alone.  Berndt's Fourier-transform paper records the
classical input that bounded variation gives `O(1/n)` Fourier decay; THM-546 is
the repo's one-far exact implementation of that principle.

## LRC(14) Proof Implication

For the pinned-base models actually audited by KPS S19 and S53, the proposed
finite decorrelated comparison is

```text
sup p0_decorr = Q(k-1) < cap_k
```

for `k=8..12`, with the maximum at a consecutive bounded base plus one
independent stranger.

This value is a model-specific decorrelated limit, not a universal majorant.
The next sharp target is therefore:

```text
p0(E) <= p0_lift(E_cluster_model) + explicit Fejer error,
followed by a proved model-specific bound on p0_lift.               (P)
```

THM-2054 closes the error term in (P) for height-separated character lifts.
It does not identify `p0_lift` with `Q(k-1)`. If that identification and bound
are proved for the intended pinned-base model, then the `H=2^20` estimate
(R1)--(R2) closes its nonresonant comparison. For other shapes the full-torus
plateau must be bounded directly.

The historical mixed-BV alternative would close once

```text
V_mix/(12*G) < cap_k - Q(k-1).
```

This gives a concrete constant-producing route.  It may be huge, but any
explicit `G` converts the true-wide branch to finite bounded-gap glue.

## Resonance Dichotomy

For more than two clusters, write the coverage function on a torus of anchor
phases.  Fourier expansion gives the exact surviving terms

```text
sum_i n_i M_i = 0.
```

The desired split is:

1. **No low-height relation:** every nonzero resonance has a large coefficient,
   so BV decay pays a factor `1/G`.
2. **Rank-one low-height relation:** route by HYP-2682/HYP-2681 phase/support
   atlases, especially AP triples like `u-2v+w=0`.
3. **Higher-rank low-height relation:** route by HYP-2676/HYP-2677
   Ruzsa/Freiman and packet-tournament finite models.
4. **Scale-invariant `d=1` rows:** normalize by THM-531 and use the finite
   true-wide ledgers already identified in HYP-2675/HYP-2679.

This is the promised bridge between Weyl/decorrelation and additive energy:
high additive energy is exactly the branch where low-height relations survive
the Fourier filter and must be modeled finitely.

## Challenged Assumption

Do not assume the tournament vertices are runners.  In this proof route the
natural vertices include:

- scale clusters;
- Fourier modes;
- primitive resonance equations;
- missed-sector profile atoms;
- sector-wall crossing events;
- packet-sign obligations;
- finite proof ledgers.

The quotient to cluster anchors preserves the decorrelation predicate and the
small-relation obstruction.  It destroys wall ownership and phase.  That loss is
why HYP-2682's AP/cube-root phase labels and HYP-2677's packet tournament labels
must be retained before scalarizing.

## Tournament Analysis

Candidate tournament vertices: scale clusters, not runners.

Pairwise observable:

```text
C_ij = the explicit BV covariance budget for clusters i,j after removing
       all low-height relations already routed to finite atlases.
```

Switch/gauge: orient `i -> j` when cluster `j` is the larger scale after
normalization and the nonresonant `1/G_ij` bound is smaller than the available
plateau margin.  If the edge is tied or fails, mark it as a low-height relation
edge and route to the finite resonance atlas instead of forcing an orientation.

Tie Hamiltonian path: scale order after collapsing exact `d=1` dilates.  A
directed cycle is not a contradiction; it is a signal that the cluster anchors
are participating in a rank-one or higher-rank relation and should be handled
by HYP-2682/HYP-2676 rather than by the nonresonant Weyl estimate.

## Next Work Items

1. For each intended cluster lift, identify the exact full-torus atom model,
   including which cluster contains the pinned zero, and prove its own
   plateau/cap bound. Do not reuse `Q(k-1)` by cardinality.
2. Use THM-2054's proved rowwise `H=2^19` error budget and record the explicit
   scale cutoff `|M|>2^19 sum|b_i|` wherever the lifted model is one of the
   correctly identified pinned-base models; otherwise recompute its plateau
   margin before selecting `H`.
3. Route every failed cutoff/low-height resonance into HYP-2682/HYP-2676 or
   the THM-2052/2053 finite relation-plane atlas.
4. Glue the remaining bounded-gap region to the finite checks already tracked
   by HYP-2675 and OPEN-Q-108.

No LRC(14) proof is claimed. The abstract analytic lemma and exact signed atom
budget are now THM-2054; this hypothesis retains model identification,
model-specific plateau control, and resonant glue.
