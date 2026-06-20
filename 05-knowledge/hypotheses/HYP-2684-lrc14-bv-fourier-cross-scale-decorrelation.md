---
id: HYP-2684
title: LRC(14) BV-Fourier cross-scale decorrelation lemma
status: OPEN; abstract resonance-filter lemma identified, LRC coverage-function variation budget still missing
source: codex-2026-06-20-S55
depends_on:
  - HYP-2675
  - THM-546
  - THM-548
  - HYP-2682
  - HYP-2681
  - HYP-2680
related:
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
`sum_{s>=1} 1/s^2 = pi^2/6`.  The missing LRC work is not this algebraic
identity; it is to prove the required mixed-variation/Fourier budget for the
actual seven-sector coverage functions and make the constants explicit.

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

The finite decorrelated comparison is no longer the obstacle.  KPS S19 and the
S53 audit show

```text
sup p0_decorr = Q(k-1) < cap_k
```

for `k=8..12`, with the maximum at a consecutive bounded base plus one
independent stranger.

So the next sharp target is:

```text
p0(E) <= p0_decorr(E_cluster_model) + explicit Weyl/BV error.
```

If the scale gap is `G` and the mixed-variation budget is `V_mix`, the
nonresonant part should close once

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

1. Define the exact cluster coverage functions `H_C` for the HYP-2675
   decorrelated model.
2. Prove a finite mixed-variation bound for `H_C` in terms of sector-wall
   counts, missed-sector profile size, or the existing `V(E)`/state-word
   ledgers.
3. Derive the explicit two-cluster threshold `G(k,C)` from
   `V_mix/(12G) < cap_k-Q(k-1)`.
4. Extend the resonance identity to `r` clusters and record the exact
   low-height relation cutoff that sends rows into HYP-2682/HYP-2676 finite
   atlases.
5. Glue the remaining bounded-gap region to the finite checks already tracked
   by HYP-2675 and OPEN-Q-108.

No LRC(14) proof is claimed.  This hypothesis isolates the analytic lemma that
would turn the now-safe decorrelated plateau into a proof of the wide branch.
