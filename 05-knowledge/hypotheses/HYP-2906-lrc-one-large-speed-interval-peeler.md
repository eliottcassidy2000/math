---
id: HYP-2906
status: PROOF FRAGMENT / sharp one-large-speed reduction
source: codex-2026-06-22-S117
tags: [lrc14, induction, scale-separation, interval-peeler, node3, open-q-108]
depends_on:
  - HYP-2904
  - HYP-2900
  - HYP-2901
  - HYP-2902
  - HYP-2895
related:
  - HYP-2905
  - HYP-2903
  - THM-526
  - THM-527
  - THM-566
  - OPEN-Q-108
results:
  - 04-computation/lrc_one_large_interval_peel_codex_s117.py
  - 05-knowledge/results/lrc_one_large_interval_peel_codex_s117.out
---

# HYP-2906: a one-large-speed interval peeler gives the ratio `n-1`

Workspace convention: `LRC(n)` means `n-1` speeds and target `1/n`.

HYP-2904 proves a global finite-comb induction lemma using the full
threshold-`1/n` safe set of the seed and its number of interval components.
There is a sharper local atom.  One connected seed-safe interval is enough.

## Lemma

Let `B` be a seed speed set with maximum speed `m`.  Suppose there is a time
`tau` and a margin `alpha>1/n` such that

```text
||b tau|| >= alpha       for every b in B.
```

Then `B union {v}` is LRC-safe at threshold `1/n` whenever

```text
v > m / (n*(alpha - 1/n)).
```

Proof.  Put `rho=(alpha-1/n)/m`.  For every `t` in the connected circle arc
`I={|t-tau|<=rho}`, every seed runner remains safe at threshold `1/n`, because

```text
||b t|| >= ||b tau|| - b|t-tau|
          >= alpha - m rho
          = 1/n.
```

The unsafe set of the added speed `v` is the comb

```text
U_v={t: ||v t|| < 1/n},
```

whose connected components have length `2/(n v)`.  If `I` were contained in
`U_v`, connectedness would force `I` to lie inside one comb tooth.  But

```text
length(I)=2(alpha-1/n)/m > 2/(n v),
```

exactly when the displayed threshold holds.  Therefore some `t in I` avoids
`U_v`, and all runners in `B union {v}` are safe.

## LRC14 corollary

Taking `alpha=1/(n-1)` from `LRC(n-1)` gives the clean symbolic threshold

```text
v > (n-1)m.
```

Thus, assuming the known smaller case `LRC(13)`, any LRC14 row with largest
speed `v` and second-largest speed `m` is safe as soon as

```text
v > 13m.
```

Equivalently, every LRC14 counterexample must be top-balanced:

```text
v_max <= 13 v_second.
```

This is much sharper than the global component-budget bound when only
existence is required.  It also explains the earlier arc-width thread
(THM-526): LRC is an existence problem, and a single arc wider than one comb
tooth produces a witness without describing its final denominator.

## Exact constants

The audit script checks the algebra exactly.

For a generic LRC14 seed with `m=13`, the seed witness supplied by LRC13 has
`alpha=1/13`, so the guaranteed seed-safe interval has length `1/1183`.  The
new danger tooth at `v=170` has length `1/1190`, so every `v>=170` is certified.

For the AP-core seed

```text
{1,2,3,4,5,6,7,8,9,10,11,13},
```

the explicit witness `tau=1/12` has `alpha=1/12` and `m=13`.  The threshold
improves to

```text
v > 78,
```

so `v>=79` is certified.  This covers the committed radical/lcm speeds
`30030`, `60060`, and `510510` with a local existence argument.  HYP-2904
certified the same AP-core family only at `v>=768`, because it was proving a
positive-measure lower bound for the whole safe set rather than one witness.

For the pure dilation hard core `{b,2b,...,12b,V}`, the generic smaller-LRC
witness gives `V>156b`.  The user's sheet-counting proof remains stronger in
the comparable regime; HYP-2906 isolates the high-ratio part without any
equidistribution or sheet counting.

## Proof impact

This inserts a strict reduction before the analytic Node-3 machinery:

```text
top speed > 13 * second speed
  -> peel by LRC13 and one interval tooth comparison
  -> LRC14 safe.
```

The remaining unproved region is therefore not arbitrary unbounded scale
separation.  It is either:

1. top-balanced (`v_max<=13 v_second`), where bounded-core / AP-hull /
   Legendre-Venn Node-2 structure has to do real work; or
2. multi-large with no single speed beating the local ratio gate, where
   HYP-2904/KPS-S31v's component or arc-count budgets and the resonant-pair
   second moment are still needed.

This also clarifies HYP-2903's routing.  Positive high-depth Bonferroni tails
do not threaten the one-large branch once a far speed is locally peelable:
the LRC witness is produced before the p0-cap obligation is invoked.  The
missing-depth parity guard is needed in the balanced/binding cap leg, not in
this scale-separated existence leg.

## Tournament Analysis

Vertices are proof carriers:

```text
connected_seed_safe_arc
single_comb_tooth_width
smaller_LRC_margin
scale_ratio_gate
global_component_budget
runner_count_only_induction
raw_runner_vertices
```

The Hamiltonian path is transitive in the audit.  The challenged assumption is
that induction vertices are runner sets or even whole safe sets.  For the
one-large branch, the preserved predicate is a connected safe interval whose
length beats one tooth of the added runner's danger comb.  The quotient
destroys global safe-set measure and component topology, which is why it proves
existence but not a uniform witness-density floor.
