# HYP-2481 - LRC14 blocking-height dominance accumulates but normalizes away

**Status:** OPEN proof program with exact finite atlas.

**Source:** codex-2026-06-13. Extends HYP-2444, HYP-2471, and HYP-2480.

**Computation:** `04-computation/lrc14_blocking_height_dominance_codex.py`; stored output `05-knowledge/results/lrc14_blocking_height_dominance_codex.out`.

## Statement

Let `h(S)` be the dilated-band blocking height:

```text
h(S) = least q with an uncovered unit numerator in (Z/q)^*.
```

For `q<h(S)`, the row blocks every unit numerator at shell `q`.  For two
speeds `v,w in S`, define a pre-height dominance relation by comparing their
dilated-band cover masks:

```text
v -> w iff sum_{14 <= q < h(S)} |U_v(q)\U_w(q)|
          >= sum_{14 <= q < h(S)} |U_w(q)\U_v(q)|.
```

The tie Hamiltonian path is the row order.  This quotient deliberately uses
speeds as cover carriers, not runners at a time.  It preserves which speeds
carry the blocked-unit cover before the first leak, and it discards continuous
owner geometry and Bprime interval widths.

The data support a sharpened form of the prompt's question:

```text
blocking height increases raw cumulative dominance,
but height-normalized dominance thins as h grows.
```

Thus the proof route is not "high rows become more dictatorial."  It is a
dichotomy:

```text
long blocking-height row
=> peelable cumulative carrier
   OR balanced-cover congruences tight enough to force the Q31/band-2 portal.
```

## Evidence

### One-stranger family

For `S(r)=7*{1,...,12} union {r}`, with primitive `r<=1092`, the atlas tested
`936` rows.  Only `72` have pre-height shells above the `q_start=14` threshold.

Correlations with height:

```text
corr(height, mean_pair_margin)            =  0.779 pearson,  0.522 spearman
corr(height, mean_pair_margin_per_shell)  = -0.410 pearson, -0.290 spearman
corr(height, mean_pair_margin_norm)       = -0.711 pearson, -0.558 spearman
corr(height, top_unique_share)            = -0.492 pearson, -0.497 spearman
```

The height buckets show the raw accumulation:

```text
h=27: n=64, mean pair margin  6.83
h=40: n= 6, mean pair margin 11.46
h=41: n= 2, mean pair margin 10.62
```

But the per-shell and per-unit versions move downward.  Higher blocking height
spreads cover work more evenly per obligation while accumulating more total
cover history.

### Random primitive rows

In `120` sampled primitive rows containing a multiple of `14`, `23` have
pre-height shells.  The same pattern is stronger:

```text
corr(height, mean_pair_margin)            =  0.942 pearson,  0.896 spearman
corr(height, mean_pair_margin_per_shell)  = -0.573 pearson, -0.725 spearman
corr(height, mean_pair_margin_norm)       = -0.729 pearson, -0.875 spearman
corr(height, top_cover_share)             = -0.581 pearson, -0.780 spearman
```

So raw dominance growth is not a one-stranger artifact.  The normalization
failure is also not a one-stranger artifact.

### Named hard packets

The five one-stranger evaders from the dilated-band probe:

```text
r in {611,702,793,962,1053}
```

have `h=40` or `h=41`.  The HYP-2470/HYP-2471 eight-core exception packets
have `h=31` or `h=33` in this atlas.

Every named packet's speed-dominance tournament is transitive:

```text
score_hist={0:1,1:1,...,12:1}
score_spread=12
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

This is not a sign of rich tournament cycling.  It is a sign that this binary
speed quotient is too coarse once cumulative cover masks are summed.  The
useful information lives in the margin scale, the unique-load ledger, and the
failure of normalized concentration.

## Relation To HYP-2471 And HYP-2480

HYP-2471 says the two HYP-2470 Q27 exceptions do not persist once the fiber
ladder is widened to Q31.  HYP-2480 explains the same exceptions as ramified
7-ideal / 13-clock packets, analogous to an Eisenstein/Newton local gate.

HYP-2481 supplies the dominance/load side of that picture.  Before the first
unit leak, high-blocking rows build a long perfect-cover streak.  The streak
has cumulative dominance, but it is not increasingly concentrated per unit.
This explains why Q27 can look closed while Q31 or band-2 shells expose the
row: the scalar residue face is balanced enough to hide the leak until the
next fiber layer, but the accumulated cover margins still identify carriers
and deletion addresses.

In the one-stranger hard rows, top cover carriers include `35`, `70`, and
sometimes the stranger; top unique-load carriers often include small core
speeds such as `7` and `49`.  That split is proof-relevant: one object carries
coverage mass, another carries private obligations.  A proof must retain both.

## Tournament Analysis

Vertices considered before choosing this quotient: runners, gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, deleted-core addresses, denominator shells,
unit numerators, owner/Bprime exits, and proof obligations.

Selected vertices: speeds inside a fixed row.

Pairwise observable:

```text
sum_q |U_v(q)\U_w(q)| - |U_w(q)\U_v(q)|, 14 <= q < h(S).
```

Switch/gauge: positive cumulative margin or row-order tie.

Tie Hamiltonian path: row order.

Preserves: cover-load dominance over fully blocked pre-height shells.

Destroys: continuous-time owner geometry, Bprime widths, and the internal
nontransitivity of unit obligations.

Challenged assumption: tournaments need not use runners as vertices.  Here
even using speeds is probably still too coarse; unit-obligation or
shell-transition vertices may be the next quotient if we want directed cycles
rather than a transitive dominance pipeline.

## Proof Program

1. Prove a peelable-carrier lemma: if a speed has enough cumulative excess or
   private load over `q<h`, then deleting/transporting that carrier exposes a
   witness, Bprime opening, or lower-core descent.
2. Prove a balanced-cover rigidity lemma: if no peelable carrier exists, then
   the cover masks must satisfy congruence/valuation constraints that force
   the ramified Q31 or band-2 portal.
3. Add leave-one-out support-criticality to the atlas: for each carrier, record
   which pre-height shells fail after deletion and how the first-leak deficit
   changes.
4. Switch tournament vertices from speeds to unit obligations, denominator
   shells, or cover-obligation types to recover nontransitive structure lost by
   the cumulative speed quotient.
5. Use this dominance/balance dichotomy as the typed budget for the next
   below-eight-core `e=5` LRC14 search.
