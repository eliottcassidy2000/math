---
id: HYP-3458
title: LRC14 AP84 coloring-recursion bridge
status: EVIDENCE / exact coloring-recursion sidecar for AP-tail bridge; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-2247/HYP-2243 coloring recursion and HYP-3454/HYP-3456/HYP-3457 AP84 clocks
tangent: T1418
technique: LTI-418
tournament_technique: LTT-318
script: 04-computation/lrc14_ap84_coloring_recursion_bridge_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_coloring_recursion_bridge_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-coloring-recursion-bridge-codex-20260629.md
related:
  - HYP-3457
  - HYP-3456
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3439
  - HYP-3438
  - HYP-3436
  - HYP-3431
  - HYP-2247
  - HYP-2246
  - HYP-2243
  - HYP-2241
  - THM-523
  - OPEN-Q-108
---

# HYP-3458: LRC14 AP84 Coloring-Recursion Bridge

## Claim

The recent AP84 bridge has become a coloring recursion in the sense of
HYP-2247, but a tame finite one rather than a Paris-Harrington fast-growing
one.

For

```text
S_m = {1,2,...,11,13,84m}
```

HYP-3456 gives the number of high `E:84m` safe gaps meeting the left low
corridor `C1=[8/49,6/35]`:

```text
N(m)=floor((504m-6)/70)-floor((96m-13)/14).
```

HYP-3458 views this as a `35`-state coloring after subtracting the Beatty
clock `floor(12m/35)`.  The boundary color is

```text
d_r = N(r)-floor(12r/35),  r=1..35.
```

It then adds the missing endpoint-address subcolor.  The HYP-3433/HYP-3454
selected address

```text
a_m = ceil(48m/7)
```

is always the first or second high gap in `C1`:

```text
rank_C1(a_m)=1  iff  7 | m
rank_C1(a_m)=2  otherwise.
```

Thus the AP-tail clock is not just a count.  It is a small recursive coloring
state:

```text
phase color        mixed for m=1..4, pure for m>=5
boundary color     d_r in {0,1,2} on residues mod 35
endpoint subcolor  rank_C1(a_m) in {1,2} on residues mod 7
```

This connects directly to HYP-2247's lesson: do not count colorings without
their extension-rank/address sidecar.

## Endpoint-Rank Lemma

The endpoint-rank subcolor has a short exact proof.  Let

```text
L_m=floor((96m-13)/14)+1
```

be the first high-gap index meeting `C1`.  Write `m=7q+s` with `0<=s<7`.
Then

```text
a_m=48q+ceil(48s/7)
L_m=48q+floor((96s-13)/14)+1
rank_C1(a_m)=a_m-L_m+1
              =ceil(48s/7)-floor((96s-13)/14).
```

The residue table is

```text
s:    0 1 2 3 4 5 6
rank: 1 2 2 2 2 2 2
```

so the endpoint witness is first exactly on the `7|m` column and second on
all other columns.  This is the concrete address-rank derivative that HYP-2247
was asking us to retain.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_coloring_recursion_bridge_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_coloring_recursion_bridge_codex_20260629.out
```

Validation through `m=350`:

```text
period_failures=[]
rank_failures=[]
phase_failures=[]
address_failures=[]
```

Boundary color word, residues `1..35`:

```text
[2,2,1,1,1,1,1,2,1,1,2,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,1,2,1,2,2,1,1,0]
```

Boundary color histogram:

```text
{0:1, 1:22, 2:12}
```

Endpoint-rank subcolor word:

```text
[2,2,2,2,2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1]
```

Endpoint-rank histogram:

```text
{1:5, 2:30}
```

Transition increments around the `35`-state cycle are

```text
[0,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,2],
```

with histogram `{0:24,1:10,2:1}` and total shift `12` per corridor, hence
`24` total escapes after mirror doubling.

## CRT Color Grid

The `35` states are better read as a `5 x 7` CRT grid.  Rows are `m mod 5`;
columns are `m mod 7`.  Entries are

```text
boundary_color / endpoint_rank_color.
```

```text
mod5=0: 0/1 1/2 1/2 1/2 1/2 1/2 1/2
mod5=1: 1/1 2/2 2/2 2/2 2/2 2/2 1/2
mod5=2: 1/1 2/2 2/2 2/2 2/2 1/2 1/2
mod5=3: 1/1 2/2 2/2 1/2 1/2 1/2 1/2
mod5=4: 1/1 2/2 1/2 1/2 1/2 1/2 1/2
```

This grid explains two older coloring echoes without overclaiming them:

1. The `mod 35 = 5*7` clock is an outer-extension coloring of the AP-tail
   boundary count.
2. The `mod 7` endpoint subcolor locates the canonical HYP-3433 homogeneous
   trace among the first two left-corridor high gaps.

## Paris-Harrington Reading

HYP-2247 says a coloring proof should retain its bad-child extension rank, not
only the number of colorings.  HYP-3458 gives the AP84 version:

```text
bad_phase_counts={q0_mixed_residues:4, q1_mixed_residues:0}.
```

The mixed endpoint color is a finite bad branch: it appears only for `m=1..4`
and disappears after the outer extension `m -> m+35`.  The pure tail is a
`35`-state automaton with retained colors, not a fast-growing PH tree.

This is useful because it gives a template for the general LRC14 proof:

```text
bad component/gate color
  -> retained endpoint/owner address
  -> outer-extension transition rank
  -> proof that bad branches die or route named debt.
```

## Proof Pull

HYP-3454, HYP-3456, and HYP-3457 already close the AP-tail endpoint interval,
escape count, and finite transients.  HYP-3458 adds the coloring-recursion
sidecar that says what a legal quotient must retain when those packets are
spliced into HYP-3439:

```text
phase color + mod-35 boundary color + mod-7 endpoint-rank subcolor.
```

The next useful proof move is to express the HYP-3438 gate law in the same
language.  HYP-3438 already says outer `(13,7)` gates are controlled by `7|m`,
inner `(11,5)` gates by `5|m`, and `35|m` leaves clean components.  HYP-3458
is the companion endpoint-address coloring: the canonical witness is first
gap exactly on the `7|m` column and second gap otherwise.

## Tournament Analysis

Tournament vertices are coloring/recursion proof carriers, not runners or raw
colors.

```text
score_hist={9:1,21:1,46:1,50:1,54:1,55:1,57:1,58:1}
directed_3cycles=0
hamiltonian_path=
  crt_35_coloring_state
  -> endpoint_rank_mod7_subcolor
  -> ph_bad_phase_derivative
  -> outer_extension_period_shift
  -> mod35_boundary_color_word
  -> transition_increment_word
  -> raw_escape_count
  -> raw_color_analogy
```

## Assumption Challenge

Alternate vertices considered:

```text
runners, residues mod 35, residues mod 7, high gaps, color classes,
bad-coloring nodes, endpoint addresses, fixed corridors, proof obligations.
```

The chosen quotient preserves AP-tail escape count, endpoint-witness location,
mixed/pure phase, and outer-extension shift.  It destroys arbitrary non-AP
wall geometry, so it is a coloring-recursion sidecar for the AP-tail bridge,
not a global LRC14 proof.
