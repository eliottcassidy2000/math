---
id: HYP-3148
title: The n=4 tournament tables expose an Erdos-870 live-core/filler/canary normal form
status: SYNTHESIS / exact n=4 table scout; not an LRC14 proof
source: codex-2026-06-27-S275
tangent: T1213
technique: LTI-274
tournament_technique: LTT-172
script: 04-computation/lrc14_erdos870_live_core_canary_scout_codex_s275.py
result: 05-knowledge/results/lrc14_erdos870_live_core_canary_scout_codex_s275.out
reflection: 07-reflections/erdos870-live-core-filler-canary-codex-s275.md
related:
  - HYP-3147
  - HYP-3146
  - HYP-3145
  - HYP-3144
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3137
  - HYP-3134
  - HYP-3133
  - HYP-3124
  - HYP-3054
  - HYP-2534
  - OPEN-Q-108
external:
  - https://github.com/davidturturean/erdos-870
---

# HYP-3148: Erdos-870 Live-Core/Filler/Canary Deletability Audit

## Claim

The two n=4 tables are not competing encodings of the same small tournament
fact.  They are two proof modes:

```text
three-bit tiling cube = fixed Hamiltonian path + live skips a,b,c
two-bit anchor       = fixed filler/canary arcs + live opposite skips x,y
```

The second is the first with the long diagonal `c` frozen.  S276/HYP-3143
already turns this into an exact-order subbasis audit: a theorem-facing packet
should be represented at its declared order and not below it.  Incoming
S274/HYP-3144 adds the adjacent pair-function guardrail: unordered quotients
preserve `a+b`/`a*b`-type data but ordered functions need sidecars.  Incoming
HYP-3145 adds the filler-core interface: fixed-path tilings are atlases and
alarms, while terminal rows should expose the small retained core.  Incoming
S274/HYP-3146 upgrades this into a shift-package policy: keep canary fiber mass
when deletion stability matters, or add finite filler/scaffold data when
gluing/classification needs a congruent quotient.  Incoming S277/HYP-3147 adds
the local n=3 edge-flip kernel and ordered-function sidecars.  This HYP-3148
continuation adds the complementary deletable-coordinate audit: a proof-facing
quotient must name its live core, deterministic fillers, canary/shift controls,
deletable coordinates, and terminal exit before class counts are allowed to
stand for proof support.

For the LRC14 packet, the proposed new fields are:

```text
live_core_bits
filler_bits
canary_bits
deletable_coordinates
class_distribution
minimal_cover_subbasis
edge_bounded_core_floor_exit
terminal_exit_or_named_debt
```

## Evidence

The scout `04-computation/lrc14_erdos870_live_core_canary_scout_codex_s275.py`
enumerates all labelled n=4 tournaments and classifies them by the four score
sequence families:

```text
T = (0,1,2,3)
S = (1,1,2,2)
+ = (0,2,2,2)
- = (1,1,1,3)
```

With a fixed Hamiltonian path `0->1->2->3`, the live skip edges are
`a=(0,2)`, `b=(1,3)`, and `c=(0,3)`.  The colored xor table on
`{E,a,b,c}` is exactly the user's table:

```text
      E  a  b  c
  E | T  +  -  S
  a | +  T  S  S
  b | -  S  T  S
  c | S  S  S  T
```

But the full three-bit cube is skewed:

```text
class_counts_full_cube = {T:1, +:1, -:1, S:5}
entropy_bits = 1.548795
minimal_all_class_covers = [('a','b')]
```

So `c` is useful as an immediate `S` witness, but it is not load-bearing for
class coverage: `{a,b}` already reaches all four classes.

The two-bit model fixes four arcs with partial score sequence `(0,1,1,2)`.
The canonical anchor is:

```text
fixed = (0,1), (0,3), (1,2), (2,3)
live x = (0,2)
live y = (1,3)
```

It gives the user's second table:

```text
      E  x  y
  E | T  +  -
  x | +  T  S
  y | -  S  T
```

and the full two-bit square is uniform:

```text
class_counts_full_square = {T:1, +:1, -:1, S:1}
entropy_bits = 2
minimal_all_class_covers = [('x','y')]
```

The scout finds `24` labelled two-bit anchors, all with disjoint live-edge
pairs; each of the three perfect matchings of `K4` contributes `8`.

## Erdos-870 Transfer

The Erdos-870 repository proves a negative answer to the question whether a
logarithmic lower bound on representation counts forces a minimal additive
subbasis.  Its public statement is an order-`k` basis with at least
`C log n` representations for large `n`, but no minimal order-`k` subbasis.
Its proof architecture separates an order-two probabilistic source,
deterministic filler reductions for `k>=4`, and a clustered canary construction
for `k=3`.

The n=4 tournament table is the small finite analogue:

- the live `x,y` square is the order-two core;
- the fixed four arcs are the deterministic filler;
- freezing `c` is the canary/shift control that removes representation skew;
- the three-bit tiling cube shows why witness abundance alone is not
  minimality.

This reframes LRC tournament evidence.  A row with many witnesses, many
tournament encodings, or many score-class representations may still have a
deletable coordinate.  Conversely, a sparse two-bit quotient can be the real
proof carrier if its filler and canary coordinates are named.

## LRC14 Use

This lane attaches to:

- HYP-3141: an edge is proof-facing only after tail, tip, orbit, commutator,
  coordinate-repair, and terminal-exit fields are present.
- HYP-3140: `Rprime` must be carried as a finite fiber-PGF conditional first
  moment before scalarizing.
- HYP-3142: the k=8 hard node needs an exact `U4` terminal sidecar, not only a
  raw score/coverage value.
- HYP-3133/HYP-3134: A000568 sits between local sector counts and full paired
  child packets, so quotients need global-consistency certificates.

The new proof target is not a theorem about n=4 tournaments.  It is a ledger
rule: before accepting any tournament quotient inside the LRC14 proof, prove
that the quotient has a nondeletable live core or explicitly record the
deletable/filler/canary debt.

## Guardrail

This is an exact n=4 finite scout and a proof-interface proposal.  It does not
prove LRC14, and it does not import any additive-basis theorem as a black box.
The value is methodological: representation counts, class counts, and score
families are proof shadows until the live coordinates and deletable coordinates
are audited.
