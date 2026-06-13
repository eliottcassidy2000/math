---
id: HYP-2077
status: PROGRESS - exact cover gaps found in frontier rows; theorem route open
source: codex-2026-06-02-S564
related:
  - HYP-2062
  - HYP-2065
  - HYP-2075
  - HYP-2073
  - HYP-2076
  - THM-369
  - THM-396
---

# HYP-2077: n=14 can be attacked as an exact even-good cover problem

## Statement

For `n=14`, split a primitive 13-speed set into even and odd speeds and define

```text
G = {t in [0,1): ||v t|| >= 1/14 for every even speed v}.
```

The even-fold reduction says `G` is the free region supplied by lower-runner
LRC.  The remaining problem is an exact cover problem:

```text
do the odd danger arcs D_v = {t: ||v t|| < 1/14} cover G?
```

HYP-2077 proposes the n=14 proof trichotomy:

```text
coarse denominator witness
or positive exact gap in G minus odd danger
or wall endpoint witness
or owner-labelled residual frontier export.
```

Equivalently, a genuine counterexample must do more than cover `G` in measure.
It must cover every positive component and also kill all boundary endpoint
witnesses.  The known wall rows cover `G` in measure but retain closed
rational witnesses.

## Evidence

`lrc_n14_evenfold_cover_angles_s564.py` computes the interval sets exactly.

Wall rows:

```text
AP_wall:
  |G| = 16/35
  |G minus D_odd| = 0
  first closed witness = 1/14
  every positive G component has a complete odd-cover chain

V_star_wall:
  |G| = 37/84
  |G minus D_odd| = 0
  first closed witness = 1/14
  every positive G component has a complete odd-cover chain
```

Thus the wall is not "odd arcs fail to cover in measure."  The wall is
"odd arcs cover in measure, but endpoint witnesses survive."

Positive-gap n=14 rows:

```text
sieve_minimal_lonely:
  |G minus D_odd| = 3/49
  first closed witness = 1/16

near_AP_apex:
  |G minus D_odd| = 426/35035
  first closed witness = 1/12

S562_packet_n14:
  |G minus D_odd| = 1/143
  first closed witness = 6/23

S562_packet_n14_lift:
  |G minus D_odd| = 1/143
  first closed witness = 3/23

no_small_pinch_proxy:
  |G minus D_odd| = 1405685802395447/11291147314967520
  first closed witness = 2/15
```

The important samples are the S562 packets and the no-small-pinch proxy.  They
are not solved by the naive small-pinch route alone, and the S562 packet is
sieve-covered for `q <= 14`; nevertheless the exact even-fold cover has
positive uncovered mass and the first closed witness lies in the fine window
or nearby denominator tier.

Incoming HYP-2075 explains why this is the right modulus scale: the complete
pinch primitive uses pair-sum moduli rather than a single apex-sensitive
denominator.  In S564 language, the first uncovered gaps in `G` should be
converted into owner-labelled pair-sum residues, not treated as arbitrary
large-denominator accidents.

## Greedy Cover Chains

For each positive component of `G`, S564 runs the canonical interval-cover
greedy algorithm using odd danger intervals.  A full cover produces a
left-to-right handoff chain.  A gap produces an explicit lonely interval.

First gaps:

```text
sieve_minimal_lonely: 1/28 .. 13/196, length 3/98
near_AP_apex:         15/182 .. 13/154, length 2/1001
S562_packet_n14:      127/1078 .. 139/1176, length 5/12936
S562_packet_n14_lift: 211/2548 .. 65/784, length 1/10192
no_small_pinch_proxy: 197/2660 .. 27/364, length 1/8645
```

This creates a new proof target.  If every component is covered, the odd arcs
must form no-gap handoff chains.  The user's hidden transitivity/no-return fact
should be applied to those handoff events:

```text
if X hands to Y and Y hands to Z, then X must precede Z and a return triangle
Z -> X -> Y -> Z is forbidden.
```

S564 has not yet retained enough ownership data to prove this.  It identifies
the next quotient: odd runner arcs should be refined into owner-labelled
danger intervals and endpoint walls.

## Tournament Analysis

Vertices were odd runners inside a fixed speed-set obligation.

Pairwise observable:

```text
exclusive odd-danger pressure on G:
measure((D_x cap G) minus D_y) versus measure((D_y cap G) minus D_x).
```

Switch/gauge:

```text
x -> y when x has larger exclusive pressure; ties follow increasing speed as
the Hamiltonian path.
```

Fingerprints:

```text
All audited rows have directed_3_cycles = 0, singleton SCCs, and one
Hamiltonian path under this quotient.
```

This transitivity is useful but also too coarse.  It preserves a pressure order
among odd runners, but it destroys the event-level handoff ownership needed to
use the no-return obstruction.  The next Tournament Analysis should use
vertices:

```text
(odd runner, danger interval, G component, endpoint owner).
```

## Assumption Challenge

Possible vertex sets considered:

```text
runners, odd runners only, G components, danger intervals, greedy handoff
events, endpoint walls, residues, CRT channels, cover arcs, fine denominators,
and proof obligations.
```

Chosen quotient:

```text
vertices = odd runners in the even-fold quotient.
```

Predicate preserved:

```text
relative odd pressure on the exact even-good set G.
```

Information destroyed:

```text
which interval endpoint is responsible for a handoff, boundary equality
witnesses, residue ownership, and local no-return triangles among interval
events.
```

Challenged assumption:

```text
the even-fold reduction is only a measure argument.
```

S564 shows it can be made exact.  The distinction between positive safe mass
and zero-measure endpoint witnesses is precisely the wall/non-wall split.

## Open Work

1. Prove an exact cover-gap lemma: outside AP/V*-type walls, some component of
   `G` has an uncovered subinterval.
2. Classify the zero-measure wall covers and prove their endpoint witnesses
   survive.
3. Refine odd runners into owner-labelled danger intervals and apply the
   no-return/cycle-exclusion fact to full cover chains.
4. Merge with HYP-2076: if the exact cover gap is found, extract its low
   denominator witness; otherwise export owner-labelled frontier mass.

## Files

- `04-computation/lrc_n14_evenfold_cover_angles_s564.py`
- `05-knowledge/results/lrc_n14_evenfold_cover_angles_s564.out`
- `07-reflections/lrc-n14-evenfold-cover-gap-route-s564.md`
