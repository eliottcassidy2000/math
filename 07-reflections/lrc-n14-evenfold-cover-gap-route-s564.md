---
source: codex-2026-06-02-S564
status: reflection + exact interval-cover audit
tags: [LRC, n14, even-fold, interval-cover, recursive-sieve, Tournament-Analysis]
---

# n=14 even-fold cover gap route

The S563 lower-bound ledger said: after the coarse sieve goes blind, try a
recursive local witness, or export positive frontier mass.  For `n=14`, the
cleanest place to test that idea is the even-fold reduction.

Define

```text
G = {t: every even speed is at distance >= 1/14}.
```

Lower-runner LRC gives this region for free in principle.  The n=14 question is
then:

```text
do the odd danger arcs cover G?
```

S564 computes this exactly as an interval problem.

## What changed

The old even-fold script measured `G` on a grid.  S564 keeps exact rational
components.  That exposes the right wall distinction.

For the AP and V* wall rows:

```text
|G minus odd danger| = 0
first closed witness = 1/14
```

So the wall is not a positive-measure gap.  It is an endpoint phenomenon.  The
odd arcs cover `G` in measure, but boundary points survive as closed lonely
witnesses.

For non-wall rows, the proof looks much friendlier:

```text
sieve_minimal_lonely: |safe| = 3/49,       witness 1/16
near_AP_apex:         |safe| = 426/35035,  witness 1/12
S562 n14 packet:      |safe| = 1/143,      witness 6/23
S562 n14 lift:        |safe| = 1/143,      witness 3/23
no-small-pinch proxy: |safe| ~= 0.124495,  witness 2/15
```

This is a real n=14 angle: prove every non-wall residual has a positive exact
gap in `G`, and separately classify the zero-measure wall covers.

The incoming HYP-2075 pair-sum-modulus result plugs directly into this:
uncovered intervals in `G` should be searched and certified through pair-sum
pinch residues, not only through small integer denominators or the single
`14`-apex lens.

## Greedy chains and no-return

On each component of `G`, odd danger intervals either cover left-to-right or
leave a first gap.  S564 records the first gaps exactly.  For the S562 packet,
one gap is:

```text
127/1078 .. 139/1176, length 5/12936.
```

The user's hidden transitivity fact belongs here.  A full odd cover of a
component is a handoff chain:

```text
X covers up to a point, Y extends it, Z extends it again.
```

If `X -> Y` and `Y -> Z`, then `X -> Z`; and the return triangle
`Z -> X -> Y -> Z` must be impossible.  At runner level the pressure tournament
is transitive in the audited rows, but that quotient is too coarse.  The next
version should make the vertices the actual handoff events:

```text
(odd runner, danger interval, G component, endpoint owner).
```

Then the no-return fact can attack the only remaining way a counterexample can
exist: every component has a complete handoff chain and every endpoint witness
is killed.

## The proof trichotomy

This suggests a sharper n=14 program:

```text
1. Coarse sieve: some q <= 14 divides no speed.
2. Even-fold exact gap: odd danger leaves positive mass in G.
3. Wall endpoint: odd danger covers G in measure, but closed endpoints survive.
4. Residual export: if a full cover chain remains, export owner-labelled
   frontier mass to the recursive sieve state.
5. Pair-sum witness extraction: when a gap appears, translate it to the
   HYP-2075 pinch-modulus primitive whenever possible.
```

This is more concrete than "try another lower bound."  It is a finite exact
object with rational endpoints, a visible wall/non-wall split, and a place for
the transitivity/no-return obstruction to act.

## Assumption challenge

I considered using runners, all speed sets, denominator tiers, G components,
odd danger intervals, endpoint walls, residues, CRT channels, cover arcs, and
handoff events as vertices.

S564 uses odd runners because that preserves the first predicate:

```text
which odd speeds uniquely pressure G?
```

But it destroys event ownership.  The next proof attempt should refine the
vertex set; the pressure order is only the shadow.

## Handoff

Build `S565` as an owner-labelled cover-chain verifier:

```text
for each G component:
  compute the greedy odd cover chain;
  label each extension endpoint by runner/residue;
  build the handoff tournament on extension events;
  count return triangles and endpoint wall kills.
```

If the only complete chains with no return defects are AP/V*-type wall chains,
that would be a serious n=14 proof route.

**Artifacts:** `04-computation/lrc_n14_evenfold_cover_angles_s564.py`,
`05-knowledge/results/lrc_n14_evenfold_cover_angles_s564.out`, HYP-2077.
