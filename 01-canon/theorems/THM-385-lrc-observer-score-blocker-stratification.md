---
id: THM-385
name: lrc-observer-score-blocker-stratification
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S517
depends_on:
  - THM-381
  - THM-384
---

# THM-385: Observer score is exactly the LRC blocker count

## Statement

Let `T_S(t)` be the observer-marked LRC tournament from THM-381 on a speed
system

```text
S = {0, v_1, ..., v_{n-1}}
```

with threshold `1/n`.  The observer is the marked vertex `0`, and observer
incident edges are oriented by

```text
0 -> i    iff    ||v_i t|| >= 1/n.
```

Define the blocker count at time `t` by

```text
B(t) = #{ i in {1,...,n-1} : ||v_i t|| < 1/n }.
```

Then

```text
outdeg_T(0) = (n - 1) - B(t)
indeg_T(0)  = B(t).
```

Consequently:

1. The observer is a source iff `B(t)=0`, recovering THM-381.
2. The observer is an almost-source with exactly one incoming edge iff exactly
   one moving runner blocks the observer.
3. More generally, the observer-score layer is exactly the distance-to-source
   layer measured in observer incident edge flips.
4. The blocker count is invariant under observer-marked tournament
   isomorphism, because marked isomorphism preserves the marked vertex score.

## Proof

For each moving runner `i`, the observer incident edge is defined independently
by the LRC threshold:

```text
0 -> i  iff  ||v_i t|| >= 1/n.
```

Thus runner `i` contributes `1` to `outdeg_T(0)` exactly when it is safe from
the observer, and contributes `1` to `indeg_T(0)` exactly when it is a blocker,
i.e. when `||v_i t|| < 1/n`.

The `n-1` moving runners partition into safe runners and blockers, so the two
displayed degree identities follow immediately.  The source and almost-source
statements are the cases `B(t)=0` and `B(t)=1`.

Changing a blocker edge `i -> 0` into `0 -> i` is exactly one observer incident
edge flip.  Therefore `B(t)=indeg_T(0)` is the number of observer incident edge
flips required to make the observer a source while holding all runner-runner
edges fixed.

Finally, an observer-marked tournament isomorphism fixes the marked vertex and
preserves adjacency.  It therefore preserves both `outdeg_T(0)` and
`indeg_T(0)`, hence preserves `B(t)`.

## Verification Record

The S517 audit script checks the theorem over exact open and wall samples:

```text
04-computation/lrc_observer_predicate_zoo_s517.py
05-knowledge/results/lrc_observer_predicate_zoo_s517.out
```

The stored audit reports zero score mismatches in bounded windows:

```text
N=4 max_speed=12 states=21040 score_mismatches=0
N=5 max_speed=10 states=25816 score_mismatches=0
N=6 max_speed=9  states=18952 score_mismatches=0
N=7 max_speed=8  states=4064  score_mismatches=0
```

It also reports zero marked-class mixing by blocker count in all four windows.
Runner-subtournament classes do mix blocker counts, confirming that deleting
the observer incident threshold edges loses this LRC stratification.

## Significance

THM-381 identifies the source layer as the LRC witness.  THM-385 identifies all
nearby observer-score layers: they are not arbitrary tournament scores but the
exact finite ledger of how many runners currently block the observer.

This gives a precise meaning to source-like predicates:

```text
source        = zero blockers
almost-source = one blocker
score layer k = k blockers
```

It also explains why HYP-1981 must remain observer-marked.  The runner-runner
half-turn sub-tournament can be useful geometry, but by itself it does not
remember the LRC distance-to-source layer.

## Related

- THM-381: observer-source marked reachability.
- THM-384: observer-adjacent gap source criterion.
- HYP-1981: source-reachability in the observer-marked A000568 quotient.
- HYP-1988: observer-score repair layers.
