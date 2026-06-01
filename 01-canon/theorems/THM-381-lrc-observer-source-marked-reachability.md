---
id: THM-381
name: lrc-observer-source-marked-reachability
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S513
depends_on:
  - THM-372
---

# THM-381: LRC loneliness is observer-source reachability in a marked tournament

## Statement

Let

```text
S = {0, v_1, ..., v_{n-1}}
```

be a Lonely Runner speed system with stationary observer `0` and threshold
`1/n`.  For a time `t`, define the observer-source LRC tournament `T_S(t)` on
vertices `{0,1,...,n-1}` as follows.

For observer-runner pairs, orient

```text
0 -> i    iff    ||v_i t|| >= 1/n,
i -> 0    otherwise.
```

For runner-runner pairs, use any complete tie-resolved Tournament Analysis
comparator, for example the half-turn circular comparator with a fixed
Hamiltonian tie path.

Then the observer is lonely at time `t` iff the marked vertex `0` is a source
of `T_S(t)`:

```text
forall i, ||v_i t|| >= 1/n
  <=>  outdegree_T(0) = n-1.
```

Consequently, for a fixed speed system, the LRC witness problem is exactly the
reachability problem asking whether the marked tournament movie visits a class
whose marked vertex is a source.  A closed-threshold counterexample is exactly
a marked tournament movie that avoids all observer-source classes.

Moreover, the number of observer-source marked tournament isomorphism classes
on `n` vertices is

```text
A000568(n-1),
```

the number of unmarked tournament isomorphism classes on the `n-1` moving
runners.

## Proof

The first equivalence is a direct unpacking of the observer-runner arcs.  By
definition, the observer is lonely at `t` exactly when every moving runner is
at circular distance at least `1/n` from the observer:

```text
forall i in {1,...,n-1}, ||v_i t|| >= 1/n.
```

The observer-source tournament orients `0 -> i` exactly under the same
condition.  Thus the displayed condition holds for all `i` iff `0` beats every
other vertex, i.e. iff `0` is a source.  The orientations among the moving
runners play no role in this source test; they only place the witness state
inside a complete marked tournament class.

Source-ness of the marked vertex is invariant under marked tournament
isomorphism.  Therefore the existence of a lonely time is equivalent to the
marked tournament trajectory visiting the source-marked subset of the marked
class space.

For the count, let `M_n^src` be the set of observer-marked isomorphism classes
of tournaments on `n` vertices in which the marked vertex is a source.  Deleting
the marked source leaves an arbitrary tournament on the remaining `n-1`
vertices.  This gives a map

```text
M_n^src -> {unmarked tournament classes on n-1 vertices}.
```

Conversely, given any tournament `U` on `n-1` vertices, add a new marked vertex
`0` and orient `0 -> u` for every vertex `u` of `U`.  This constructs an
observer-source marked tournament.  The two operations are inverse on
isomorphism classes: any marked isomorphism must fix the marked source and
therefore restricts to an isomorphism of the runner sub-tournaments, while any
runner sub-tournament isomorphism extends uniquely by fixing the marked source.

Thus `|M_n^src|` is exactly the number of unmarked tournament isomorphism
classes on `n-1` vertices, namely `A000568(n-1)`.

## Verification Record

The construction was independently audited in

```text
04-computation/lrc_observer_source_tournament_s511.py
05-knowledge/results/lrc_observer_source_tournament_s511.out
```

The stored audit reports zero mismatches between observer-source status and
direct LRC safety on the tested families:

```text
n=4 initial (0,1,2,3):        mism 0
n=4 (0,1,2,5):                mism 0
n=5 initial (0,1,2,3,4):      mism 0
n=5 (0,2,3,5,7):              mism 0
n=6 initial (0,1,2,3,4,5):    mism 0
n=6 (0,1,3,4,5,9):            mism 0
```

The source-target count was checked in

```text
05-knowledge/results/lrc_qr_paley_speed_tournament_s511.out
```

for

```text
n=4..10: A000568(n-1) = 2,4,12,56,456,6880,191536.
```

The computation is not needed for the proof, but it records that the
implemented class-count convention agrees with the theorem.

## Significance

This theorem formalizes the exact part of HYP-1981.  The earlier A000568
projection threads found that raw half-turn tournament class does not determine
LRC safety.  The fix is to change the observer incident arcs to the LRC
threshold arcs and keep the observer marked.  In that marked class space,
loneliness is exactly source reachability, and the target set has the clean
size `A000568(n-1)`.

The remaining open problem is not this equivalence.  It is to characterize the
reachable source-marked classes for arithmetic runner clocks and prove that
every primitive speed system reaches one of them.

## Related

- HYP-1977: LRC projection defect over A000568.
- HYP-1978: marked A000568-style quotient.
- HYP-1979: sparse marked chamber walk.
- HYP-1981: source-reachability in the observer-marked A000568 quotient.
- `07-reflections/lrc-as-source-reachability-in-marked-A000568-s511.md`
