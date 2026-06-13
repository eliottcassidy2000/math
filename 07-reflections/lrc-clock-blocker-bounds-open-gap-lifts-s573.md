---
source: codex-2026-06-03-S573
status: exact bounded audit + correction to S572
tags:
  - lonely-runner
  - clock-blockers
  - second-gap
  - open-gap
  - antipodal-transversal
  - unit-action
  - divisibility
  - n-clock
---

# Clock-Blocker Bounds and the First Open-Gap Lifts

The useful sharpening is not another slogan about the AP.  It is a list of
clocks that any strict sub-edge row must already defeat.

For `k=n-1` moving runners, write the proposed edge as

```text
E_n = 2/(2n-1).
```

S572 saw only floor-tight rows below `E_n` in the original small primitive
boxes.  S573 expands those boxes and finds that this floor-only reading is not
global.  The better result is a clock-blocker ledger.

## Three Necessary Gates

If `M(S) < E_n`, then `S` must pass all three gates below.

### D: small denominator blockers

For every

```text
2 <= q <= n-1,
```

some speed must be divisible by `q`.  Otherwise the clock

```text
t = 1/q
```

makes every runner land at a nonzero `q`-th root, so the score is at least

```text
1/q >= 2/(2n-1).
```

This is a clean bound because it is independent of the `2n-1` summand graph.
It says any strict sub-edge row must already look like a divisibility cover of
the small clocks.

### U: unit-shell coverage

At the odd modulus

```text
C = 2n-1,
```

every unit antipodal shell `{a,C-a}` must be hit.  If a unit shell is missed,
S553 chooses `k=a^{-1}` and the time `k/C` witnesses the edge.

This is the summand-graph/unit-action bound from S571.

### N: n-clock blockers

For every tick `j/n`, at least one runner must land at distance `0` or `1/n`.
If all runners are at distance at least `2/n`, then that `j/n` clock scores

```text
2/n > 2/(2n-1).
```

So the floor clock is not only an equality witness.  It is also a blocker
ledger: strict sub-edge rows have to spoil every `n`-clock tick.

## Exact Expansion

`04-computation/lrc_second_gap_bounds_s573.py` applies those gates first and
then exact-checks only the survivors.  The expanded exact boxes are:

```text
k=3, B=60
k=4, B=40
k=5, B=32
k=6, B=26
k=7, B=22
k=8, B=20
k=9, B=18
k=10, B=17
```

The survivor counts are:

```text
2623, 5676, 3237, 1446, 4654, 238, 145, 111.
```

That is the practical win: a huge primitive box becomes an exact-checkable
clock-obligation ledger.

## The Correction: Open-Gap Rows Exist

The expanded audit finds real rows in the interval

```text
1/n < M(S) < 2/(2n-1).
```

The first one appears at `n=7`:

```text
S = (1,5,6,11,16,17)
M(S) = 5/33
t = 10/33
```

This is a perfect antipodal transversal modulo `13`, with flip-set `{2}`.
So the `{2}` flip does not simply die after `n=6`; lifted representatives let
it rise above the floor while staying below the old edge.

At `n=8`, S573 finds:

```text
(1,2,3,4,5,7,18)     M=3/23
(1,3,4,5,7,13,18)    M=3/23
```

These are not perfect transversals modulo `15`; they are nonunit-hole rows
missing `{6,9}` and doubling `{3,12}`.  Again, they are safely above the LRC
floor `1/8`, but below `2/15`.

So the corrected picture is:

```text
2/(2n-1) is not a global second value over lifted integer rows.
```

It remains a crucial witness edge for unit-shell and lower-residue reductions,
but the integer lift direction creates intermediate values.

## What Bound Actually Got Sharper

The old proof target was too blunt:

```text
below 2/(2n-1) => floor-tight.
```

The new target is sharper and more plausible:

```text
below 2/(2n-1)
=> D_q blockers + U_a blockers + N_j blockers
=> lifted floor/open-gap families
=> M(S) >= 1/n.
```

This preserves LRC while giving up the false spectral-gap statement.

## Addition, Multiplication, Odd, Even

The user's four-way split gets cleaner:

```text
addition       -> antipodal shell ledger U_a
multiplication -> unit-inverse visibility for missed shells
divisibility   -> small denominator blockers D_q
even/floor     -> n-clock blockers N_j
```

Odd `2n-1` still removes the midpoint and makes shells honest.  But the integer
representatives live above the residue quotient, and multiplication by lifts can
create new pair-sum witnesses like `10/33` and `4/23` between the floor and the
unit-shell edge.

## Tournament Analysis / Assumption Challenge

The vertices should be proof obligations:

```text
D_q, U_a, N_j
```

rather than runners.  A runner-vertex tournament would destroy exactly the
predicate that matters: which clocks are blocked.

One useful tournament gauge is:

- **Vertices:** obligations `D_q`, `U_a`, `N_j`.
- **Observable:** number of speeds blocking the obligation, then smallest
  blocking speed, then whether the blocker is also active at the exact witness.
- **Switch:** orient toward the more fragile obligation.
- **Tie Hamiltonian path:** denominator/order first (`D`, then `U`, then `N`),
  with increasing modulus label.

This preserves the proof predicate: a counterexample or near-counterexample has
to be a simultaneous hitting set for all obligations.

## Short Version

S573 turns "look below the second edge" into a finite blocker ledger.  It also
corrects the repo: there are genuine open-gap lifted rows below `2/(2n-1)`.
The next theorem should classify those lifts, not deny their existence.

Incoming Burnside/Fourier work (HYP-2085 and HYP-2087) points at the same next
object from the dual side: enrich binary lonely time words with labelled clock
obligations.  The `D_q/U_a/N_j` ledger is exactly that missing label layer.
