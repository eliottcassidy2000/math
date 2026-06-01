---
source: oracle-2026-06-01-S24
status: exploratory synthesis (Tournament Analysis: the runner clock)
tags:
  - tournament-analysis
  - lonely-runner
  - metagraph
  - tournament-clock
  - circular-tournaments
  - wall-crossings
---

# The Tournament Clock: Runner Systems as Closed Walks in G_n

Continuing Tournament Analysis (S22/S23). S23 built the stack
`positions → metrics → comparators → tournaments → invariants → wall-crossings`
and deferred the pattern-hunt. This session runs the hunt for one canonical
comparator and finds that a runner system is, very cleanly, **a closed walk
through a fixed tiny corner of the tournament metagraph `G_n`** — and that the
walk's geometry reads off the speed structure.

## The phase comparator and the clock

Take `n` runners with integer speeds `s_0 < … < s_{n-1}` on the unit circle,
positions `x_i(t) = frac(s_i t)`. Lift to a tournament by the **phase
comparator**

```
i → j   iff   frac(x_i − x_j) = frac((s_i − s_j) t) ∈ (0, 1/2)
```

— "i leads j by less than half a lap." Ties (`frac ∈ {0,1/2}`) break by the base
path `0→1→…→n−1` (the speed-rank labels — the basketball "1..5" of S23). This is
the *unique* threshold giving a tournament: `α = 1/2` makes exactly one of
`i→j, j→i` hold. As `t` sweeps `[0,1)`, `T(t)` is piecewise constant and traces a
**closed walk**. Each edge `(i,j)` flips exactly at `t = m/(2|s_i−s_j|)`, so the
cells are exact rationals (`tournament_clock_s24.py`).

## Finding 1 — the menu of tournaments is FIXED per n (speeds only pick the walk)

Across every speed family tried, the set of `H`-values (directed Hamiltonian
path counts) of the visited tournaments is the *same finite set*, depending only
on `n`, not on the speeds:

```
n=5 : H ∈ {1, 9, 11, 15}     (4 iso-classes, of 12 tournaments on 5)
n=6 : H ∈ {1, 17, 23, 41, 45}
n=7 : H ∈ {1, 33, 47, 51, 105, 123, 137, 151, 175}
```

The phase comparator can only ever produce the **"circular" tournaments** — those
realizable by `n` points on a circle under the half-turn rule — a *fixed small
sub-family* of `G_n`. The speeds do not change *which* tournaments are possible;
they change only *which are visited, how often, and in what cyclic order*. This
is the right way to read S23's "a metric is not automatically a new lens": the
phase lens has a fixed image in `G_n`, and Tournament Analysis of runners is the
study of **walks inside that fixed image.**

**What the n=5 menu is (exactly).** The 4 realizable classes form a *chain
ordered by spread*:

```
H=1   score (0,1,2,3,4)   transitive        — runners bunched in a semicircle
H=9   score (1,1,2,3,3)   near-transitive
H=11  score (1,2,2,2,3)   intermediate
H=15  score (2,2,2,2,2)   REGULAR, ROTATIONAL (R_5) — runners evenly spread
```

The top of the chain is the unique regular tournament on 5 vertices, the
rotational/quadratic-residue `R_5` (max `H = 15`). It is realised exactly when
the runners are near the **regular pentagon** — and the cell containing `R_5`
contains the lonely-runner witness times `t = a/5` (positions
`{0,a/5,…,4a/5}`, all at distance `≥ 1/5` from the observer). So the *most-spread
tournament cell hosts the loneliness times*; bunching (transitive, `H=1`) and
even spread (regular, `H=15`) are the two ends, and the clock oscillates between
them.

n=6 confirms the fixed menu: 6 classes (of 56), `H ∈ {1,17,23,23,41,45}`, again a
spread chain from transitive `(0,1,2,3,4,5)` up to the **near-regular** top
`(2,2,2,3,3,3)` — and here two *distinct* classes share that top score with
`H = 41` vs `45` (no regular tournament exists on even `6`, so the top splits).

**Geometry constrains the image (basketball contrast).** The circular menu is a
property of *circle geometry*, not of tournaments. Feeding *arbitrary* pairwise
data through the discrete flux comparator (the basketball pass-matrix lift)
reaches tournaments **outside** the circular menu: two sample teams gave `H = 5`
and `H = 13`, neither in `{1,9,11,15}`. So:

> the comparator's geometry fixes its reachable set in `G_n`. Points-on-a-circle
> can only build 4 of the 12 tournaments on 5 vertices; arbitrary relations build
> all 12.

This is the precise content of "a metric is not automatically a new lens" (S23):
a *geometric* comparator has a constrained, structured image; only unconstrained
data fills `G_n`.

## Finding 2 — transitive ⟺ empty semicircle ⟺ a ½-gap (the LRC bridge)

The transitive tournament (`H = 1`, score `0,1,…,n−1`) is visited iff **all
runners lie inside some open semicircle** — i.e. there is an empty half-circle, a
*gap of length ≥ 1/2*. Verified exactly: across the tested systems, "cell is
transitive" matched "max circular gap > 1/2" with **zero mismatches**.

So the comparator's *transitive moments are exactly the loneliness moments at
threshold `1/2`* (the `n=2` lonely-runner gap). The half-turn comparator is the
`α = 1/2` member of a family; finer LRC thresholds `1/n` correspond to richer
comparators (the two-neighbour bracket of S22). **The comparator threshold *is*
the lonely-runner gap threshold** — one knob unifying the two programs.

S26b correction: the next paragraph is reliable for the `n=5` chain and for the
endpoint `H=1`, but not as a global scalar theorem.  Later samples through `n=9`
show that, once `H>1`, max-gap ranges overlap and pointwise monotonicity fails.
The corrected statement is: `H=1` exactly detects the open-semicircle state, and
larger `H` values are circular-tournament class features correlated with spread,
not a scalar max-gap meter.

For the `n=5` chain, bucketing clock cells by `H` and recording the runners'
max circular gap (`arith 0..4`) gives:

```
H=1   max_gap = 0.750   (transitive — bunched, an empty semicircle)
H=9   max_gap = 0.417
H=11  max_gap = 0.292
H=15  max_gap ≈ 0.26    (regular — most even spread)
```

At all `n`, `H=1` still means the clean pecking-order/empty-semicircle state.
For `H>1`, the Hamiltonian-path count reads a coarser circular arrangement
class.  Few Hamiltonian paths still indicate bunched/transitive-side behavior,
and many paths usually indicate even spread, but the statement is correlation
and class structure rather than pointwise max-gap monotonicity.

## Finding 3 — the LRC-extremal set is the MINIMAL clock

The arithmetic / lonely-runner extremal speeds `(0,1,…,n−1)` give the *smallest*
clock: fewest cells, fewest iso-classes (4 at both n=5 and n=6), and exactly
**2 transitive cells**. Spread/coprime speeds give many cells and more
iso-classes (e.g. n=6 coprime `0,1,5,7,11,13`: 104 cells, 6 classes, 14
transitive cells). So:

> the extremal lonely-runner configuration = the minimal, most resonant
> tournament-clock orbit.

This is the dot the user asked to connect: the speed set that is hardest for LRC
(the tight `1,…,n−1`) is the one whose tournament walk is *simplest* — maximal
resonance = minimal walk. The number of cells scales with the total spread
`Σ|s_i − s_j|` (the clock "ticks" once per `m/(2|s_i−s_j|)`), and arithmetic
speeds minimise the number of *distinct* wall times by sharing denominators —
resonance literally collapses wall-crossings.

## Finding 4 — the walk is a loop in G_n; its signature is the new invariant

Per runner system we now have a closed walk in `G_n` with a computable
**signature**: `(#cells, #distinct iso-classes, H-multiset with multiplicities,
#transitive cells, the cyclic class-sequence)`. Two speed sets with the same
signature are "Tournament-Analysis equivalent." Conjecturally the signature is a
`GL`-type invariant of the speed set's resonance lattice (the integer relations
among `s_i − s_j`), since wall-times and their coincidences are governed exactly
by those differences. This makes Tournament Analysis a bridge: **runner
resonances ↔ closed-walk homology in the metagraph.**

## The basketball anchor and the "corresponding tournament"

The discrete model: 5 starters, `i→j` if `i` passes to `j` more than vice-versa,
ties broken by the jersey path `1→2→3→4→5`. That tournament's `H`, score
sequence, and iso-class place the team as a *point* in `G_5`; a season is a *walk*
in `G_5` (one point per game), exactly the discrete analogue of the runner clock.
The user's "corresponding tournament" is then naturally: the runner system's walk
has a most-visited class / a canonical-time class — the team's "identity
tournament" — and Tournament Analysis compares systems by their identity class
and walk signature, not by raw metrics.

## Patterns / next

1. Prove Finding 1: the phase-comparator image in `G_n` is exactly the circular
   (point-on-circle half-turn) tournaments; compute that family's size for `n≤8`.
2. Make the comparator-threshold ↔ LRC-threshold family explicit: for `α = 1/k`
   brackets, the "all-in-arc-`1/k`" cells should be the loneliness cells at
   threshold `1/k`.
3. Signature as resonance invariant: test whether speed sets with the same
   difference-relation lattice share the clock signature.
4. Homology of the walk in `G_n`: does the extremal set's loop generate a
   distinguished cycle of the metagraph?

## Artifacts

```
04-computation/tournament_clock_s24.py             (the clock + invariants)
04-computation/tournament_clock_identify_s24.py    (universality + identification)
05-knowledge/results/tournament_clock_s24.out
05-knowledge/results/tournament_clock_identify_s24.out
```
