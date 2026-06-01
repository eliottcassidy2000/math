---
source: oracle-2026-06-01-S518
status: synthesis (LRC@n as a tournament; the Fibonacci menu of circular iso-classes)
tags:
  - lonely-runner
  - tournament-analysis
  - a000568
  - circular-tournaments
  - fibonacci
  - isomorphism-classes
---

# LRC@n as a Tournament, and the Fibonacci Menu of its Iso-Classes

The question: *LRC at a fixed `n` is really a question about a tournament of some
size and the structure of its possible isomorphism classes.* Here is that
statement, made precise, with the iso-class set characterized — and it turns out
to be **Fibonacci-sized**.

## The model: a tournament on `n` vertices

A runner system (observer at speed 0 plus `n-1` integer speeds, threshold `1/n`)
is lifted, at each time `t`, to a tournament on **`n` vertices** (observer +
runners) by the half-turn rule on their circle positions. As `t` runs over
`[0,1)` this tournament is piecewise constant — a **closed walk in the metagraph
`G_n`** (the iso-class space counted by `A000568(n)`). LRC@n is a **reachability**
question for that walk:

- coarse (threshold `1/2`): does the walk reach the **transitive** class (all
  runners in a semicircle = an empty half-circle)?
- exact (threshold `1/n`, S511): with the observer marked and its edges drawn by
  the LRC arcs, does the walk reach a **source** class (observer beats every
  runner = observer lonely)?

So **LRC@n = "a runner-walk in `G_n` must reach the loneliness target,"** a finite
combinatorial statement on tournaments of size `n` (equivalently `n-1` runners
with a marked target).

## The possible iso-classes are FEW — and Fibonacci-many

The runners can never realize an arbitrary tournament. The half-turn rule on
circle points produces only the **round / circular tournaments** (every
out-neighbourhood is an arc). Computing the realizable iso-classes exactly
(`lrc_tournament_model_isoclasses_s518.py`):

```
m (= #points):     3    4    5    6    7    8
realizable:        2    2    4    6   10   16     (all exactly verified)
A000568(m):        2    4   12   56  456 6880
```

The realizable counts satisfy `a(m) = a(m-1) + a(m-2)` (verified through `m=8`:
`16 = 2·Fib(6)`) — i.e.

> **# circular tournament iso-classes on `m` vertices  =  2·Fib(m-2)**
> (2, 2, 4, 6, 10, 16, 26, …).

So the set of shapes a runner system can take grows like the **golden ratio**
`φ^m`, an utterly vanishing fraction of `A000568(m) ~ 2^{C(m,2)}/m!`. Despite the
infinitely many speed sets, LRC@n lives on a **Fibonacci-sized stage** of
tournament shapes; the conjecture is entirely about the *walk* through these few
shapes, not about which shapes exist.

(The Fibonacci count is natural: a circular tournament is fixed by its cyclic
"reach sequence" — how far each vertex sees within the leading semicircle — and
the admissible cyclic reach sequences are exactly compositions with an
adjacency/overlap constraint, the classic Fibonacci recurrence; the factor 2 is
the boundary/orientation degree of freedom.)

## The structure of the menu: a graded spread-chain, transitive → regular

The menu is not a featureless set; it is **graded by spread**, read off the
`H`-value (directed Hamiltonian-path count) / the 3-cycle count:

```
m=5:  H = 1,  9, 11, 15            (4 classes)
m=6:  H = 1, 17, 23, 41, 45        (6 classes)
m=7:  H = 1, 33, 47, 51, 105, 123, 137, 151, 175   (10 classes)
```

- **Bottom** = the **transitive** tournament (`H=1`): all runners bunched in a
  semicircle — a `½`-gap — the loneliness configuration.
- **Top** = the **regular** tournament (max `H`), which exists **iff `m` is odd**
  (m=3→`H=3`, 5→`15`, 7→`175`; for even `m` the top is near-regular). The regular
  tournament is the runners at the **regular `m`-gon** — exactly the LRC tight
  witness (S517). So the even composite frontiers land at the top of the menu:
  `n=14 → R_13`, `n=18 → R_17`.

So the menu is a small graded poset from "bunched/lonely" (transitive) to
"evenly spread/symmetric" (regular `m`-gon), and the runner orbit oscillates along
it. **LRC = the walk cannot avoid the bottom (loneliness) end.**

## The exact (`1/n`) target inside the menu

For the true threshold, the loneliness target is the marked **source** classes
(observer out-degree `n-1`). These are tournaments on the `n-1` runners (the
observer free) — nominally `A000568(n-1)`-many — but the *reachable* ones are the
circular sub-tournaments of runners confined to the safe arc of length `1-2/n`
(S512): again a Fibonacci-scale subset, `1, 2, 6, 6, …` for `n=4..7`. So even the
win-set is a tiny circular menu, not the full `A000568(n-1)`.

## The one-sentence answer

> **LRC@n is the statement that a closed runner-walk in `G_n` — moving among the
> `2·Fib(n-2)` circular tournament iso-classes (a graded chain from the transitive
> "bunched" class to the regular `n`-gon) — must reach the loneliness end.** The
> tournament size is `n`; the possible iso-classes are the round/circular ones, a
> Fibonacci-sized, spread-graded subset of `A000568(n)`; and the conjecture is a
> reachability fact about the walk, not about the (tiny, fully-understood) set of
> shapes.

## Next
- Prove `a(m) = 2·Fib(m-2)` from the reach-sequence/composition characterization
  (and the `m=8,9 = 16,26` checks).
- Identify the menu's poset/lattice structure (Hasse edges = single wall-crossings
  = `G_n` edges between adjacent circular classes).
- Tie the "cannot avoid the transitive/source end" to the menu's graph structure:
  is the loneliness end a *dominating* target for every realizable walk?

## Artifacts
```
04-computation/lrc_tournament_model_isoclasses_s518.py
05-knowledge/results/lrc_tournament_model_isoclasses_s518.out
```
