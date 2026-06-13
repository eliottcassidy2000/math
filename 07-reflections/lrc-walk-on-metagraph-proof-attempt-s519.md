---
source: oracle-2026-06-01-S519
status: proof attempt (walk-on-metagraph for n=14) — honest negative + the real obstruction
tags:
  - lonely-runner
  - n14
  - tournament-metagraph
  - feedback-vertex-set
  - proof-attempt
---

# Walk-on-the-Metagraph: an Honest Attempt at n=14, and Where It Really Stalls

Goal: use "LRC@n = a closed runner-walk in the tournament metagraph must reach a
lonely class" (S511/S518) to attack n=14, reframing the mapping if needed. This is
an honest report: the clean reframe **fails**, and pinning down *why* locates the
irreducible content of LRC precisely.

## The faithful framing

A runner system gives, as `t` runs over `[0,1)`, a piecewise-constant tournament
`T(t)` that changes by ONE arc-flip at each wall — a **closed, t-ordered walk** in
the menu metagraph `M_n` (vertices = the `2·Fib(n-2)` circular classes, S518; edges
= single wall-crossings). LRC@n ⟺ **every realizable closed walk hits a lonely
class** (transitive at threshold `½`; the marked-source class at the exact `1/n`).

## Reframe attempt: lonely classes as a feedback vertex set

Sufficient condition for "every closed walk hits lonely": the lonely classes are a
**feedback vertex set** (FVS) — removing them makes `M_n` acyclic. Then no cycle,
hence (one hopes) no closed walk, can avoid them. Computed `M_n` and tested
(`lrc_menu_metagraph_fvs_s519.py`):

```
n=5: M_5 = 4 vertices, 3 edges  -> a TREE (no cycles); transitive a leaf.
n=6: M_6 = 6 vertices, 7 edges  -> 2 independent cycles; transitive NOT a FVS.
n=7: M_7 = 10 vertices,13 edges -> 4 independent cycles; transitive NOT a FVS.
```

So for `n ≥ 6` the half-turn "transitive" (the `½`-gap / bunched-in-a-semicircle
target) is **not** a feedback vertex set — there are cycles in `M_n` avoiding it.
This matches reality: not every runner system ever bunches into a semicircle
(evenly-spread systems don't), so the `½`-threshold target genuinely isn't always
reached. The coarse resolution cannot prove LRC this way.

## Why the reframe does not trivialize LRC (the real obstruction)

Two things kill the pure-graph approach, and naming them is the point:

1. **Closed walks are too free.** A closed *walk* (unlike a simple cycle) may
   backtrack: on any graph, a closed walk can avoid a set `S` as long as `V∖S`
   contains a single edge (walk back and forth on it). So "every closed walk hits
   `S`" is almost never a graph property — even on the `n=5` tree a closed walk can
   sit in a subtree away from transitive. The FVS criterion is therefore *not*
   equivalent to "every closed walk hits lonely."

2. **Realizable walks are arithmetically special.** The runner walk is not an
   arbitrary closed walk: it is the **t-ordered cyclic sequence of cells**, with
   each edge `(i,j)` flipping exactly `2|v_i−v_j|` times (the holdback count, S25),
   in the precise order dictated by the speeds. The realizable walks are a
   measure-zero, arithmetically-constrained family inside all closed walks.
   **Characterizing which closed t-walks are realizable is equivalent to LRC.**

So the metagraph picture is *faithful and clarifying* but **not reducing**: it
moves the difficulty from "find a lonely time" to "show the realizable t-ordered
closed walk can't dodge the lonely class," and that realizability constraint is
the whole conjecture. The hoped-for "small mapping tweak → finite graph property"
is not available through FVS/acyclicity, because the obstruction is the *walk
family*, not the *graph*.

## What survives, and the honest next lever

- **A genuine new object:** the circular-menu metagraph `M_n` — small, explicit
  (`2·Fib(n-2)` vertices), with computed structure (tree at `n=5`, `2,4,…` cycles
  after). Worth its own study (its cycle space = the "non-bunching" walks).
- **The one place the walk is pinned to arithmetic:** the orbit *always* visits
  the regular-polygon times `t = a/n` (it is a continuous loop over `[0,1)`), and
  there loneliness ⟺ **no speed divisible by `n`** — exactly the sieve (THM-369).
  A counterexample dodges only by carrying a multiple of `n` (and of every
  `m≤n`). So the walk meets the sieve precisely at the `n`-gon cells; the residual
  question is the *fine* (`q>n`) regime between those cells — the coarse/fine split
  of S18.
- **The real reframe that could work** is not graph-topological but
  *order-combinatorial*: model the realizable walk as the **interval-exchange /
  circular wall-crossing order** of the speed differences, and show that order
  cannot keep the observer non-source for a full lap. That is the bounded-ansatz
  program (S514) in disguise — the cofactor bound on `t = j/(14 s)` is exactly a
  bound on how fine the wall-order can get.

## Verdict
The walk-on-the-metagraph concept is the right *language* (it makes LRC a
reachability statement on a tiny, fully-understood set of shapes), but it does
**not** by itself prove n=14: the difficulty is intrinsic and lives in *which
closed walks the runner arithmetic realizes*, equivalent to LRC. The honest path
forward couples this language with the sieve (at the `a/n` cells) and the
bounded-ansatz/interval-exchange structure (between them) — not with a pure
metagraph property.

## Artifacts
```
04-computation/lrc_menu_metagraph_fvs_s519.py
05-knowledge/results/lrc_menu_metagraph_fvs_s519.out
```
