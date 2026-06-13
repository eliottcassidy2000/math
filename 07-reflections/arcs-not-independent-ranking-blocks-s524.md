---
source: oracle-2026-06-01-S524
status: reflection + verified computation (THM-354 on all tilings n<=6; round=extremes)
tags: [tiling-model, ranking, SCC, THM-354, round-tournament, A000016, A000568, dependence, condensation]
---

# Flipping an arc is not independent: a tournament is a ranking of anti-transitive blocks

**User's thesis.** `2^C(n,2)` treats each arc as an independent switch. But a
tournament is a RANKING (Hamiltonian path) composed recursively of sub-rankings
that either align with the path (transitive) or run against it (anti-transitive).
The tiling model captures this hidden dependence; an arc flip is less independent
than it looks.

This is exactly right, it is a **theorem**, and it explains where the
LRC-realizable necklace (S523) sits.

## 1. The thesis is THM-354

Fix the base Hamiltonian path `0->1->...->(n-1)`. The non-path arcs are the
`C(n-1,2)` **tiles**; a tile `(i,j)` is *aligned* (`i->j`, with the ranking) or an
*upset* (`j->i`, against it). Now the classical fact (THM-354, re-verified here on
**every** tiling for n<=6):

> the strong components of any tournament are intervals of every Hamiltonian
> path, and `#good-cuts = n - #SCC`, where a good cut is a path position crossed
> by a backward (upset) arc.

So the ranking is a **transitive chain of `#SCC` strongly-connected
anti-transitive blocks**, glued at the `#SCC - 1` clean (upset-free) cuts. That is
the user's "recursive sub-rankings, aligned or against the current," precisely.
The arcs are not free: once the block decomposition is set, *between* blocks every
arc is forced aligned, and *within* a block the arcs must keep it strongly
connected. Independence is an illusion of the coordinate choice.

## 2. The hidden dependence, quantified (`tile_dependence_blocks_s524.py`)

Over the full tile cube (base path 0->..->n-1):

- **Collapse.** `2^C(n-1,2)` tilings surject onto **exactly** `A000568(n)` iso
  classes (8->4, 64->12, 1024->56 for n=4,5,6) with wildly uneven fibers
  (n=6: 43 down to 1; fiber size = `H/|Aut|`, CLAUDE.md). If arcs were independent
  axes the classes would not pile up like this — the pile-up *is* the dependence.
- **No arc is a free axis.** For each tile, the fraction of contexts in which
  flipping it changes the iso class is **never 1** (every tile is sometimes a
  *silent* flip) and **varies by position** in a symmetric staircase pattern
  (n=6: corner/skip-1 tiles 0.979, central tiles 0.926). A genuine independent
  switch would be expressive in 100% of contexts, uniformly. It is neither.
- **Most tournaments are one block.** `#SCC` distribution at n=6:
  `{1:903, 2:101, 3:15, 4:4, 6:1}`. The fully transitive tournament (`#SCC=n`)
  is a single tiling; the overwhelming majority are a single strongly-connected
  anti-transitive block (`#SCC=1`). "Anti-transitive against the current" is the
  generic state, not the exception.

## 3. Where the LRC-realizable necklace lives: the block-structure extremes

Connecting to S523 (the LRC clock-movie realizes only ROUND tournaments =
`A000016` necklace). In the tile/block language the round slice is sharp:

- **round-with-base-path count = `(n-1)!`** (6, 24, 120 for n=4,5,6).
- **`#SCC` of a round tournament is `1` or `n` only** (n=6 round: `{1:119, n:1}`;
  n=5: `{1:23, 5:1}`). So:

> **A round (LRC-realizable) tournament is either fully transitive (0 good cuts,
> all aligned) or a single pure strong block (`n-1` good cuts, maximally
> against-the-current as one irreducible piece). It never has an intermediate
> chain of blocks.**

This is THM-374 (round points in a semicircle => transitive; wrapping around =>
strongly connected) read through THM-354. The user's recursive block structure has
a full ladder `#SCC = 1, 2, ..., n` for general tournaments, and **the LRC clock
can only realize the two ends of that ladder.** The whole middle — tournaments
that are genuine transitive chains of several non-trivial anti-transitive blocks —
is exactly what the clock-movie cannot reach at open times. That is another face
of "realizable = a 2% necklace": the necklace is the block-structure boundary.

## 4. The synthesis across the thread

```
arcs as independent switches      2^C(n,2)            (false essence)
  -> ranking + tiles               base path + upsets  (THM-354 coordinates)
       -> blocks                    chain of #SCC strong anti-transitive blocks
            -> LRC-realizable        round = #SCC in {1,n} only = A000016 necklace
                 -> extremal          regular polygon R_{n-1} = roots of unity (S522)
```

Each arrow removes false freedom. The "minor hidden dependence among seemingly
independent variables" the user named is, concretely: arcs are constrained by the
condensation (block) structure; the realizable-by-LRC tournaments are pinned to
the two extremes of that structure; and the symmetric extreme of the strong-block
end is the regular polygon. The conjecture's difficulty is not spread over
`2^C(n,2)` switches — it is concentrated on one necklace of irreducible blocks.

## Open (→ HYP-2000)
- Prove `round <=> #SCC in {1,n}` for all n (n<=6 verified); and
  `round-with-fixed-path = (n-1)!`.
- The intermediate `#SCC = 2,...,n-1` tournaments are LRC-*unrealizable* at open
  times: do they ever appear as BOUNDARY (tie-resolved) lonely classes? If not,
  the lonely menu is itself confined to the block extremes.
- Per-tile silent-flip rate as a function of staircase position: is the symmetric
  profile (0.979 corners ... 0.926 center) a known wiggly-class statistic
  (CLAUDE.md), and does it have a closed form?

## Anchor
`04-computation/tile_dependence_blocks_s524.py` (+ `.out`): THM-354 on all tilings
n<=6; collapse to A000568; per-tile context-dependence; round = `#SCC in {1,n}`,
count `(n-1)!`. Builds on THM-354, HYP-1998 (S523), HYP-1995 (S522).
