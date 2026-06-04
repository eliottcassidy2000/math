---
source: monad-compute-2026-06-04-S4
status: EXHAUSTIVE closure of the HYP-2193 finite window for H=21 (compute result)
tags: [H21, H-impossibility, moon-theorem, strong-tournament, finite-check, exhaustive, isomorph-free, canonical-augmentation]
---

# Closing the H=21 finite window by exhaustion (m = 9, 10, 11, 12)

## The setup I inherited

HYP-2193 (claude-2026-06-03-S617) reduced the last open permanent H-gap, `H≠21`,
to a **finite** check:

1. `H` is multiplicative over strong components; `21 = 3·7` and `7` is **not** a
   strong H-value (THM-029) ⟹ if `H(T)=21` then some **strong** component `S`
   has `H(S)=21`.
2. `H(S)=21 = I(Ω,2) = 1 + 2α₁ + 4α₂ + …` ⟹ `α₁ ≤ 10`; 3-cycles are odd cycles
   ⟹ `c₃(S) ≤ α₁ ≤ 10`.
3. Moon (1968): a strong tournament on `m` vertices has `c₃ ≥ m−2` ⟹ `m ≤ 12`.
4. THM-079 Part G already killed `H=21` **exhaustively for m ≤ 8** (2^28 at n=8).

So the *entire* remaining obligation was: **strong tournaments on
`m ∈ {9,10,11,12}` with `c₃ ≤ 10`.** S617 left this as "exhaust, or prove
`α₁ ≥ 11` for strong `m ≥ 9`", with only *sampling* evidence (min H = 75, 153 at
m = 9, 10).

## What this session did: exhausted it

This node has no nauty / no C compiler / no numba — pure Python. But the window
is small once you use the right pruning:

- **Monotone pruning.** `c₃` never decreases when you add a vertex (a new vertex
  only creates new triangles), so *every induced subtournament* of a `c₃≤10`
  tournament is itself `c₃≤10`. Build vertex-by-vertex and prune partial
  `c₃ > 10`; nothing is lost.
- **Isomorph-free generation.** Generate iso-classes level by level: from each
  rep on `k` vertices try all `2^k` orientations of a new vertex, keep those with
  `c₃≤10`, canonicalize, dedup. Validated against **A000568** (1,1,2,4,12,56,456
  for n=1..7) with the cap removed — exact match, so the canonical form and the
  generation are correct.
- **Cheap canonicalization.** Low-`c₃` tournaments are *near-transitive*, so the
  score sequence is nearly `0,1,…,m−1` and colour-refinement gives singleton
  classes almost always ⟹ the canonical labelling is essentially forced.

The iso-class counts with `c₃≤10` grow slowly (the constraint tightens as the
mean `c₃ = C(m,3)/4` climbs): 339 (m=7), 1,105 (m=8), 2,575 (m=9), 5,277 (m=10),
9,989 (m=11), [m=12 …]. For each **strong** class I computed `H` by Held-Karp.

## Result

| m | iso-classes c₃≤10 | strong | H=21? | min H |
|---|---|---|---|---|
| 8 | 1,105 | 467 | none | 45 |
| 9 | 2,575 | 605 | **none** | 75 |
| 10 | 5,277 | 709 | **none** | 125 |
| 11 | 9,989 | 560 | **none** | 225 |
| 12 | 17,947 | 256 | **none** | 375 |

(m≤8 reproduced as a cross-check against THM-079 Part G.) **`H=21` occurs
nowhere in the entire finite window.** The exhaustive min `H` (75, 125, 225, 375
for m=9,10,11,12) sharpens S617's sampling (min H = 75 at m=9) and grows fast:
every strong `c₃≤10` tournament has `H` *far* above 21 — because a strong
tournament on `m ≥ 9` is vertex-pancyclic (Moon), so it carries many long odd
cycles (5,7,9-…) on top of its triangles, pushing `α₁ ≫ 10`. Note at m=12 the
strong count *drops* (256) because `c₃ ≤ 10` and Moon `c₃ ≥ m−2 = 10` force
`c₃ = 10` exactly — only the Moon-minimal strong tournaments survive.

Two independent engines agree: v1 (`h21_finite_check_monad_s4.py`, direct
`c₃`-count pruning) and v2 (`h21_finite_check_v2_monad_s4.py`, DFS incremental-
`c₃` pruning, ~10× faster) give identical counts for m ≤ 10; both reproduce
A000568 with the cap removed. v2 carried the run to m=12.

## Why it finishes the proof

`H=21` requires `α₁ ≤ 10`, hence `c₃ ≤ 10`. So the *only* strong tournaments
that could possibly give `H=21` on `m ∈ {9,…,12}` are exactly the ones
enumerated here — and none does. Combined with THM-079 Part G (`m ≤ 8`),
H-multiplicativity over strong components, THM-029 (`7` non-strong), and Moon's
`c₃ ≥ m−2`, this closes the finite window: **`H(T) ≠ 21` for all tournaments.**
With THM-343 (`H≠7`), `{7, 21}` is the complete set of permanent H-gaps.

## Meta

The pattern that recurs in this project: a "search over all tournaments on
`n ≥ 9`" that looks hopeless (2^36–2^66) becomes a few-thousand-class
enumeration once a *structural* constraint (`c₃ ≤ 10`, near-transitivity) is
pushed into the *generation* rather than applied as a post-filter. The Moon
bound does the heavy lifting: it converts an asymptotic question into a bounded
one, and near-transitivity makes the bounded set tiny.
