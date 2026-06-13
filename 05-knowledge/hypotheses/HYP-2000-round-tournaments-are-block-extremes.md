---
id: HYP-2000
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S524
related:
  - HYP-1998
  - HYP-1995
  - THM-354
  - THM-374
---

# HYP-2000: arcs are not independent — LRC-realizable (round) tournaments are the extremes of the SCC block ladder

**Frame (user's thesis, = THM-354).** A tournament is a ranking (Ham path) =
transitive chain of `#SCC` strongly-connected ANTI-transitive blocks; in the
tiling model (fixed base path, non-path arcs = tiles, aligned vs upset),
`#good-cuts = n - #SCC`. Arcs are not independent switches.

**VERIFIED (`tile_dependence_blocks_s524.py`, n=4,5,6):**
1. THM-354 holds on ALL `2^C(n-1,2)` tilings (assert passed).
2. The tile cube surjects onto exactly `A000568(n)` iso classes (8->4, 64->12,
   1024->56) with uneven fibers = `H/|Aut|` -> the collapse IS the dependence.
3. Per-tile arc-flip is NEVER expressive in 100% of contexts (every tile is
   sometimes silent) and varies by staircase position (n=6: 0.979 corners ...
   0.926 center) -> no arc is a free axis.
4. `#SCC` over all tilings concentrates at 1 (n=6: 903/1024 strongly connected);
   transitive (`#SCC=n`) is a single tiling.

**KEY CLAIM (verified n<=6, conjectured general):**
> A ROUND (= LRC-realizable, S523) tournament has `#SCC in {1, n}` ONLY: it is
> fully transitive (0 good cuts, all aligned) or a single pure strong block
> (`n-1` good cuts). It never has an intermediate chain of blocks. Also
> `#{round tournaments containing a fixed Hamiltonian path} = (n-1)!`
> (6, 24, 120 for n=4,5,6).

This is THM-374 (round + semicircle => transitive; else strongly connected) in
block language. Consequence: the general SCC ladder is `#SCC = 1,...,n`, and the
LRC clock-movie realizes ONLY the two ends. The `A000016` realizable necklace
(S523) = the block-structure boundary; the intermediate multi-block tournaments
are LRC-unrealizable at open times.

**OPEN:**
- (A) Prove `round <=> #SCC in {1,n}` for all n, and `round-with-fixed-path=(n-1)!`.
- (B) Do intermediate `#SCC = 2..n-1` tournaments ever appear as BOUNDARY
  (tie-resolved) lonely classes (S523 boundary layer)? If never, the lonely menu
  is confined to the block extremes.
- (C) Closed form for the per-tile silent-flip profile by staircase position
  (relate to CLAUDE.md wiggly-class self-loop rates).

**Why it matters.** Pins LRC's difficulty off the `2^C(n,2)` switch picture and
onto a single necklace of irreducible (strongly-connected) blocks + the transitive
point. Reflection: `07-reflections/arcs-not-independent-ranking-blocks-s524.md`.
