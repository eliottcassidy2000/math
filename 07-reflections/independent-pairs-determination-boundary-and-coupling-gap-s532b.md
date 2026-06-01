---
source: oracle-2026-06-01-S532b
status: complement to S532 (the iso-determination boundary; coupling gap; channel gcd-decomposition)
tags:
  - lonely-runner
  - independent-pairs
  - matching
  - a000568
  - coupling
  - channels
  - modular-decomposition
---

# Independent Pairs: the Determination Boundary (n≤4) and the Coupling Gap

The user's metric for the multi-channel generalization — *the amount and state of
independent pairs* — is exactly right, and it has a sharp boundary. An **independent
pair** is an edge of a matching (two vertex-disjoint arcs). The concurrent S532
(`independent_pairs_channels_s532.py`, HYP-2012) verified the user's `n=4` claim
(with a fixed frame, the 2 matching-arc bits index all 4 iso classes; 8 of 16 frames
give the bijection), tied the count `⌊n/2⌋` to opus-S524's `n=14` CRT (7 pairs), and
showed the `n=4` parity law is the one-channel fusion of the 2 pairs. This note adds
the **boundary** and a **quantitative coupling measure** (and does not reduplicate).

## Independent pairs determine the iso class — exactly iff n ≤ 4

Flipping a set `S` of arcs with the rest fixed can index all `A000568(n)` iso
classes only if `2^|S| ≥ A000568(n)`. For *independent* arcs (a matching),
`|S| ≤ ⌊n/2⌋`. So independent pairs alone can determine the iso class only if
`2^{⌊n/2⌋} ≥ A000568(n)` — and (`independent_pairs_boundary_s532b.py`):

```
n :  A000568   2^floor(n/2)   indep suffice?   coupling gap = ceil(log2 A000568) - floor(n/2)
3 :     2          2              YES                 0
4 :     4          4              YES                 0
5 :    12          4              no                 +2
6 :    56          8              no                 +3
7 :   456          8              no                 +6
8 :  6880         16              no                 +9
```

> **Independent pairs index the entire iso-class set exactly iff `n ≤ 4`** (equality,
> coupling gap `0`). Verified directly at `n=5`: the best 2-arc matching over the best
> frame reaches only `4` of the `12` classes — the `2^2` ceiling — so the other `8`
> require **coupled** (vertex-sharing) arcs.

The **coupling gap** `ceil(log2 A000568(n)) − ⌊n/2⌋` (`0,0,2,3,6,9` for `n=3..8`) is a
new, clean quantity: the number of genuinely-coupled bits *beyond* the independent-pair
channels. It is `0` precisely at the cases (`n≤4`) where the structure is purely
"independent," and grows thereafter — the multi-channel content made numerical.

## Why n=4 is the boundary: the diameter channel is born there

The skip-shell channels of the regular `n`-gon decompose as
`skip j → gcd(n,j) cycles of length n/gcd(n,j)`. A channel is a **perfect matching**
(maximally independent — disjoint arcs) iff its cycle length is `2`, i.e. `j = n/2`,
i.e. the **diameter shell of even `n`**:

```
n=4: j2 = MATCHING (the 2 diameters)         <- first diameter channel
n=6: j3 = MATCHING (3 diameters);  j2 = 2×C3
n=7: j1,j2,j3 all = single C7                <- ODD n: NO matching channel at all
n=8: j4 = MATCHING (4 diameters)
```

So the user's "`2` independent pairs" at `n=4` are exactly the **2 diameters**, and
`n=4` is the *smallest even `n` with a diameter (perfect-matching) channel*. Odd `n`
has no matching channel — every shell is one or more odd cycles, all coupled. This is
the same boundary as "the inside is born at `n=4`" (S529): the first genuine interior
chord is the diameter, and the diameter is the first independent-pair channel.

## The synthesis: independent pairs = the factoring channels

Independent (vertex-disjoint) arcs are precisely the **factoring** degrees of freedom.
S531 proved `H` is *multiplicative over vertex-disjoint modules* and *coupled along
nesting*; here the same dichotomy reappears at the level of iso-class determination:

- **independent pairs (matching, `⌊n/2⌋`)** = the free / factoring channels — they fully
  pin the structure only through `n=4`;
- **the coupling gap (`≥ 5`)** = the nested/shared-vertex part — the *inside debt* (S529),
  the *nested coupling* (S531), the genuinely multi-channel obstruction.

So "the amount and state of independent pairs" is the correct coordinate for the *free*
part of the multi-channel problem, and the conjecture's difficulty is exactly its
complement — the coupling gap, which vanishes only for `n ≤ 4`. The whole arc is
consistent: LRC is easy where the structure is all independent pairs (the diameter
matching, `n ≤ 4`), and hard exactly where coupled channels appear and the gap opens.

## Verdict / next
- New: the **determination boundary** (independent pairs suffice iff `n≤4`) and the
  **coupling gap** `ceil(log2 A000568) − ⌊n/2⌋` (0,0,2,3,6,9), with the channel
  gcd-decomposition showing the diameter matching is born at `n=4` (even `n` only).
- Converges with S532/HYP-2012 (independent pairs = `⌊n/2⌋` = `n=14` CRT) and S531
  (independent = factoring) and S529 (coupling = inside debt).
- Concrete next: is the coupling gap `ceil(log2 A000568) − ⌊n/2⌋` an *achievable*
  lower bound (can `⌊n/2⌋` independent + (gap) coupled arcs really index all classes)?
  And does the LRC inside debt's resonance order match the coupling gap channel-count?

## Artifacts
```
04-computation/independent_pairs_boundary_s532b.py
05-knowledge/results/independent_pairs_boundary_s532b.out
```
Related: S532 / HYP-2012 (independent-pairs = channels), S531 (modular H / independent
= factoring), S529 (inside debt), S530 (diameter pair = speed n/2), opus-S524 (n=14 CRT).
