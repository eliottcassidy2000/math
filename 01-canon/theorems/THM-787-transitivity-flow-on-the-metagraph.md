---
id: THM-787
title: The flow of transitivity on the metagraph — the transitive node drains through a SINGLE blue line; blue lines never touch pure black; the blue Δx parity law with max exactly 8(n−2); the black sea carries the downhill flux
status: PROVED (blue incidence, unique transitive line, and THM-785 promotion of blue parity/max laws) + VERIFIED-EXACT n = 3..7 (A000568-certified class census and black concentration)
source: opus-2026-07-14-S304 (owner directive: trace the flow of transitivity from the strict ordering to the circulant; quantify the left-right symmetry and imbalance)
depends_on: []
related: [HYP-6850, THM-781, THM-785, the geometric-alignment frame (spine/ribs/sea), merged-metagraph-invariants]
verification: 04-computation/metagraph_transitivity_flow_opus_S304.py
  (+ 05-knowledge/results/metagraph_transitivity_flow_opus_S304.out)
---

# THM-787 — the flow of transitivity on the metagraph

> **Concurrent completion by THM-785.**  THM-785 proves
> `x=n(n^2-1)/3-8C3`, the exact line flux
> `Delta C3=d0-dlast-1`, and the closed blue step distribution for every size.
> Consequently the parity and maximum laws conjectured in the original
> version of §3 are now theorems for all `n`.  THM-785 also performs the
> converse-merged **node** census and directional category-flow audit.  The
> phase histograms below count unmerged classes; for example, 368 pure-black
> classes at `n=7` form 184 converse-merged pure-black nodes.  This is a level
> distinction, not a numerical conflict.

**The axis.** For an iso class with score multiset (s_i), set
x = Σ_i (2s_i − (n−1))² — an integer, maximal at the transitive class
(x_max = Σ_i (2i − (n−1))², i = 0..n−1), minimal at the maximally distributed
score class (0 at odd n, n at even n), invariant under reversal — so it
descends to the merged metagraph G_n/Z_2. All x-levels lie in one step-8
arithmetic progression, and EVERY line (both colours, all n ≤ 7) has
Δx ≡ 0 (mod 8). Classification of all 2^{C(n−1,2)} tilings certified exact by
the A000568 counts (2, 4, 12, 56, 456) and the fiber×Aut = H identity.

## (1) Blue lines never touch pure black (PROVED, one line)

A blue line's tiling t and its full-flip t̄ are both grid-symmetric
(isGridSym(t̄) = isGridSym(t)), so BOTH endpoint classes contain a
grid-symmetric tiling — neither can be pure black. ∎
(Verified: the only blue endpoint-type pairs occurring are
(pureblue, mixed) and (mixed, mixed) for n ≥ 4; (pureblue, pureblue) occurs
only at n = 3.)

## (2) The transitive node drains through a SINGLE line, and it is blue (PROVED)

The transitive class has exactly one tiling (fiber 1 = H/|Aut| = 1/1... its
fiber is a single tiling: the all-forward tiling), which is grid-symmetric; so
the transitive node carries exactly ONE line, necessarily blue, landing at
x_max − 8(n−2) — verified: 20→4, 40→16, 70→38, 112→72 for n = 4..7 — on a
MIXED class. **The strict ordering's entire connection to the rest of the
metagraph at the d = m layer is one blue pipe of drop exactly 8(n−2).**

## (3) The blue Δx parity law (PROVED for all n by THM-785)

Blue |Δx| spectra: n=3: {8}; n=4: {0,16}; n=5: {8,24}; n=6: {0,16,32};
n=7: {8,24,40}. Hence:

> **Odd n: every blue line has |Δx| ≡ 8 (mod 16) — blue lines are NEVER level.
> Even n: every blue line has |Δx| ≡ 0 (mod 16) — level blue lines exist.
> In both cases max |Δx| = 8(n−2), attained (among others) by the transitive
> pipe.** Blue |Δx| histograms: n=5: {8:6, 24:2}; n=6: {0:12, 16:16, 32:4};
> n=7: {8:160, 24:80, 40:16}.

This is the owner's predicted LEFT-RIGHT SYMMETRY made exact: the blue
spectrum is rigidly quantized in an arithmetic progression of step 16
 (three values at `n=5,6,7`, with the number growing thereafter), and at even
`n` the level value Δx = 0 makes blue lines mirror-glue
equal-x classes, while at odd n blue lines always TRANSPORT (never stay level)
— the parity of n flips the blue geometry between "mirror" and "conveyor."

## (4) The black imbalance (VERIFIED-EXACT; the owner's predicted asymmetry quantified)

- **Type flow:** black lines by endpoint types at n=7:
  (mixed, mixed) 876; (mixed, pureblack) 5,044; (pureblack, pureblack) 10,208 —
  the black SEA. At n=6: 10 / 156 / 314. Pure black classes receive transitivity
  ONLY through black lines (by (1)) — the owner's "then along black lines into
  pure black" is exact and exclusive.
- **Level fraction:** black lines are level (Δx = 0) at a stable ≈ 22% rate
  (108/480 at n=6; 3,584/16,128 at n=7; exactly 1/3 at n=5) — blue at odd n: 0%.
- **Where the types live:** pure blue sits at the TOP of the axis (n=7:
  x ∈ {112, 104, 96, 72}); mixed concentrates LOW (n=7 histogram peaks at
  x = 16 with 17 classes and spans down to x = 0); pure black dominates the
  MIDDLE (n=7 peak at x = 32 with 74 classes). The flow picture: blue pipes
  hand transitivity from the top into the mixed layer; the black sea then
  carries it across the middle into the distributed end.
- **The downhill flux concentrates in the late middle, toward the distributed
  end:** black downhill
  line counts by level-pair peak at (32→24): 1,504 and (40→32): 1,206 at n=7 —
  an order of magnitude above the fluxes near the transitive end
  ((104→72): 2; (96→64): 8).  THM-785's exact cut-current refinement shows
  that this is not monotone all the way to the endpoint: in `C3` coordinates
  the current rises from 4 at cut `1|2` to 4,290 at `10|11`, then falls to 70
  at `13|14`.  “Late-middle concentration” is therefore the precise claim;
  “increases at every cut” would be false.
- **An n-parity deviation worth recording:** at n = 5 the distributed end
  (x = 0) is PURE BLUE; at n = 7 the x = 0 classes are all MIXED and pure blue
  stops at x = 72. The circulant end is blue-pure only at small odd n — the
  quasi-regular classes at n = 7 already carry non-grid-symmetric tilings.

## (5) Summary tables (exact, n = 3..7)

| n | classes | SC | merged | pure blue | mixed | pure black | blue lines | black lines |
|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 2 | 2 | 0 | 0 | 1 | 0 |
| 4 | 4 | 2 | 3 | 1 | 1 | 2 | 2 | 2 |
| 5 | 12 | 8 | 10 | 3 | 5 | 4 | 8 | 24 |
| 6 | 56 | 12 | 34 | 2 | 10 | 44 | 32 | 480 |
| 7 | 456 | 88 | 272 | 4 | 84 | 368 | 256 | 16,128 |

(Blue lines = `2^(r-1)` for
`r=(m+floor((n-1)/2))/2`: 1, 2, 8, 32, 256.  The anti-diagonal reflection has
`floor((n-1)/2)` fixed tiles, so `r` is its tile-orbit count.  This corrects
the original shorthand `2^(ceil(m/2)-1)`, which fails already at `n=5`.)

## (6) Next checks (named)

- n = 8 (6,880 classes; 2^21 tilings): THM-785 predicts the full blue spectrum
  `{0:640,16:960,32:384,48:64}`; the invariant-certification method scales
  (A000568(8) = 6,880 validates).
- **Delivered by THM-785:** blue reflection gives `Delta C3=2d0-n`; its exact
  binomial support has maximum `n-2` and parity `n`.  Multiplication by eight
  proves the stated `Delta x` parity and maximum without a finite-size bound.
- The merged-metagraph version of the flux tables, and the spine/ribs/sea
  overlay (SC-SC vs SC-NS vs NS-NS refined by line colour).
