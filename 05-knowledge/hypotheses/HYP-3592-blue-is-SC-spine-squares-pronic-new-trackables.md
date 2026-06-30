---
id: HYP-3592
title: BLUE/BLACK tiling patterns -- BLUE(n)=#grid-symmetric tilings = 2^{e}, e = k^2 (odd n=2k+1, SQUARES) / k(k-1) (even n=2k, PRONIC); and the sharp identity BLUE-CLASSES = SC (a class has a grid-symmetric tiling IFF it is self-complementary; pure-black = NS exactly), so blue lives on the SC spine with odd grid-sym multiplicity (= S675b's odd-boundary); the genus-obstruction is the BINDING SUBSET within the blue spine (the doublet at N=14), not the blue count
status: VERIFIED (n=3..9 for BLUE(n)=2^{sq/pronic}; n=3..6 for blue-classes=SC and pure-black=NS, by exhaustive tiling enumeration + canonical iso classes). New invariants/trackables defined.
source: klein-2026-06-29-S13
depends_on:
  - HYP-3591   # the genus is the odd boundary (this corrects/refines it)
related:
  - merged-line-parity-even-odd-s675b   # codex: black=Eulerian, blue=odd-boundary; this gives the exact support (SC) + count (squares/pronic)
  - THM-584    # SC = R-fixed spine (= the blue support)
  - THM-578    # the doublet = the genus-1 binding atom within the blue spine
  - HYP-3586   # cusps=Klein, hardness=genus, nu_2=0<=>Paley (the parity facts the square/pronic alternation watches)
results:
  - 04-computation/blue_black_tiling_patterns_klein.py
  - 05-knowledge/results/blue_black_tiling_patterns_klein.out
---

# HYP-3592 — blue is the SC spine; squares & pronic; new trackables

## Verified patterns

1. **`BLUE(n) = #grid-symmetric tilings = 2^{e(n)}`**, with `e(n) = k^2` for odd `n=2k+1` (SQUARES) and
   `e(n) = k(k-1)` for even `n=2k` (PRONIC). Sequence `BLUE = 2,4,16,64,512,4096,65536` (n=3..9), exponents
   `1,2,4,6,9,12,16`. The blue tilings are the fixed sub-cube `Fix(grid) ≅ Q_{e(n)}` of the anti-diagonal
   staircase reflection.
2. **BLUE-CLASSES = SC, exactly.** `#{iso classes containing >=1 grid-symmetric tiling} = SC(n)` (= 2,2,8,12
   for n=3..6); equivalently **pure-black classes = NS classes exactly**. A class has a grid-symmetric
   tiling IFF it is self-complementary. (Reason: transpose = complement for tournaments; the grid involution
   is the tiling-level transpose; grid-fixed = transpose-fixed = SC. Sharpens CLAUDE.md's "transpose-self
   never pure-black" to the exact count.)
3. Within the SC spine: pure-blue split `= 2,1,3,2` (n=3..6), the rest mixed.

## The correction to HYP-3591

Blue support = SC = `2,2,8,12,88` is NOT the genus `0,0,1,2,2`. So blue is the **arena** (the SC spine,
`R`-fixed), not the obstruction count. S675b's "odd-boundary" = the **odd multiplicity** of grid-sym tilings
per SC class (0 over NS). The genus-obstruction = the **binding subset** of the blue spine (the doublet at
`N=14`); the genus counts independent binding atoms (1 for LRC14), the blue/SC spine is where they sit.
WHERE (blue=SC spine) vs HOW MANY (genus) are different invariants on the same arena.

## New trackables (defined)

- **`BLUE(n) = 2^{e}`**, `e = square/pronic` (parity-alternating). The square side is the LRC-relevant (apex,
  Paley) parity.
- **blue=SC identity** (and pure-black=NS); persistence to track at n>=7.
- **blue multiplicity `mu_diamond(C) = #grid-sym tilings in class C`**: 0 if NS, ODD if SC (the odd-boundary
  spectrum); pure-blue when `mu = #tilings(C)`.
- **black Eulerian cycle rank** (S675b open Q, defined, uncomputed): `b_1` of the black non-grid-sym carrier.
- the **square/pronic parity** as a diagnostic, watched against `nu_2=0<=>Paley` and the genus parity (HYP-3586).

Pairs with mac-mini S31 (the `4cos^2(3pi/7)` cusp value as an EVEN object): the binding doublet atom sits in
the blue/SC spine, and its gap is the cusp constant.
