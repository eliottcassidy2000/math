---
id: HYP-9540
title: "Packing-number line laws of the tiling metagraph: (A) at every order the triple words of tight lines (nu(t) + nu(t+1) = floor(n/3)) are exactly the tight colourings; (B) black lines join nu-classes i <= j iff i >= 1 and i + j >= floor(n/3); blue lines join the same pairs except (1, 1) from n = 6 on, plus the transitive line (0, floor(n/3)); (C) the triple words of BLUE tight lines are exactly the palindromic tight colourings (so at n = 3k + 2 the blue boundary pairs are the {j, k - j} with j even or k odd); (D) the blue-line (1,3) excess of THM-4951 (5) persists at every odd n >= 11."
status: >
  OPEN.
  (A): PROVED at n = 3k + 2 (THM-4951 (3)); FINITE-EXACT for n <= 12; realizability checked at n = 13 (first-hit,
  pre-registered T5 CONFIRMED).
  (B): FINITE-EXACT for n <= 9 (complete census). The boundary i + j = k is realised at every n (THM-4951 (2)), and
  the pairs (0, k) and (1, k) at every n >= 4 (the transitive line, and the single reversed tile (0, 3)).
  (C): PROVED at n = 3k + 2 (THM-4951 (3d)); FINITE-EXACT for n <= 12.
  (D): one data point (n = 11, (1,3)-share 0.52148). OPEN.
session: opus-2026-11-09-S87
source: 05-knowledge/results/ternary_metagraph_20261109/README.md
related:
  - THM-4951
  - THM-643
  - THM-4922
  - THM-791
---

# HYP-9540 - packing-number line laws of the tiling metagraph

**Setting.** As in THM-4951:
- nu is the maximum number of vertex-disjoint cycles, and k = floor(n/3).
- The triple word of a tiling is x_a = [tile (a, a + 2) reversed], a = 0..n-3.
- A colouring of [0, n-3] is tight if P(A) + P(B) = k, with P the largest 3-separated subset.

## (A) Word law

**Statement.** For every n >= 3, a word is the triple word of a tight line iff it is a tight colouring.
- One direction is proved at every n (THM-4951 (1)).
- The converse is proved at n = 3k + 2 (multi-staircase).
- The converse is FINITE-EXACT for n <= 12 (50 of 50 words at n = 12), and every colouring is realised at n = 13 (20
  of 20).

**Why it should hold.** The tight colourings at n = 3k + 1 appear to be the block words for block decompositions of
[0, 3k - 2] into k - 1 blocks of 3 and one block of 2. Their number is (k + 1) 2^(k-2), checked for k <= 13. A
multi-staircase with one window of 4 vertices should realise each. n = 3k allows more local patterns, such as 1010
at n = 6.

**Falsifier.** A tight colouring at some n = 3k or 3k + 1 that no tight line realises.

## (B) Allowed pairs

**Statement.** For n >= 4 and nu-classes i <= j:
- some BLACK line joins them iff i >= 1 and i + j >= k;
- some BLUE line joins them iff (i, j) = (0, k), or i >= 1, i + j >= k and, when n >= 6, (i, j) != (1, 1).

The pair (0, k) is the transitive line alone, which is blue.

**What is known.**
- Necessity is proved: Session26's bound, the uniqueness of nu = 0, and the blue (1, 1) exclusion (THM-4951 (1)).
- Sufficiency is proved on the boundary i + j = k (THM-4951 (2a)).
- (1, k) is proved at every n >= 4. Reverse the single tile (0, 3): every cycle uses 3 -> 0, so nu(t) = 1, while
  t + 1 keeps all consecutive 3-cycles on {3b, 3b + 1, 3b + 2}.
- The census n <= 9 has exactly these pairs in each colour.

**Falsifier.** An allowed pair with no line at some n >= 10.

## (C) Blue word law

**Statement.** For every n >= 3, a word is the triple word of a BLUE tight line iff it is a palindromic tight
colouring. Then the blue tight nu-pairs are the pairs {P(A), P(B)} over palindromic tight colourings.

**What is known.**
- Necessity is proved: sigma reverses the triple word.
- Sufficiency is proved at n = 3k + 2 (THM-4951 (3d)). There the blue tight pairs are exactly {j, k - j} with j even
  or k odd, so for example no blue (1, 3) tight line exists at n = 14.
- Even j is realised at every n >= 6j + 4 by the double staircase.
- FINITE-EXACT for n <= 12: 1, 1, 1, 1, 1, 1, 3, 2, 2, 3 palindromic tight colourings for n = 3..12, each realised by
  a blue tight line. At n = 12 the only blue pair besides (0, 4) is (2, 2) (296 lines). The palindromic colourings
  with P(A) = 1 are not tight there.

**Falsifier.** A palindromic tight colouring with no blue tight line.

## (D) The (1,3) excess

**Statement.** At odd n >= 11 the blue lines with (H mod 4, H mod 4) = (1, 3) are more than half.

**Data.**
- The share is exactly 1/2 at n = 5, 7, 9 and 0.52148 at n = 11.
- By THM-4951 (5d) the excess is 2^(c-2)(|D_1| - |D_0|), where D is the set of y with B(y) = B(y + 1).

**Test.** n = 13 has 2^30 values of y. A sampled estimate needs the mod-2 transfer DP of THM-4922 (4), step 5, per
sample.

**Falsifier.** A sampled or exact (1,3)-share <= 1/2 at n = 13.
