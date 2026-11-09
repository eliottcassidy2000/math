---
id: HYP-9540
title: "Packing-number line laws of the tiling metagraph: (A) at every order the triple words of tight lines (nu(t) + nu(t+1) = floor(n/3)) are exactly the tight colourings - proved for n != 0 (mod 3), open at n = 3k; (B) black lines join nu-classes i <= j iff i >= 1 and i + j >= floor(n/3); blue lines join them iff i >= 1 and i + j > floor(n/3), or the pair is that of a palindromic tight colouring (including the transitive line (0, floor(n/3))); (C) the triple words of BLUE tight lines are exactly the palindromic tight colourings - proved for n != 0 (mod 3); (D) the blue-line (1,3) excess of THM-4951 (5) is nonzero at every odd n = 2c + 1 >= 11 with sign (-1)^(c+1) (the earlier clause 'the excess persists' is REFUTED at n = 13)."
status: >
  OPEN.
  (A): PROVED for n != 0 (mod 3) (THM-4951 (3a), (3e), (3f): block laws and the run-window theorem). At n = 3k it is
  FINITE-EXACT for n <= 12 (50 of 50 words at n = 12, pre-registered T1).
  (B): FINITE-EXACT for n <= 9 (complete census of all lines). The blue clause is also FINITE-EXACT for n <= 12 (the
  audit's census of all blue lines).
  - The boundary i + j = k is realised in black at every n (THM-4951 (2a)). In blue exactly the pairs of palindromic
    tight colourings are possible (THM-4951 (3d)). Through (C), at n != 0 (mod 3) exactly those occur.
  - (0, k) and (1, k) occur at every n >= 4: the transitive line, and the single reversed tile (0, 3), which is black
    for n >= 5.
  - The first version of the blue clause ("the same pairs as black, except (1, 1) from n = 6 on") is REFUTED. At
    n = 12 no blue line has nu-pair (1, 3): such a line would be tight, with a palindromic tight word, and the three
    palindromic tight colourings of [0, 9] have P-pairs (2, 2), (2, 2) and (4, 0). The audit found this (MISTAKE-742).
  (C): PROVED for n != 0 (mod 3) (THM-4951 (3d), (3f)); FINITE-EXACT for n <= 12.
  (D): the original clause (the (1,3) excess persists at every odd n >= 11) is REFUTED at n = 13: pre-registered U1,
  and the audit's independent sampler. The sign law replacing it was stated after U1-U3 and passed the pre-registered
  U4 at n = 17. Its data are n = 11 exact and n = 13, 15, 17 sampled. OPEN.
  Audited with THM-4950 and THM-4951: PASS WITH CORRECTIONS (M1, M2 applied; MISTAKE-742).
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

**What is known.**
- One direction is proved at every n (THM-4951 (1)).
- The converse is proved at n != 0 (mod 3) (THM-4951 (3a), (3e), (3f)):
  - the tight colourings there are block words;
  - all their runs have length >= 2;
  - so the run-window tilings realise them.
- At n = 3k the converse is FINITE-EXACT for n <= 12.

**The open case n = 3k.**
- The tight colourings may have runs of length 1. There are (k^2 + 9k - 2) 2^(k-4) of them up to complement (fitted).
- The run-window tiling realises 5 of 5, 16 of 17, 47 of 50, 128 of 136 and 332 of 352 at n = 6, 9, 12, 15, 18.
- For n = 3k <= 18 it fails exactly when an isolated 0 lies between two 1-runs of length >= 2: every failing word
  contains 1110111. This is the follow-up audit's FINITE-EXACT characterisation. The first statement here ("an
  isolated 0 between two 1-runs") was wrong: 1010 at n = 6 is realised (MISTAKE-742).
- The failing words are still realised by other tight lines at n = 9 and 12. So a proof at n = 3k needs a
  construction for the 1110111 pattern.

**Falsifier.** A tight colouring at some n = 3k that no tight line realises.

## (B) Allowed pairs

**Statement.** For n >= 4 and nu-classes i <= j:
- some BLACK line joins them iff i >= 1 and i + j >= k;
- some BLUE line joins them iff one of:
  - i >= 1 and i + j > k;
  - i + j = k and {i, j} = {P(A), P(B)} for a palindromic tight colouring. This includes (0, k), from the all-ones
    word.

**Notes.**
- The blue boundary pairs are fixed by (C). For n <= 26 they are exactly the {j, k - j} with j even or k odd
  (audit, palindromic_pairs.out). At n = 3k + 2 this is proved (THM-4951 (3d)).
- The pair (1, 1) is excluded in blue by this rule from n = 6 on, so no separate exception is needed. At n = 6, 7, 8 it
  is a tight pair with k = 2 and i = 1 odd. From n = 9 on, i + j = 2 < k.

**What is known.**
- Necessity is proved: Session26's bound, the uniqueness of nu = 0, and the palindrome condition on blue tight words.
- Sufficiency on the boundary is proved: in black at every n (THM-4951 (2a)), in blue at n != 0 (mod 3) (THM-4951
  (3f)).
- (1, k) is proved at every n >= 4. Reverse the single tile (0, 3): every cycle uses 3 -> 0, so nu(t) = 1, while
  t + 1 keeps all consecutive 3-cycles on {3b, 3b + 1, 3b + 2}.
- The complete census n <= 9 (all lines) and the audit's census of all blue lines for n <= 12 have exactly these
  pairs. At n = 10, 11, 12 every non-tight allowed pair occurs in blue.

**Falsifier.** An allowed pair with no line of the stated colour at some n >= 10, or a forbidden pair that occurs.

## (C) Blue word law

**Statement.** For every n >= 3, a word is the triple word of a BLUE tight line iff it is a palindromic tight
colouring.

**What is known.**
- Necessity is proved: sigma reverses the triple word.
- Sufficiency is proved at n != 0 (mod 3): the run-window tiling of a palindrome is blue (THM-4951 (3f)).
- At n = 3k + 2 the blue tight pairs are exactly {j, k - j} with j even or k odd, so for example no blue (1, 3) tight
  line exists at n = 14.
- The double staircase with parameter j realises the pair {2j, k - 2j} in blue at every n >= 6j + 4.
- FINITE-EXACT for n <= 12: there are 1, 1, 1, 1, 1, 1, 3, 2, 2, 3 palindromic tight colourings for n = 3..12, each
  realised by a blue tight line. At n = 12 the only blue pair besides (0, 4) is (2, 2) (296 lines).

**Falsifier.** A palindromic tight colouring at some n = 3k with no blue tight line.

## (D) Sign law of the (1,3) excess

**Statement.** For odd n = 2c + 1 >= 11, the (1,3)-share of blue lines is 1/2 + s_n, with s_n != 0 and
sign(s_n) = (-1)^(c+1). So s_n > 0 at n = 3 (mod 4) and s_n < 0 at n = 1 (mod 4).

**Data.**

| n | 11 (exact) | 13 | 15 | 17 |
|---|---|---|---|---|
| s_n | +0.021484 | -0.000906 (z = -54) | +0.001689 (z = +76) | -0.000111 (z = -8) |
| D-density | 0.0898 | 0.0274 | 0.0238 | 0.0092 |

- s_n = 0 exactly at n = 5, 7, 9, where D is empty.
- By THM-4951 (5d), s_n = (|D_1| - |D_0|) / 2^(N + 1), with D = {y : B(y) = B(y + 1)}.

**History.** The first version of (D), "the excess persists at every odd n >= 11", was refuted at n = 13 by the
pre-registered U1. The audit's independent sampler confirmed it. The sign law was formulated after U1-U3 and confirmed
by the pre-registered U4 at n = 17.

**Falsifier.** A sign violation at n = 19 or 21, or s_n = 0 exactly at some odd n >= 11.
