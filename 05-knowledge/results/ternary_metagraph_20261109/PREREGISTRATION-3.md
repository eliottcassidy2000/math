# Pre-registration 3: tight lines of the packing-number classification at orders 12-14

opus-2026-11-09-S87. Written and pushed while `tight_lines.c` was running at n = 12, before any of its output was
read. The n = 13 and n = 14 realizability searches have not been written yet.

## Objects

- **Packing number.** nu(T) is the maximum number of vertex-disjoint cycles, which equals the maximum number of
  vertex-disjoint 3-cycles.
- **The line bound.** On every line of the tiling metagraph (THM-643 frame), nu(t) + nu(t + 1) >= k = floor(n/3). The
  mechanism is Session26's consecutive-triples argument, which that session stated for the feedback vertex number.
- **Tight lines.** A line is tight if nu(t) + nu(t + 1) = k.
- **Triple word.** The triple word of a tiling is x_a = [tile (a, a + 2) reversed] for a = 0, ..., n - 3. A line is
  recorded by the representative with x_0 = 1.
- **Tight colouring.** A colouring of the positions 0..L-1 (L = n - 2) into A and its complement B is tight if
  P(A) + P(B) = ceil(L/3) = k. Here P(X) is the largest subset of X with pairwise differences >= 3.
- **Necessity (proved).** The triple word of every tight line is a tight colouring, because consecutive 3-cycles on
  3-separated triples are disjoint.
- **Small orders (proved).** At n = 3k + 2 the tight colourings are exactly the block-monochromatic words on the
  aligned blocks {3b, 3b + 1, 3b + 2}.

## In-sample data (n <= 11, computed by `tight_lines.c`, which reproduces the full n <= 9 census)

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|
| tight lines | 1 | 1 | 1 | 57 | 31 | 13 | 3673 | 1802 | 541 |
| blue | 1 | 1 | 1 | 1 | 1 | 1 | 39 | 10 | 19 |
| tight colourings (x_0 = 1) | 1 | 1 | 1 | 5 | 3 | 2 | 17 | 8 | 4 |

- For n = 6..11 the set of triple words of tight lines equals the set of tight colourings.
- For k >= 2 the only nu-pairs on tight lines have been (0, k), the transitive line, and (1, k - 1).

## Predictions

- **T1 (word law at n = 12).** The triple words of the tight lines at n = 12 are exactly the 50 tight colourings of
  [0, 9] with x_0 = 1.
- **T2 (pairs at n = 12, k = 4).** The tight lines at n = 12 realise the nu-pairs (0, 4), (1, 3) and (2, 2), each with
  at least one line.
- **T3 (count at n = 12).** The number of tight lines at n = 12 lies in [10^5, 10^6].
  - Point estimate: 2.3 * 10^5, from the in-class growth ratios 57-64 for n = 0 (mod 3).
- **T4 (blue at n = 12).** At least one blue tight line exists at n = 12 besides the transitive line.
- **T5 (realizability at n = 13, 14).** Each tight colouring is realised by at least one tight line: 20 colourings at
  n = 13, and the 8 block-monochromatic words at n = 14. This will be tested by a first-hit search, one per
  colouring, with the triple word fixed.

Outcome: to be appended below without editing the predictions.

## Outcome (appended after the runs; the predictions above are unedited)

The n = 12 count used `tight_lines.c` in 12 slices split at vertex 8. The first-hit searches at n = 13, 14 used
`tight_lines.c` in first-hit mode. The data are in tight_lines.out and t5_first_hits.out.

- **T1 CONFIRMED.** The tight lines at n = 12 have exactly 50 distinct triple words, and they are the 50 tight
  colourings of [0, 9].
- **T2 CONFIRMED.** The nu-pairs at n = 12 are (0, 4): 1 (blue, transitive line), (1, 3): 96,254 (all black) and
  (2, 2): 467,058 black + 296 blue.
- **T3 CONFIRMED.** The count is 563,609 tight lines, inside [10^5, 10^6]. The point estimate 2.3 * 10^5 was low by a
  factor 2.4: the n = 0 (mod 3) class grows by 153x from n = 9 to n = 12.
- **T4 CONFIRMED.** There are 297 blue tight lines: the transitive line and 296 of pair (2, 2). There are no blue
  (1, 3) tight lines at n = 12.
- **T5 CONFIRMED.** Every tight colouring is realised: 20 of 20 at n = 13 and 8 of 8 at n = 14. At n = 14 the
  palindromic word 111000000111 was first hit by a blue line.

During the runs, the n = 3k + 2 half of the word law was proved for every k (THM-4951 (3); multi-staircase). The
n = 14 part of T5 is therefore also a theorem. The n = 12 and n = 13 parts remain finite checks.
