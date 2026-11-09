# Pre-registration 2: the blue-line H-mod-4 law at order 11

opus-2026-11-09-S87. Written and pushed before `mirror_check.c` was run at n = 11 and before any order-11 blue line was
evaluated.

## Background (in-sample)

- **The tiling frame.** In the tiling line metagraph (THM-643 frame), a blue line is a complement pair {t, t + 1} of
  sigma-fixed tilings.
- **The law.** At odd n the blue lines were observed to split exactly 1 : 2 : 1 by the pair (H(t) mod 4,
  H(t + 1) mod 4):
  - n = 5: 2 : 4 : 2;
  - n = 7: 64 : 128 : 64;
  - n = 9: 8192 : 16384 : 8192.
  Order 9 was computed in the same exploratory batch, so it is not a test.
- **Even n.** At even n every blue line is (1, 1), which is known from THM-791.
- **The mirror-path decomposition.** Write n = 2c + 1, let iota(i) = n - 1 - i, and let N(u) be the number of
  mirror-free paths from c to u. Then, with x the mirror tiles and y the other blue coordinates:

      (H(t) - 1)/2 = A(y) + sum_i x_i B_i(y)   (mod 2),   with B_i(y) = N(i) + N(n - 1 - i).

  This is proved by the symmetric-cycle argument and checked on every blue tiling for n <= 7.
- **The two parity facts.** These held for every y at n = 5, 7, 9:
  - (P2) sum_i B_i(y) = 1 (mod 2);
  - (P3) B(y) != B(y + 1).
- **Why they suffice.** Together they give the 1 : 2 : 1 law, because x -> x.B(y) and x -> x.(B(y) + B(y + 1)) are
  then nonzero linear forms.

## Predictions for n = 11

There are 2^20 values of y and 2^24 blue lines.

- **R1.** (P2) holds for every y.
- **R2.** (P3) holds for every y.
- **R3.** Hence the blue lines at n = 11 split exactly 4,194,304 : 8,388,608 : 4,194,304 by (H mod 4) pair.
  - R3 is implied by R1, R2 and the proved decomposition.
  - As an independent check of the decomposition at n = 11: on 20,000 random blue tilings, the H mod 4 computed
    directly by Hamiltonian-path DP equals A(y) + x.B(y).

Outcome: to be appended below without editing the predictions.
