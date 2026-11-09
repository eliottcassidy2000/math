# Pre-registration 4: the blue-line (1,3) excess at orders 13 and 15 (HYP-9540 (D))

opus-2026-11-09-S87. Written and pushed before `mirror_sample.c` was written or run.

## Objects

Notation as in THM-4951 (5): n = 2c + 1. The blue lines are classified by the pair (H(t) mod 4, H(t + 1) mod 4).

- By THM-4922 (4), (H - 1)/2 = A(y) + x.B(y) (mod 2), where y are the coordinates on the 2-element sigma-orbits of
  tiles and x the c mirror tiles. Hence:
  - D = {y : B(y) = B(y + 1)};
  - D_e = {y in D : A(y) + A(y + 1) + 1 = e};
  - the (1,3)-share of blue lines is 1/2 + (|D_1| - |D_0|) / 2^(N + 1), with N the number of y-coordinates.
- Exact values so far:
  - D is empty at n = 5, 7, 9, so the share is 1/2.
  - At n = 11, |D| / 2^20 = 0.08984 and (|D_1| - |D_0|) / 2^20 = 0.04297, so the share is 0.52148.

## Method (to be written)

- `mirror_sample.c` computes the mirror-free path parities N(u) from c by a DP over the sets of used mirror pairs.
- It is validated exhaustively at n = 9, where D is empty, and at n = 11, where |D| = 94,208 and
  |D_1| - |D_0| = 45,056 must be reproduced exactly.
- It then samples uniformly random y: at least 10^7 at n = 13 and at least 10^6 at n = 15.

## Predictions

- **U1 (n = 13).** The estimated (1,3)-share exceeds 1/2 by more than 5 standard errors.
- **U2 (n = 13).** The estimated density |D| / 2^N is > 0.01. The degeneracy does not vanish.
- **U3 (n = 15).** The estimated (1,3)-share exceeds 1/2 by more than 5 standard errors.

Outcome: to be appended below without editing the predictions.

## Outcome (appended after the runs; the predictions above are unedited)

`mirror_sample.c` reproduces the exact values exhaustively: D is empty at n = 5, 7, 9; at n = 11, |D| = 94,208,
|D_1| = 69,632 and |D_0| = 24,576. A 100,000-sample calibration run at n = 13 was made before the full runs. The full
runs had 12 slices with fixed seeds (mirror_sample.out).

- **U1 REFUTED.** At n = 13 (24,000,000 samples) the (1,3)-share is 1/2 - 0.000906, with standard error 0.000017
  (z = -53.6). The excess changes sign: (1,3) lines are slightly FEWER than half.
- **U2 CONFIRMED.** At n = 13 the estimated density |D| / 2^30 is 0.02744.
- **U3 CONFIRMED.** At n = 15 (12,000,000 samples) the (1,3)-share is 1/2 + 0.001689, with standard error 0.000022
  (z = +75.8).

## Addendum (written and pushed after the outcome above, before any n = 17 run)

The signs so far are + at n = 11, - at n = 13 and + at n = 15. This suggests sign = (-1)^(c+1) with n = 2c + 1.

- **U4 (n = 17, c = 8).** The (1,3)-share is below 1/2 by more than 5 standard errors. Test: at least 2,400,000
  uniform samples.

Outcome: to be appended below without editing U4.
