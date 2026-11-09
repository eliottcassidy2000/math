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
