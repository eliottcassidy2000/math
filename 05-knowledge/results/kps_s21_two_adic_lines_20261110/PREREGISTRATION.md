# Pre-registration: unit-offset seeds at order 9 (kind-pasteur-2026-11-10-S21)

Written and pushed BEFORE the order-9 run. Chan-Pak (every odd number other than 7 and 21 is a Hamiltonian-path
count of a tournament) is OPEN and stays open whatever the outcome.

## Objects

For a rooted tournament (L, v): H = H(L), d = H(L - v), E = number of Hamiltonian paths of L ending at v,
U = number of covers of L by an ordered pair of paths (P1, P2) with P1 possibly empty and P2 starting at v.
The known fixed-root tail law (O/Z tails: Sessions 195-205, THM-4888, THM-5085 at x = 2) gives, when t new
vertices are appended that beat v and lose to everything else,

    H(L o 0^t) = c0 + M 2^t   (t >= 1),      c0 = H - E + 2d - U,   M = U - d.

A *unit-offset root* has c0 = 1. It puts the whole ray {1 + M 2^t : t >= 1} into the spectrum. Write M = c 2^j with
c odd: the ray is the tail s >= j + 1 of the 2-adic line L_c = {1 + c 2^s : s >= 1}.

## In-sample (orders <= 8, program scripts/seeds.c on the class lists of classgen.c)

Unit-offset roots: 2, 5, 7, 16, 31, 68, 207 at orders 2..8. Odd c with a seed at order <= m: 1, 1, 2, 3, 7, 24, 89
for m = 2..8. Among odd c <= 63 the ten values 9, 13, 19, 27, 33, 35, 41, 53, 59, 63 have no tail seed at order <= 8
(9 has a module seed at order 6).

## Predictions for the order-9 run (191,536 classes, every root)

- **P1.** The number of distinct odd c with a unit-offset tail seed of order <= 9 lies in [250, 420].
- **P2.** At least 6 of the ten odd c <= 63 listed above acquire a tail seed at order 9.
- **P3.** c = 13 acquires a tail seed at order 9.
- **P4.** c = 33 (the weight-three row d = 5: 1 + 2^b + 2^(b+5)) acquires a tail seed at order 9.
- **G1 (gate, forced by THM-1370, not a prediction).** No ray of the run contains 7 or 21: the least s0 for c = 3
  stays 2 and the least s0 for c = 5 stays >= 3. A violation means the program is wrong.

Refutation criteria: each of P1-P4 fails exactly when its sentence is false for the run. Outcomes will be appended
below without editing the text above.

## Outcomes (appended after the run; the text above is unchanged)

Run: `classgen 9` (191,536 classes), `seeds 9` on every class and root; 1,333 unit-offset roots at order nine.

- **P1 FAILED.** 636 odd c have a tail seed of order <= 9 (predicted 250-420; the count grew by 7.1, not 3.5).
- **P2 CONFIRMED.** 8 of the ten: 9, 13, 27, 35, 41, 53, 59, 63. Still without a tail seed: 19, 33.
- **P3 CONFIRMED.** c = 13: order nine, (H, d, E, U) = (183, 61, 35, 269), M = 208 = 13 * 16.
- **P4 FAILED.** c = 33 has no unit-offset root of order <= 9.
- **G1 holds.** Least s0: 2 for c = 3, 3 for c = 5.

Odd c <= 129 without a seed (tail of order <= 9 or module of order <= 8): 19, 33, 71, 79, 83, 97, 101, 103, 107,
113, 115.

# Second pre-registration: order ten through the cut matrix (written and pushed BEFORE the run)

Engine `scripts/unitsets.c`: for a base B, a vertex set O is the out-neighbourhood of a unit-offset root of
L = B + v iff the sum of the two-path-cover couplings W_ef over the pairs of O equals H(B) - 1; the multiplier is
the sum of the responses U_B(y), y in O. Validated before this run: on all bases of orders 7 and 8 the sets of
(out-degree, M) agree with the direct runs of `seeds.c` at orders 8 and 9 (140 and 1,007 pairs).

Run: every class B of order 9 (191,536), every O with 2 <= |O| <= 4. This finds the unit-offset roots of order ten
with root out-degree 2, 3 or 4; out-degrees 0 and 1 give only fans; out-degree >= 5 is NOT searched.

- **Q1.** At least 6 of the eleven odd c <= 129 without a seed (19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115)
  acquire a seed.
- **Q2.** c = 19 acquires a seed.
- **Q3.** c = 33 acquires a seed.
- **Q4.** The number of odd c <= 1025 with a seed rises from 377 to at least 450 (of 513).
- **G2 (gate, forced by THM-1370).** No ray contains 7 or 21.
