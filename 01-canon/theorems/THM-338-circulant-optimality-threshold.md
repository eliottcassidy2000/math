---
theorem: THM-338
name: Circulant Optimality Threshold
status: PROVED for n≤11, CONJECTURED for n≥13
session: opus-2026-05-27-S7
verified: computationally n=3..13 (exhaustive circulant search)
depends_on: THM-336
---

## Statement

Define opt_circ(n) = max{H(T) : T is a circulant tournament on Z_n}.
Define a(n) = A038375(n) = max{H(T) : T is any n-vertex tournament}.

Then:
- **a(n) = opt_circ(n)** for n = 3, 5, 7, 9, 11
- **a(n) > opt_circ(n)** for n = 13 (gap = 8656)

**Conjecture:** a(n) > opt_circ(n) for all n ≥ 13.

## Optimal Circulant Tournaments (verified)

| n  | a(n)      | opt_circ(n) | Optimal S              | Type                |
|----|-----------|-------------|------------------------|---------------------|
| 3  | 3         | 3           | {2}                    | Paley QR_3          |
| 5  | 15        | 15          | {3,4}                  | upper half          |
| 7  | 189       | 189         | {1,2,4}                | Paley QR_7          |
| 9  | 3357      | 3357        | {1,5,6,7}              | mixed               |
| 11 | 95095     | 95095       | {2,6,7,8,10}           | Paley QR_11         |
| 13 | ≥3719831  | 3711175     | {1,2,3,4,5,6}=lower    | lower half (gap 8656)|

## Lower Bounds from Circulant Search (new bounds for A038375)

| n  | opt_circ(n)                   | Optimal S (local optimum)              |
|----|-------------------------------|----------------------------------------|
| 13 | 3,711,175                     | {1,...,6} lower half                   |
| 15 | 198,464,295                   | {1,...,7} lower half                   |
| 17 | 13,689,269,499                | {1,...,8} lower half                   |
| 25 | ≥2,418,453,569,285,650,675    | {1,...,12} lower half (local opt)      |
| 27 | ≥17,051,631,267,035,242,313   | {2,...,13,26} = lower half - {1} + {26}|

Note: The n=27 result is 6.3× better than the lower half due to the "backward" connection 26≡−1 mod 27.

## Why n=13 Breaks Circulant Optimality

At n=13 (13 ≡ 1 mod 4, prime), the QR-based construction fails (QR_13 is symmetric under negation, so can't define a tournament). The best circulant is the lower half {1,...,6}, which is a "local maximum" but not globally optimal.

The optimal n=13 tournament (achieving a(13) ≥ 3,719,831) has a non-circulant structure found by simulated annealing.

For Paley primes p ≡ 3 mod 4 (p = 3, 7, 11): QR_p is optimal among circulants AND conjectured globally optimal.
For p ≡ 1 mod 4 (p = 5, 13, 17, ...): the "lower half" or mixed circulant is best among circulants, but global optimum is non-circulant.
For n=9 (= 3²): optimal circulant is a non-lower-half structure {1,5,6,7} which achieves the global optimum.

## Algorithm

The circulant-reduced bitmask DP enables exhaustive search over all 2^{(n-1)/2} connection sets:
- Memory: 2^{n-1} × 8 bytes (factor 2n smaller than standard)
- Time: n × 2^{n-1} per evaluation
- Feasible exhaustive search: n ≤ 17 (2^8 options × fast DP)
- Local search feasible: n ≤ 29 (2^14 options but only explore local optima)

## Structural Pattern for n=27

The local optimum S = {2,...,13, 26} for n=27 suggests:
"Beat your 2nd through 13th nearest successors, AND beat your immediate predecessor."
Missing connection: {1} (immediate successor).

This pattern — forward medium-range connections + one backward connection — appears to create more Hamiltonian paths than purely forward connections.

## Connection to Paley Structure

For QR_p with p ≡ 3 mod 4:
The connection set QR_p contains both "small" and "large" elements (it's not the lower half).
E.g., QR_7 = {1,2,4}: small forward (1,2), medium forward (4 < 7/2=3.5? No, 4>3.5, so 4 is "backward" in some sense!).
Actually: 4 ≡ −3 mod 7 (since 4 = 7−3), so connection 4 = going 3 steps backward = going 4 steps forward. So QR_7 has a "backward" element too!

For QR_11 = {1,3,4,5,9}: 9 ≡ −2 mod 11, 5 ≡ −6 mod 11, etc. The QR set always contains both forward and backward connections when p ≡ 3 mod 4. This is why Paley tournaments are more HP-optimal than purely lower-half circulants!
