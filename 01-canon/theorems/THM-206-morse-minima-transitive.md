# THM-206: Morse Minima of H are Exactly the Transitive Tournaments

**Status:** VERIFIED (exhaustive n <= 6), forward direction PROVED
**Source:** opus-2026-03-14-S89

## Statement

A tournament T on n vertices is a local minimum of H on the tournament hypercube Q_{C(n,2)} if and only if T is transitive (H(T) = 1).

Equivalently: H(T) <= H(T^e) for all edges e if and only if H(T) = 1.

## Proof

**(=>)** If T is transitive, H(T) = 1. Since H(T') >= 1 for all tournaments T' (Redei), we have H(T) <= H(T') for ALL T'. In particular H(T) <= H(T^e) for every edge e. So T is a global minimum, hence a local minimum. ∎

**(<=)** Verified exhaustively for n = 3, 4, 5, 6. For every non-transitive tournament T at these n, there exists an edge e with H(T^e) < H(T).

**Reverse direction evidence (very strong):** Exhaustive verification reveals a much stronger property:

**CONJECTURE (Upset Monotonicity):** If i->j is an arc with score(i) < score(j), then H(T) > H(T^{ij}). That is, flipping any "upset" arc always decreases H.

This was verified for ALL tournaments at n=3,4,5,6 (100% success rate, zero exceptions). It immediately implies the reverse direction: every non-transitive tournament has at least one upset arc (since the score sequence is not {0,1,...,n-1}), and flipping it decreases H.

Additional verified properties (n=3..6):
- Every 3-cycle has at least one arc whose flip decreases H
- In each 3-cycle, flipping the arc where the highest-score vertex loses always decreases H

## Verification

- n=3: 6 local minima = 6 = 3! transitive tournaments
- n=4: 24 local minima = 24 = 4! transitive tournaments
- n=5: 120 local minima = 120 = 5! transitive tournaments
- n=6: 720 local minima = 720 = 6! transitive tournaments

## Corollary

The number of local minima of H on Q_{C(n,2)} is exactly n!.

## Additional Findings

### Local Maxima
- n=3: 2 local maxima (both 3-cycles, H=3, regular score (1,1,1))
- n=4: 24 local maxima (all H=5, near-regular score (1,1,2,2))
- n=5: 64 local maxima (H=15, scores (2,2,2,2,2) and (1,2,2,2,3))
- n=6: 1200 local maxima — first n with non-global local maxima:
  - 480 at H=45 (global max, score (2,2,2,3,3,3))
  - 720 at H=37 (local but not global, score (1,2,2,3,3,4))

### Gradient Flow
Under steepest ascent, all tournaments flow to the global maximum for n=3,4,5. At n=6, 88.3% flow to H=45, 11.7% flow to H=37.

### Edge Sensitivity
|dH| is always even when flipping one arc (trivial from Redei: odd - odd = even). The Walsh derivative formula: dH_e(T) = 2 * sum_{S: e in S} H_hat[S] * chi_S(T).
