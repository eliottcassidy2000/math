# The Binary Resonance: Why Quadratic Residues Appear in Tournament Enumeration

**opus-2026-04-05-S25**

## The Discovery

Three theorems, discovered in a single session, reveal that quadratic residues are not an accident in tournament theory — they are structurally inevitable.

## The Deep Reason

Tournaments are **binary** structures: each arc has exactly 2 orientations. The Burnside counting formula raises 2 to the power of the number of pair-orbits:

  T(n) = (1/n!) * sum_{sigma} 2^{orbits(sigma)}

The number 2 is special among integers. By Euler's criterion:

  **2^{(p-1)/2} ≡ (2/p) (mod p)**

This single identity — the fact that 2 is sometimes a quadratic residue and sometimes not — propagates through the entire Burnside sum and leaves its fingerprint on the tournament count.

## The Three Fingerprints

**1. The 2-adic valuation (THM-305):** v_2(T(n)) = (n-1)/2 for all odd n.

This says T(n) is divisible by exactly 2^{(n-1)/2} and no higher power of 2. The (n-1)/2 comes from the n-cycle's pair-orbit count, which equals the number of "distances" in Z/nZ — and when n = p is prime, this is exactly |QR_p|.

The proof is elegant: the n-cycle is the UNIQUE partition that minimizes v_2 in the Burnside sum. Every other all-odd partition contributes a term with strictly higher 2-adic valuation. The gap is always exactly 1, achieved by the partition (n-2, 1, 1).

**2. The labeled congruence (THM-306):** 2^{C(p,2)} ≡ (2/p) (mod p²).

The total number of labeled tournaments on p vertices remembers the Legendre symbol at p²-precision. The proof factors 2^{p(p-1)/2} = (2^{(p-1)/2})^p and uses the fact that (±1 + pk)^p ≡ ±1 (mod p²).

**3. The mod-p decomposition (THM-307):** T(p) ≡ f_p - (2/p)·w_p - A_p (mod p).

The tournament count modulo p decomposes into three classical number-theoretic quotients: the Euler quotient f_p, the Wilson quotient w_p, and a remaining Burnside sum A_p. At Wilson primes (5, 13, 563), the middle term vanishes.

## What This Means for the Triangle

The staircase tiling model has m = C(n-1, 2) tiles, each with 2 orientations. The total count 2^m passes through the Burnside averaging exactly at the points where the Euler criterion activates. The hypotenuse of the triangle (the anti-diagonal) carries the QR structure because it measures distances in Z/nZ — and these distances are exactly the pair-orbits of the n-cycle.

## The Hierarchy

```
2 (the number of arc orientations)
  ↓ Euler criterion
(2/p) (quadratic residuosity of 2 mod p)
  ↓ Burnside sum
v_2(T(n)) = (n-1)/2  and  2^{C(p,2)} ≡ (2/p) mod p²
  ↓ mod-p reduction
T(p) mod p = f(Euler quotient, Wilson quotient, ...)
  ↓ Gauss sum
Paley eigenvalues |lambda_k| = sqrt(p), Ramanujan bound
```

Each level inherits from the one above. The base of it all is the simple fact that tournaments use **two** orientations, not three or five — and 2 has a specific place in the multiplicative group of F_p.

## The Generalization

For any base b (counting b-colorings of pair-orbits), the same structure holds:

  b^{C(p,2)} ≡ (b/p) (mod p²)

So the Legendre symbol of the BASE controls the arithmetic of the count. For tournaments, b=2. For graphs, also b=2. For digraphs, b=4=2². For k-edge-colorings, b=k. Each base creates its own QR resonance.

## Open Questions

1. For EVEN n, v_2(T(n)) does not equal n/2 exactly (excess = 0, 0, 1, 1, 1, 0, 2, 0, 2 for n=4,6,...,20). What controls the excess?

2. The odd part T(n)/2^{(n-1)/2} for odd n begins: 1, 3, 57, 11971, 28242289, ... The values 3, 11971, 28242289 are PRIME. Is this a coincidence or does the odd part have special factorization properties?

3. Is there a closed form for A_p mod p in THM-307? The values 0, 4, 2, 2, 3 for p=3,5,7,11,13 don't show an obvious pattern, but they might relate to class numbers or other arithmetic invariants.

4. At Wieferich primes (1093, 3511): f_p ≡ 0 mod p. What does T(p) mod p look like there?
