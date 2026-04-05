# THM-308: General Base QR Congruence

**Status:** PROVED
**Found by:** opus-2026-04-05-S26
**Verified in:** `04-computation/binary_resonance_cross_field_s26.py`

## Statement

For any integer b >= 2 and any odd prime p with gcd(b, p) = 1:

**b^{C(p,2)} ≡ (b/p) (mod p²)**

where (b/p) is the Legendre symbol.

This generalizes THM-306 (which was the b = 2 case) to ALL bases.

## Proof

Identical to THM-306: b^{C(p,2)} = b^{p(p-1)/2} = (b^{(p-1)/2})^p. By Euler's criterion, b^{(p-1)/2} = (b/p) + pk for some integer k. Since (b/p) = ±1 and p is odd, ((b/p) + pk)^p ≡ (b/p) (mod p²). QED.

## Classification of Bases

**Perfect square bases** (b = 1, 4, 9, 16, 25, ...): (b/p) = 1 for all p. These always satisfy b^{C(p,2)} ≡ 1 (mod p²). **Trivial QR resonance.**

**Non-square bases**: (b/p) varies with p, creating **nontrivial QR resonance**.

| Base b | Discriminant | Resonance pattern |
|--------|-------------|-------------------|
| 2 | p mod 8 | (2/p) = 1 iff p ≡ ±1 mod 8 |
| 3 | p mod 12 | (3/p) = 1 iff p ≡ ±1 mod 12 |
| 5 | p mod 5 | (5/p) = 1 iff p ≡ ±1 mod 5 |
| 6 = 2×3 | p mod 24 | (6/p) = (2/p)(3/p) |
| 7 | p mod 28 | (7/p) governed by QR of 7 |

## Combinatorial Applications

| Structure | Base | Labeled count mod p² |
|-----------|------|---------------------|
| Tournaments | 2 | ≡ (2/p) |
| Graphs | 2 | ≡ (2/p) |
| Digraphs | 4 = 2² | ≡ 1 always (trivial!) |
| Oriented graphs | 3 | ≡ (3/p) |
| k-edge-colorings | k | ≡ (k/p) |
| q-state Potts model | q | ≡ (q/p) |

## Corollary: Digraph Triviality

The labeled digraph count 4^{C(p,2)} ≡ 1 (mod p²) for ALL odd primes p. Digraphs have NO QR discrimination at the labeled level because 4 is a perfect square.

## Verified Cases

All 70 combinations of b ∈ {2,...,10} and p ∈ {3,5,7,11,13,17,19,23} verified. See computation script for full table.

## Related

- THM-305: v_2(T(n)) = (n-1)/2 (the tournament-specific consequence)
- THM-306: The b = 2 case
- THM-307: T(p) mod p decomposition
