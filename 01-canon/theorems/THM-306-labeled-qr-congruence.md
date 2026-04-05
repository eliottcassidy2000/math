# THM-306: Labeled QR Congruence

**Status:** PROVED
**Found by:** opus-2026-04-05-S25
**Verified in:** `04-computation/qr_tournament_enumeration_s25.py`

## Statement

For every odd prime p:

**2^{C(p,2)} ≡ (2/p) (mod p²)**

where (2/p) is the Legendre symbol and C(p,2) = p(p-1)/2.

In words: the total number of labeled tournaments on p vertices is congruent to the Legendre symbol of 2 modulo p-squared.

## Proof

Write C(p,2) = p(p-1)/2. Then:

2^{p(p-1)/2} = (2^{(p-1)/2})^p

By the Euler criterion: 2^{(p-1)/2} ≡ (2/p) (mod p).

So 2^{(p-1)/2} = (2/p) + p·k for some integer k.

**Case (2/p) = 1:** 2^{(p-1)/2} = 1 + pk. Then:
  (1 + pk)^p = 1 + C(p,1)·pk + C(p,2)·(pk)² + ...
  ≡ 1 + p²k (mod p²)  [since C(p,1)·pk = p²k, higher terms have p³|...]
  ≡ 1 = (2/p) (mod p²). CHECK.

**Case (2/p) = -1:** 2^{(p-1)/2} = -1 + pk. Then:
  (-1 + pk)^p = (-1)^p + C(p,1)·(-1)^{p-1}·pk + ...
  = -1 + p²k + O(p³)   [since p is odd, (-1)^p = -1, (-1)^{p-1} = 1]
  ≡ -1 = (2/p) (mod p²). CHECK.

In both cases: (2^{(p-1)/2})^p ≡ (2/p) (mod p²). QED.

## Interpretation

The Legendre symbol (2/p) equals:
  +1 if p ≡ ±1 (mod 8)
  -1 if p ≡ ±3 (mod 8)

So the labeled tournament count 2^{p(p-1)/2} satisfies:
  ≡ 1 (mod p²)   when p ≡ ±1 (mod 8)
  ≡ -1 (mod p²)  when p ≡ ±3 (mod 8)

## Verified Cases

| p | C(p,2) | (2/p) | p² | 2^{C(p,2)} mod p² |
|---|--------|-------|----|--------------------|
| 3 | 3 | -1 | 9 | 8 = -1 mod 9 |
| 5 | 10 | -1 | 25 | 24 = -1 mod 25 |
| 7 | 21 | 1 | 49 | 1 mod 49 |
| 11 | 55 | -1 | 121 | 120 = -1 mod 121 |
| 13 | 78 | -1 | 169 | 168 = -1 mod 169 |
| 17 | 136 | 1 | 289 | 1 mod 289 |
| 19 | 171 | -1 | 361 | 360 = -1 mod 361 |
| 23 | 253 | 1 | 529 | 1 mod 529 |
| 29 | 406 | -1 | 841 | 840 = -1 mod 841 |
| 31 | 465 | 1 | 961 | 1 mod 961 |
| 37 | 666 | -1 | 1369 | 1368 = -1 mod 1369 |
| 41 | 820 | 1 | 1681 | 1 mod 1681 |
| 43 | 903 | -1 | 1849 | 1848 = -1 mod 1849 |
| 47 | 1081 | 1 | 2209 | 1 mod 2209 |

All verified.

## Generalization

For any prime q (not just 2), the number of q-colorings of pairs is q^{C(p,2)}, and:

q^{C(p,2)} = (q^{(p-1)/2})^p ≡ (q/p) (mod p²)

This gives a family of congruences parametrized by both q and p.

## Related

- THM-305: v_2(T(n)) = (n-1)/2 (the unlabeled consequence)
- THM-307: T(p) mod p decomposition
