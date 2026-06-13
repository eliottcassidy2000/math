# THM-307: Burnside-QR Decomposition of T(p) mod p

**Status:** PROVED
**Found by:** opus-2026-04-05-S25
**Verified in:** `04-computation/qr_tournament_enumeration_s25.py`

## Statement

For every odd prime p:

**T(p) ≡ f_p - (2/p)·w_p - A_p (mod p)**

where:
- f_p = (2^{(p-1)/2} - (2/p))/p is the **Euler quotient of 2**
- w_p = ((p-1)! + 1)/p is the **Wilson quotient**
- A_p = (1/p) sum_{lambda != (p), (1^p)} (p!/z_lambda) 2^{e(lambda)} is the **remaining Burnside sum**
- (2/p) is the Legendre symbol

## Proof

The Burnside formula gives p! * T(p) = sum over all-odd partitions lambda of n:

**Term 1** (identity, lambda = (1^p)):
  Contribution: 2^{C(p,2)}.
  Mod p²: ≡ (2/p) by THM-306.

**Term 2** (p-cycle, lambda = (p)):
  Contribution: (p-1)! * 2^{(p-1)/2}.
  By Wilson: (p-1)! ≡ -1 (mod p), so (p-1)! = -1 + p*w_p.
  By Euler: 2^{(p-1)/2} ≡ (2/p) (mod p), so 2^{(p-1)/2} = (2/p) + p*f_p.
  Product mod p²: (-1 + pw_p)((2/p) + pf_p) = -(2/p) + p((2/p)w_p - f_p) + p²(...)
  So mod p²: ≡ -(2/p) + p((2/p)w_p - f_p).

**Term 3** (all others):
  For lambda != (1^p) and lambda != (p): v_p(z_lambda) = 0 (no part = p, no multiplicity >= p).
  So p | (p!/z_lambda), hence p | each remaining term. Sum = p * A_p.

**Combining:**
  p! * T(p) = (2/p) + [-(2/p) + p((2/p)w_p - f_p)] + p*A_p + O(p²)
            = p((2/p)w_p - f_p + A_p) + O(p²)

Since p! * T(p) must be divisible by p!:
  T(p) = [p((2/p)w_p - f_p + A_p)] / p! = ((2/p)w_p - f_p + A_p) / (p-1)!

By Wilson: (p-1)!^{-1} ≡ -1 (mod p). Therefore:
  T(p) ≡ -((2/p)w_p - f_p + A_p) ≡ f_p - (2/p)w_p - A_p (mod p). QED.

## Special Cases

**Wilson primes** (p = 5, 13, 563): w_p ≡ 0 (mod p), so the middle term vanishes:
  T(p) ≡ f_p - A_p (mod p)

**Wieferich primes** (p = 1093, 3511): q_p(2) ≡ 0 (mod p), hence f_p ≡ 0 (mod p):
  T(p) ≡ -(2/p)w_p - A_p (mod p)

## Verified Cases

| p | (2/p) | f_p%p | w_p%p | A_p%p | formula | T(p)%p | match |
|---|-------|-------|-------|-------|---------|--------|-------|
| 3 | -1 | 1 | 1 | 0 | 2 | 2 | ✓ |
| 5 | -1 | 1 | 0* | 4 | 2 | 2 | ✓ |
| 7 | 1 | 1 | 5 | 2 | 1 | 1 | ✓ |
| 11 | -1 | 3 | 1 | 2 | 2 | 2 | ✓ |
| 13 | -1 | 5 | 0* | 3 | 2 | 2 | ✓ |

*Wilson prime

## Relationship Between Quotients

The Euler quotient and Fermat quotient are related:
  f_p ≡ (2/p) * q_p(2) * 2^{-1} (mod p)

where q_p(2) = (2^{p-1} - 1)/p is the Fermat quotient of 2. Verified for all p <= 23.

## Related

- THM-305: v_2(T(n)) = (n-1)/2
- THM-306: 2^{C(p,2)} ≡ (2/p) mod p²
