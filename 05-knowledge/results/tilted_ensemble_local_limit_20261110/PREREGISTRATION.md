# Pre-registration (mac-mini-2026-11-10-ensemble)

Written and committed before the order-10 run of `betarank.c` (fast mode: rank of `beta_R` only, no out-set loop).
In-sample data: orders 2-9 (section 4.1 of README.md).

**P1 (HYP-9659 (B), out of sample).** Over all 9,733,056 tournament classes of order 10 the rank of
`beta_R = (mu_R / 2 mod 2)` has

- rank 10 in a fraction in `[0.405, 0.430]` (uniformly random alternating 10 x 10 matrix: 0.4197),
- rank 8 in a fraction in `[0.550, 0.575]` (random: 0.5591),
- rank at most 6 in a fraction in `[0.020, 0.032]` (random: 0.0213).

The intervals allow the in-sample drift (top rank slightly under-represented, low ranks over-represented by a factor
of about 1.4-1.7 at orders 8-9).

**P2 (control, must hold by Theorem C (1)).** No pair with `mu_R(i,j)` odd; least `mu_R(i,j)` equal to 2.

**Refutation.** Any of the three fractions of P1 outside its interval refutes HYP-9659 (B) as stated; a failure of P2
would refute Theorem C.
