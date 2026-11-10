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

*Annotation added after the run (2026-11-10, audit R3; the intervals above are unchanged):* the random rank-8
probability is `0.5590477`, so "random: 0.5591" above should read `0.5590`. Outcome of P1-P2: CONFIRMED
(`out/betarank_classes_10_summary.out`; fractions `0.41890 / 0.55940 / 0.02170`).

---

# Second pre-registration: P3-P6 (mac-mini-2026-11-10-ensemble)

Written and committed before the order-10 run of `biasweights.c` (exact mod-4 bias `b_R` of every vertex cube, with
labelled weights `m!/|Aut(R)|`). In-sample data: orders 2-9 (`out/biasweights_summary_2_9.out`; section 4.1a of
README.md). The program reproduces the audit's labelled rank tables and discrepancies at orders 8 and 9 digit for
digit. Notation: `b_R = E_X (-1)^(q_R(X))`, `rho_m = 2^m E_R[b_R^2]` (by (4.3) the mod-4 correlation of two one-vertex
deletions at order `m+2`, in units of the twin term), `D_n = P(H = 1 mod 4) - P(H = 3 mod 4)` over labelled
tournaments of order `n`.

In-sample: `rho_m = 1.244, 1.156` at `m = 8, 9`; frequency of non-zero bias at corank 2: `1, 0.522, 0.326` at
`m = 4, 6, 8`; `E_R|b_R| = 0.05028, 0.03155` at `m = 8, 9`; `D_7, ..., D_10 = 0.1041, 0.04153, 0.007522, 0.002423`.

**P3 (HYP-9659 (B''), out of sample).** `rho_10` lies in `[1.01, 1.12]`. In particular it is smaller than `rho_9`.

**P4.** Among forms of rank 8 at order 10 (corank 2) the labelled frequency of non-zero bias lies in `[0.255, 0.290]`
(uniformly random quadratic form: 1/4).

**P5.** `D_11` is positive and lies in `[0.0003, 0.0011]`.

**P6.** `E_R|b_R|` at `m = 10` (the bound of Theorem C (4) for order 11) lies in `[0.0220, 0.0238]`.

**Control (must hold by Theorem C).** No odd `mu`, no odd `a_i`; for every class the exact Gauss sum over the `2^10`
out-sets is `0` or `+-2^(10 - rank/2)`, and it is `0` exactly when `q_R` does not vanish on the radical of `beta_R`;
the labelled weights sum to `2^45`.

**Refutation.** `rho_10 >= rho_9 = 1.156` refutes the decay of the excess in (B''); any of P3-P6 outside its interval
is reported as FAILED; a failure of the control would refute Theorem C or the program.

*Outcome of P3-P6 (run after commit 23d6c4033f was pushed; `out/biasweights_10_parts.out`,
`out/biasweights_summary_10.out`).* Control: no failure in 9,733,056 classes; labelled weights sum to `2^45`.
P3: `rho_10 = 1.040628` - CONFIRMED. P4: `0.26407` - CONFIRMED. P5: `D_11 = 260546929/549755813888 = 0.000473932`
- CONFIRMED. P6: `E_R|b_R| = 0.022563` - CONFIRMED.

*Annotation added after the second audit (2026-11-10, `audit2/AUDIT2.md` m8, m11; nothing above is rewritten).* P3
and P6 are not independent of P4: given the labelled rank fractions of order 10, which were known from the first
audit before this pre-registration (`audit/b_rankonly_10_parts.out`, same commit), `rho_10` and `E_R|b_R|` are affine
functions of the non-zero-bias frequencies, dominated by the corank-2 term; the P4 interval maps to
`rho_10 in [1.020, 1.099]` and `E_R|b_R| in [0.02225, 0.02347]`. The independent confirmations are P4 and P5, and
"In-sample data: orders 2-9" above should have added the order-10 rank law. In P6, "the bound of Theorem C (4)"
means the bound inside its proof (`E_R|b_R|`); the theorem as stated has `E_R[2^(-rank/2)]` (`0.05077` at `m = 10`).
The second audit recomputed all four outcomes with independent code (identical).
