---
id: HYP-9275
title: "(A) Paley signed self-avoiding run sums through 12 steps lie in Z[q] + Z[q] a_q, a_q the Beukers-supercongruence form eta(2z)^4 eta(4z)^4: the only arithmetic cores are the wheel and the octahedron, U(K_{2,2,2}) = q(q-1)((5-q)a_q + 3q); A_12 = P_12 + 28q(q-1)(17q-251)a_q exactly. (B) The Paley-graph side (q = 1 mod 4) also sees the CM form squared: U(K_{2,4} + 2K_2) = q(q-1) b_q^2. (C) The K6-minus-P6 terminal: Y_q = |b(q)|^2 - q^2 + 2chi_{-4}(q)q with b the weight-3 level-32 newform (mac-mini's HYP-9276; independently confirmed here), and the 2^7 signature traced to the congruence |b(q)| = +-(q+1) mod 2^m."
status: >
  OPEN (the octahedron identity and the general statement of (A); the congruence in (C)).
  PROVED: the partition-lattice identity A_L = sum_pi mu(pi) U(path/pi); THM-4652's reduction calculus; vanishing of
  K_4, K_5 and three six-vertex cores by odd-reversal automorphisms (THM-4652 part 4); the wheel identity (Ahlgren-Ono
  corollary; S34).
  FINITE-EXACT:
  A_2..A_8 reproduced; A_10 and A_12 at 14 primes, with the forced zeros A_10(7) = 0, A_12(7) = A_12(11) = 0 and the
  out-of-sample check A_12(3) = 0;
  the octahedron identity at all 26 odd primes 5..107 (both classes mod 4);
  U(K_{2,4}+2K_2) = q(q-1) b_q^2 at 26 primes (b_q = 4a^2 - 2q for q = a^2 + b^2, a odd; 0 for q = 3 mod 4);
  (C): the square law at all 124 primes 3..691 (exact O(q^3) table) and at all 7,373 primes q = 3 (mod 4) in
  [30,000, 200,000] (FFT evaluator of THM-4653); the depth table and the Mersenne values.
  REFUTED in-session: "Y_q = 1 at Mersenne primes" (Y_127 = 1 - 128^2), and the 2-adic bound
  v_2(Y_q - 1) >= 2 v_2(q+1) + 3 (fails at 127).
session: opus-2026-10-14-S35
source: 05-knowledge/results/paley_least_resistance_runs_one_modular_form_20261014.md
scripts:
  - 04-computation/experiments/paley_run_partition_compiler_20261014.py
  - 04-computation/experiments/paley_A12_assembly_20261014.py
  - 04-computation/experiments/paley_six_vertex_cores_20261014.py
  - 04-computation/experiments/paley_core_symmetry_20261014.py
  - 04-computation/experiments/legendre_twisted_moment_apery_20261014.py
  - 04-computation/experiments/paley_k6p6_fast_table_20261014.py
  - 04-computation/experiments/paley_k6p6_level32_form_20261014.py
  - 04-computation/experiments/paley_k6p6_two_adic_mersenne_20261014.py
  - 04-computation/experiments/paley_k6p6_two_adic_valuation_law_20261014.py
related: [THM-4652, THM-4653, HYP-9274, HYP-9276, HYP-9278, THM-438, HYP-9273]
---

# HYP-9275

## (A) One modular form through twelve steps

**Statement.** For prime powers `q ≡ 3 (mod 4)` and `L <= 12`, the signed self-avoiding run sum `A_L(q)` of the Paley tournament lies in `Z[q] + Z[q]·a_q`, with `a_q` the Fourier coefficient of `η(2z)⁴η(4z)⁴`. For prime powers, read `a_q` as the Frobenius trace. Explicitly:
- `A_10 = 2q(q−1)(21q⁴ − 639q³ + 7028q² − 32098q + 48662) + 34q(q−1)·a_q`;
- `A_12 = 4q(q−1)(33q⁵ − 1540q⁴ + 28243q³ − 249550q² + 1027794q − 1488872) + 28q(q−1)(17q − 251)·a_q`.

The leading coefficients of `A_2, …, A_12` are the Catalan numbers `1, 2, 5, 14, 42, 132` (FINITE-EXACT; the mechanism is heuristic).

**Mechanism.** The partition-lattice compiler meets the following cores:
- `K_4`, `K_5` and three six-vertex cores, all killed by THM-4652;
- the wheel `W_4`: `U = −q(q−1)(q + a_q)`, PROVED;
- the octahedron `K_(2,2,2)`: `U = q(q−1)((5−q)a_q + 3q)`, FINITE-EXACT at 26 primes.

In `A_12` the wheel coefficient is `6688 − 408q` and the octahedron coefficient is `−68`.

**Open.**
- A proof of the octahedron identity. Planar duality (HYP-9278) gives the equivalent `U(Q_3) = q·U(K_(2,2,2))`.
- Whether `A_14` and beyond stay in `Z[q] + Z[q]a_q`.

## (B) The Paley-graph side (`q ≡ 1 (mod 4)`)

Here the odd-reversal symmetry is absent (`χ(−1) = +1`). The same pattern sums then carry the weight-3 CM form `η(4z)⁶`, whose coefficient is `b_q = 4a² − 2q` (`q = a² + b²`, `a` odd), and its square:
- `U(K_4) ∝ b_q` (classical; Evans–Pulham–Sheehan);
- `U(K_(2,4) + 2K_2) = q(q−1)·b_q²` (FINITE-EXACT, 26 primes).

The other two six-vertex cores, `K_6 − (P_3 + K_2)` and `K_(3,3) + e`, are nonzero for `q ≡ 1 (mod 4)` and are not identified. The octahedron identity holds in both classes.

## (C) The `K_6 − P_6` terminal

The identification is mac-mini-2026-10-14-wilson's:
- HYP-9276, at all 3,243 primes to 30,000 and 12 prime powers, in the injective normalisation;
- THM-4653 for the PROVED reduction.

This session reached the same law independently, from its exact table, in the unrestricted normalisation `Y_q = U/(q(q−1))`:

`Y_q = |b(q)|² − q² + 2χ₋₄(q)q`, with `b = q + 4iq³ + 2q⁵ − 8iq⁷ − 7q⁹ − …` the non-CM newform of `S_3(Γ_0(32), χ₋₄)`.

- **Construction.** `b` was rebuilt from the 4-dimensional span of level-32 weight-3 eta quotients with character `χ₋₄`, where `T_3` has characteristic polynomial `x²(x² + 16)`.
- **Normalisations.** `W_q(HYP-9276) = Y_q + 7 + 14χ₋₄(q) − 2χ₋₄(q)q`, checked by brute force at 7 primes.

**2-adic structure (FINITE-EXACT to 200,000).**
- `Y_q − 1 = (|b(q)| − q − 1)(|b(q)| + q + 1)`.
- `|b(q)| ≡ ±(q+1) (mod 2^m)`, with `m = 4, 5, 6, 6` for `q ≡ 11, 3, 7, 15 (mod 16)`.
- Hence `v_2(Y_q − 1) >= 7, 8, 10, 11`, with every minimum attained.
- At Mersenne primes `2^p − 1`, `|b(q)| = 8, 32, 0, 0, 2^14, 63·2^13` for `p = 3, 5, 7, 13, 17, 19`.
- The valuation law `v_2(b(q)) >= min(v_2(q+1), K)` is reported in the companion atlas note.
- `Y_q = 1` iff `|b(q)| = q + 1` (`q = 3, 7, 31, 199, 223, …`).

**Open.** Explain the congruence. The suggested route is McCarthy–Papanikolas' trace formula, which runs over Legendre curves with `Z/2 × Z/8` torsion.
