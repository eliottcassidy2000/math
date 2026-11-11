---
id: THM-5132
title: "The pentagon theorem. For a tournament on n vertices let q_1 be the number of vertices of odd score and q_2 the number of vertices with an odd number of odd-score out-neighbours (s = A1, t = As over F_2). PROVED: (1) c_3 = C(q_1,2) + floor(n/4), c_3 + c_5 = C(q_2,2) + floor(n/4), c_5 = C(q_1,2) + C(q_2,2) (mod 2); (2) on a vertex cube the increment of c_3 + c_5 is the triangle form C(|y|,2) + u.y in the coordinates y = Ax + |X| 1; (3) the cube sum G(R) = sum_X (-1)^((c_3+c_5)(R+v_X)) has a closed form in the kernel of A^T and the solutions of A^T u = 1 over F_2 (weights mod 4), in particular G(R) >= 0 for every tournament R; (4) joint law of (s, t) for a random labelled tournament; (5) closed forms at every order: E(-1)^(c_3+c_5) = 2^(-(n-2)/2) and E(-1)^(c_5) = 2^(-(n-2)) for n = 2 mod 4, and for odd n 2^(-(n-1)/2)(2 + (2/n)2^(-(n-1)/2) - 2^(-(n-2))) and 2^(-(n-1))(3 - (2/n)2^(-(n-3)/2)), with the value of n - 1 at n = 0 mod 4. Hence H mod 4 is decided by q_2 mod 4 for n <= 6. FINITE-EXACT: the race for every set of cycle lengths, n <= 10. The direct continuation to the heptagon fails. Chan-Pak OPEN"
status: >
  PROVED: (1)-(5) (source note, section 6.2: an induction on the order for (1), a change of variables for (2),
  Poisson summation over the image of x -> (Ax, |x|) for (3), a count of tournaments with prescribed out-degree
  parities for (4), and for (5) the sums given by (4) together with a finite-dimensionality argument and an exact
  check for n <= 130). Machine checks: (1) and (4) on every labelled tournament of order <= 6; (1), (3) and the
  auxiliary identities on all 9,932,003 classes of orders 2-10 and on 120,000,000 random labelled bases of orders
  11-40, no failure; the exact labelled sums of the censuses equal the closed forms of (5) for n <= 11.
  KNOWN, credited: tr A^k = k c_k for k = 3, 4, 5 (THM-118); the odd-cycle count is spectral iff n <= 6 (THM-500), which
  is the exhaustive fact d(6) = 2 of THM-1980; the score-parity vector is uniform on its coset for every n
  (Session72 SCORE_PARITY_QUOTIENT.md; inside a switching class for odd n, THM-1470); the classical score formula
  c_3 = C(n,3) - sum C(s_v,2), which gives the first formula of (1) directly; H = 1 + 2 c_odd mod 4 (THM-466,
  THM-478); the vertex-cube quadratic form of the full count and its layers (THM-5127); E_5[{5}] = 7/32 (THM-466).
  RELATED: A^2 1 over the integers is the Kendall-Wei second-order score (THM-896; THM-555: it does not determine
  c_5); A over F_2 and de Caen's rank bound (THM-4857); a partial law for H mod 4 at n = 5, 6 (HYP-6935 [ii]).
  NEW here (the audit searched the repository as well): the weight mod 4 of the second-order score mod 2 and the
  formulas of (1) for c_3 + c_5 and c_5, the second equality of (2), (3), (4), (5), the corollary for n <= 6.
  FINITE-EXACT: the table of E (-1)^(sum_{k in K} c_k) for every set K of odd lengths and n <= 10 (all classes with
  labelled weights); the subset form Y_k = sum_l C(n-l, k-l) c_l = (p_k - C(n,k))/2 mod 2 is PROVED (double counting).
  Pre-registered: R5, R6 (extrapolations from n = 6, 8) FAILED at n = 10; R7, R8 (value at n = 14, 18; equal values at
  n = 11, 12 and 15, 16) CONFIRMED by plain sampling and afterwards PROVED as instances of (5).
  NEGATIVE: the direct continuation, tr(A^3 + A^5 + A^7) = C(|A^3 1|, 2) + floor(n/4), is false (whether another
  vector works is open, HYP-9663 (H) (3)); the cube sum is not the same for every base (0, 8, 16 at order 5).
  NOT CLAIMED: novelty in the literature; any law for cycle lengths >= 7 beyond the table; the Chan-Pak conjecture.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/mod4_race_square_term_20261111/audit/AUDIT.md): every part of
  the theorem re-derived by hand and re-checked with independent programs (all labelled tournaments of order <= 7,
  all classes to order 10, exact census at order 11), no error; minor corrections of credit and wording applied
  (m5-m7, m13, m16, m18). MISTAKE-760 concerns THM-5131.
session: mac-mini-2026-11-11-race
source: 05-knowledge/results/mod4_race_square_term_20261111/README.md
related:
  - THM-5131
  - HYP-9663
  - THM-5127
  - THM-118 (trace-cycle-identity)
  - THM-500 (the-second-spectral-boundary-odd-cycle-count-non-spectral-from-n7)
  - THM-502
  - THM-1980 (the-two-adic-edge-of-the-redei-count)
  - THM-1470 (even-tournaments-are-the-tournament-two-graph-theorem)
  - THM-466 (the-2-adic-digits-of-H-are-odd-cycle-collection-counts)
  - THM-478 (digit-one-degree-law)
verification:
  - 05-knowledge/results/mod4_race_square_term_20261111/k35cube.c
  - 05-knowledge/results/mod4_race_square_term_20261111/k35_tables.py
  - 05-knowledge/results/mod4_race_square_term_20261111/pentagon_closed_forms.py
  - 05-knowledge/results/mod4_race_square_term_20261111/tower_test.py
  - 05-knowledge/results/mod4_race_square_term_20261111/cycrace.c
  - 05-knowledge/results/mod4_race_square_term_20261111/cycle_tables.py
  - 05-knowledge/results/mod4_race_square_term_20261111/c35mc.c
  - 05-knowledge/results/mod4_race_square_term_20261111/c35_tables.py
  - 05-knowledge/results/mod4_race_square_term_20261111/out/
---

# THM-5132 - the pentagon theorem

Source, proofs and tables: [source note](../../05-knowledge/results/mod4_race_square_term_20261111/README.md),
section 6. The Chan-Pak conjecture is **OPEN**.

## Notation

`T` a tournament on `n` vertices; `A` its adjacency matrix over `F_2` (`A_ab = 1` iff `a -> b`); `c_k` the number of
directed `k`-cycles; `s = A1` (indicator of the vertices of odd score), `t = As` (indicator of the vertices with an
odd number of out-neighbours of odd score); `q_1 = |s|`, `q_2 = |t|`. `R + v_X`: a new vertex `v` with out-set `X`
added to `R`. `E_n[K] = E (-1)^(sum_{k in K} c_k)` over uniformly random labelled tournaments on `n` vertices.
`(2/n)` is the Kronecker symbol.

## Statements (PROVED)

**(1) Score formulas.** Modulo 2,

    c_3 = C(q_1, 2) + floor(n/4),      c_3 + c_5 = C(q_2, 2) + floor(n/4),      c_5 = C(q_1, 2) + C(q_2, 2).

**(2) Vertex cube.** For `R` on `m` vertices, `x` the indicator of `X`, `sigma = |X| mod 2`, `y = Ax + sigma 1`,
`u_1 = As`, `u_0 = As + s + 1`:

    (c_3 + c_5)(R + v_X) - (c_3 + c_5)(R) = x^T (A + A^3)(1 + x) = C(|y|, 2) + u_sigma . y      (mod 2).

**(3) Cube sum.** Let `G(R) = sum_X (-1)^((c_3 + c_5)(R + v_X))`, `K_0 = ker A^T`, `K_1 = {u : A^T u = 1}`, and
`n_ij` the number of elements of `K_i` of weight `j (mod 4)` (only `n_00, n_01, n_10, n_13` can be non-zero). Then

| `m mod 4` | `G(R)` |
|---|---|
| 0 | `2^(m/2) |K_0|` |
| 1 | `2^((m-1)/2) (|K_0| + n_13 - n_10)` |
| 2 | `2^(m/2) (n_00 + n_13)` |
| 3 | `2^((m+1)/2) (n_00 + n_10)` |

So `G(R) >= 0` for every tournament; more precisely `G(R) = 0` or `G(R) >= 2^ceil(m/2)`, and `G(R) = 0` only for
`m = 1 (mod 4)`.

**(4) Joint law.** For a uniformly random labelled tournament, `s` is uniform on `{|s| = C(n,2) mod 2}` (KNOWN:
Session72 `SCORE_PARITY_QUOTIENT.md`);
given `s` with support `O`, `1 <= |O| <= n - 1`, `t` is uniform on
`{t : |t on O| = C(|O|, 2), |t off O| = C(n - |O|, 2) (mod 2)}`; `t = 0` if `O` is empty, `t = 1` if `O = V`.

**(5) Exact races.**

| `n` | `E_n[{3}]` | `E_n[{3,5}]` | `E_n[{5}]` |
|---|---|---|---|
| odd | `2^(-(n-1)/2)` | `2^(-(n-1)/2) (2 + (2/n) 2^(-(n-1)/2) - 2^(-(n-2)))` | `2^(-(n-1)) (3 - (2/n) 2^(-(n-3)/2))` |
| `2 (mod 4)` | `2^(-(n-2)/2)` | `2^(-(n-2)/2)` | `2^(-(n-2))` |
| `0 (mod 4)` | `2^(-(n-2)/2)` | the value at `n - 1` | the value at `n - 1` |

**Remark (conditional form).** For `1 <= |O| <= n - 1`,
`E[(-1)^(c_3 + c_5) | s = 1_O] = +-2^(2 - ceil(|O|/2) - ceil((n - |O|)/2))`; for `n = 2 (mod 4)` it equals
`2^(-(n-2)/2)` for every `O`, so `(-1)^(c_3 + c_5)` is uncorrelated with every function of the score parities.

Values: `E_n[{3,5}] = 1/2, 1/2, 13/32, 1/4, 67/256, 67/256, 263/2048, 1/16, 1007/16384` and
`E_n[{5}] = 1, 1, 7/32, 1/16, 11/256, 11/256, 23/2048, 1/256, 49/16384` at `n = 3, ..., 11`.

**Corollary (the race up to order 6).** For `n <= 6`, `H = 1 + 2 (C(q_2, 2) + floor(n/4)) (mod 4)`, and
`D_n = E_n[{3,5}]`: `D_5 = 13/32`, `D_6 = 1/4`.

## The table for all cycle lengths (FINITE-EXACT, `n <= 10`)

| `n` | `{3}` | `{5}` | `{3,5}` | `{7}` | `{3,7}` | `{5,7}` | `{3,5,7}` | all |
|---|---|---|---|---|---|---|---|---|
| 5 | 1/4 | 7/32 | 13/32 | | | | | 13/32 |
| 6 | 1/4 | 1/16 | 1/4 | | | | | 1/4 |
| 7 | 1/8 | 11/256 | 67/256 | 811/8192 | 1553/8192 | 461/8192 | 853/8192 | 853/8192 |
| 8 | 1/8 | 11/256 | 67/256 | 1277/32768 | 1361/32768 | 437/32768 | 1361/32768 | 1361/32768 |
| 9 | 1/16 | 23/2048 | 263/2048 | 16343/2^21 | -1453/2^21 | 9019/2^18 | 35881/2^20 | 63103/2^23 |
| 10 | 1/16 | 1/256 | 1/16 | 35881/2^24 | 9019/2^22 | 9019/2^18 | 35881/2^20 | 650471/2^28 |

At `n = 10` further: `{9}`: `1050641/2^28`, `{3,7,9}`: `533081/2^28`. Not explained: `E_10[{7}] = E_10[{3,5}]
E_10[{3,5,7}]`, `E_10[{3,7}] = E_10[{3,5}] E_10[{5,7}]`, and the equal values of `{5,7}`, `{3,5,7}` at `n = 9, 10`.

## What is not true

- The cube sum of `c_3 + c_5` is not the same for every base (Theorem T's property): `0, 8, 16` at order 5.
- `tr(A^3 + A^5 + A^7) = C(|A^3 1|, 2) + floor(n/4) (mod 2)` fails for 42% of random tournaments of orders 3-14
  (only this direct continuation is refuted).
- "The pentagon is invisible at even order" (true at `n = 6, 8`) is false at `n = 10`.

## Falsifier

A tournament violating (1) or (3), or an order at which a labelled census disagrees with (5).
