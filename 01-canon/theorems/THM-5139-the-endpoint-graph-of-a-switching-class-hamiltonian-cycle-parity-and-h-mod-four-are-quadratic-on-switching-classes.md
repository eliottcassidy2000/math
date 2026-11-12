---
id: THM-5139
title: "The endpoint graph of a switching class. For a tournament T let e_T(v) be the parity of the number of Hamiltonian paths that start or end at v. (G) On a switching class e is affine: e_{T^X} = e_T + Lambda 1_X for every vertex set X, where Lambda is a symmetric F_2-matrix with zero diagonal and zero row sums (an even graph) depending only on the class, with Lambda(u,w) = 1 + <a + b, e_Q> + a^T Lambda_Q b for Q = T - u - w and a, b the out-neighbourhoods of u, w in Q. (C1) hc(T^X) = hc(T) + sum_{a in X} (1 + e_T(a)) + sum_{pairs in X} Lambda(a,b) (mod 2) for the number hc of Hamiltonian cycles. (C2) The parity of the number of cycles of each length, and (H - 1)/2 mod 2, are quadratic functions of the switching set; for any three vertices the eight switchings have sum H = 0 (mod 4). This proves the observation (QF) of HYP-9667. (C3) Recursion for the Hessian M of H mod 4; (C4) at order 6, M(a,b) = 1 + d_a + d_b and every switching class has exactly 20 members with H = 1 and 12 with H = 3 (mod 4); (C5) the classes of (source) + P_19, D_19, P_23, D(P_11) are constant in every odd cycle parity; (C6) the imbalance of a class is 0 or +- 2^(n-1-r/2)"
status: >
  PROVED: (G) by induction on the order from the deletion identity e_T(u) = 1 + sum_{v in N^+(u)} e_{T-u}(v) (mod 2)
  (Redei) - source note, section 2.2; (C1)-(C4), (C6) are corollaries proved there; (C5) is (C2) applied to
  recorded finite data (all single and double switchings of one member of each class,
  05-knowledge/results/character_of_two_records_20261111/out/cycle_parities_pairs_24.out and the order-20 row of
  the table of that note).
  Machine checks: every labelled tournament of orders 3, 4, 5, 6, 7 (2,097,152 at order 7): Lambda symmetric with
  zero diagonal and even rows, equal at T and at every single switching; hc(T) + hc(T^a) = 1 + e_T(a); the second
  difference of hc equals Lambda; every triple sum of hc is even and of H is divisible by 4 - no violation
  (qf_exhaustive.c). Random tournaments of orders 3-10: the affine law, (C1), the recursion, the formula for
  H mod 4 (orders <= 8), (C3) (orders 5-9), (C4) (endpoint_graph.py). Histogram of classes by the number of members
  with H = 1 (mod 4), orders 3-7.
  KNOWN, credited: switching as a Z_2 gauge symmetry of tournaments, its invariants and the class mean of H
  (THM-474, THM-479, THM-1405, THM-1460, THM-1467, THM-1470, THM-4524, THM-4531, Session34, Session272;
  Babai-Cameron 2000); the real Walsh expansion of H on a class, of degree at most 2 floor(N/3), and its Ising
  form for N <= 5 (THM-4524); (QF) as an observation, exhaustive for N <= 6 (HYP-9667, THM-5136, this machine,
  2026-11-11); the quadratic form of a vertex cube (THM-5127 (C)), a different cube; H = 1 + 2 (number of odd
  cycles) (mod 4) (Grinberg-Stanley, arXiv:2307.05569, Theorem 7.1; THM-002; Session29); the character sum of a
  quadratic function over F_2 (Dickson).
  NEW to the repository (not found by the author's survey or by the audit): Theorem (G), the endpoint graph
  and its recursion; (C1); the proof of (QF) and of the quadratic law for every cycle length; (C3)-(C5). In the
  literature nothing was found, with the caveat that Forcade 1973 was not consulted.
  NOT CLAIMED: any value of H produced or excluded; the Chan-Pak conjecture is OPEN. Lambda is not a bijection
  from classes to even graphs (consistent with THM-4531, THM-4722).
  AUDIT: PASS (05-knowledge/results/switching_quadratic_census_gauge_20261111/audit/AUDIT.md; independent, own
  programs): every step of the proof of (G) and of (C1)-(C6) re-derived by hand; no violation on every labelled
  tournament of orders 3-7 (Lemmas 1, 2, the affine law for all X, the recursions, (C1) for all X, the cycle-length
  law, the formula for H mod 4, (C3), (C4), the class histograms, the image sizes 1, 4, 22, 91, 2972); the data
  used by (C5) cover all single and double switchings. Minor: credit Grinberg-Stanley (arXiv:2307.05569, Thm 7.1)
  for H = 1 + 2 c_odd (mod 4); Forcade 1973 (Parity of paths and circuits in tournaments) could not be consulted
  and has to be checked before (C1) is called new in the literature. MISTAKE-766.
session: mac-mini-2026-11-11-gauge
source: 05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md
related:
  - HYP-9669
  - HYP-9667
  - THM-5136
  - THM-5127
  - THM-4524
  - THM-474
  - THM-479
  - THM-1467
  - THM-4531
  - THM-002
---

# THM-5139 - the endpoint graph of a switching class; H mod 4 is quadratic on a class

Source: [the note](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md), section 2.
Independent audit: [AUDIT.md](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/audit/AUDIT.md), PASS
(no major finding; credit and literature caveat added).

## Statement

`T` a tournament on `V`; `T^X` its switching at `X`; `e_T(v)` the parity of the number of Hamiltonian paths starting
or ending at `v`; `hc` the number of Hamiltonian cycles.

- **(G)** `e_{T^X} = e_T + Lambda_T 1_X` for every `X`; `Lambda_T` is symmetric over `F_2`, has zero diagonal and
  zero row sums, and `Lambda_{T^X} = Lambda_T`. For `n >= 4`: `Lambda_T(u,w) = 1 + <a + b, e_Q> + a^T Lambda_Q b`,
  `Q = T - u - w`, `a = N^+(u) cap Q`, `b = N^+(w) cap Q`; `Lambda` is the triangle for `n = 3`.
- **(C1)** `hc(T^X) = hc(T) + sum_{a in X} (1 + e_T(a)) + sum_{{a,b} in X} Lambda_T(a,b) (mod 2)`.
- **(C2)** For every `l`, the parity of the number of `l`-cycles of `T^X` is a quadratic function of `X`, and
  `(H(T^X) - H(T))/2 = sum_{a in X} lambda_a + sum_{{a,b} in X} M(a,b) (mod 2)` with
  `M = sum_{|W| odd} Lambda_{T[W]}`; for any three vertices, `sum_{S in the triple} H(T^S) = 0 (mod 4)`.
- **(C3)** `M_T(a,b) = sum_{v in R : exactly one of a -> v, b -> v} E_R(v) + iota_a^T M_R iota_b`, `R = T - a - b`, `E_R(v)` the parity of the
  number of paths of `R` on an odd number (at least 3) of vertices with `v` as an endpoint (`n >= 4`).
- **(C4)** Order 6: `M(a,b) = 1 + d_a + d_b`; every class has 20 members with `H = 1` and 12 with `H = 3 (mod 4)`.
- **(C5)** Every member of the classes of `(source) + P_23`, `(source) + D(P_11)` has an even number of cycles of
  each odd length and `H = 1 (mod 4)`; every member of the classes of `(source) + P_19`, `(source) + D_19` has
  `c_(2k+1) = C(9,k) (mod 2)` and `H = 3 (mod 4)`.
- **(C6)** On a class, `#{H = 1} - #{H = 3}` (mod 4) is `0` or `+- 2^(n-1-r/2)`, `r = rank M` on `F_2^V/<1>`.

## Proof

Section 2.2 of the note (Lemmas 1, 2, Theorem G, Corollaries 1-6).
