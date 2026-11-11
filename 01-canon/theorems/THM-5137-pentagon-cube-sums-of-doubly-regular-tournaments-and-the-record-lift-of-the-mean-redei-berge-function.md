---
id: THM-5137
title: "Two applications. (D) Pentagon sums of doubly regular tournaments, from the kernel formula of THM-5132: every regular tournament has an even number of 5-cycles; a doubly regular tournament R of order n = 7 (mod 8) has rank (n+1)/2 over F_2 with a doubly even kernel, and ALL 2^n one-vertex extensions of R have c_3 + c_5 even; for n = 3 (mod 8) the adjacency matrix is orthogonal over F_2 and the cube sum is 2^((n+1)/2), the least possible. (R) Records (the method of Chen-Mogielnicki-Ono, arXiv:2607.12873, applied to tournaments): the mean Redei-Berge function over all tournaments on n vertices is n! 2^(-n) q_n (Schur's one-row Q-function), the sprout sequence with seed (2+t)/(2-t); its canonical noncommutative lift counts odd-cycle packings by max-ordered type, n! 2^(l(alpha)-n)/pi_u(alpha); labelled tournaments in which the largest vertex reaches every vertex are equinumerous with connected labelled graphs (A001187), and tournaments and graphs are equidistributed by max-ordered type."
status: >
  PROVED: (D) (source note, section 5.1: THM-5132 (1) and (3) with the identities A A^T = (t+1) I + t J, A + A^T =
  J + I over F_2 and Sylvester's rank inequality); (R) (source note, section 6: multiplicativity of the Redei-Berge
  function over ordered sums, Lemmas 3.4 and 4.1 and Proposition 6.1 of the cited paper, the exponential formula).
  Machine checks: (D) for P_q, q = 3, 7, 11, 19, 23, 27, 31, 43, 47, and the Cayley-Dickson sign tournaments CD_k,
  k = 2..7, by the kernel formula, by elimination of the quadratic form and (orders <= 11) by counting 3- and
  5-cycles in all extensions; (R) on every labelled tournament of order <= 5, the identity
  sum_{alpha with odd parts} 2^l(alpha)/pi_u(alpha) = 2 for n <= 22, the exponential formula for n <= 12.
  KNOWN, credited: THM-5132 (this machine, this morning; audited); double regularity of the Paley tournaments and
  of the Cayley-Dickson tower (THM-448, THM-451, THM-5133 (G)); rank (n+1)/2 and the quadratic-residue codes of
  Paley tournaments (THM-481, THM-4857); U_T and the odd-cycle formula (Grinberg-Stanley; THM-002, THM-073, THM-078);
  the scalar forms sum E[H_n] t^n/n! = (2+t)/(2-t) (HYP-1400, HYP-1528) and E[path-cover polynomial] with generating
  function ((2+t)/(2-t))^x (kind-pasteur S21, audit A4); W_X in noncommuting variables (INV-051,
  Mitrovic-Stojadinovic); T = 1/(1 - S) for tournaments and strong tournaments (THM-1862, THM-1936, THM-4817) and
  the relation of T with connected graphs (THM-4746); the question answered by (R) 2 (kind-pasteur S22,
  SOURCES_AND_DOMINATION.md section 7).
  KNOWN in the literature (CITED): (D) 1 is an immediate consequence of the classical identity c_5 + 2 c_4 =
  n(n^2-1)(n^2-9)/160 for regular tournaments (Rowlinson 1986; Savchenko 2016).
  NEW here (the audit searched the repository as well): (D) 1 as a statement in this repository, (D) 2-3; the symmetric-function form of the mean and its
  identification with q_n and the hook Schur functions; the record interpretation of the lift; the reach-set
  reading l_n = A001187(n) and the equidistribution by max-ordered type.
  NOT CLAIMED: novelty in the literature (the mean of U_T is a short computation; connected graphs and tournaments
  sharing an exponential generating function is classical); any consequence for the Chan-Pak conjecture.
  Pre-registered: R4, R5 CONFIRMED (cube sums of P_q and CD_k) and now instances of (D).
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-atoms
source: 05-knowledge/results/character_of_two_records_20261111/README.md
related:
  - THM-5136
  - HYP-9667
  - THM-5132
  - THM-448 (drt-mersenne-doubling-tower)
  - THM-5133
  - THM-481
  - THM-4857
  - THM-073
  - THM-078
  - THM-4817
  - THM-4746
  - THM-589
---

# THM-5137 - pentagon sums of doubly regular tournaments; records

Source: [the note](../../05-knowledge/results/character_of_two_records_20261111/README.md), sections 5.1, 6, 7.
Independent audit: AUDIT-PENDING (pre-audit checkpoint)

## (D) Pentagon sums

Notation of THM-5132: `G_35(R) = sum_X (-1)^((c_3 + c_5)(R + v_X))`, the sum over the `2^m` one-vertex extensions of a
tournament `R` on `m` vertices; `K_0 = ker A^T`, `K_1 = {u : A^T u = 1}` over `F_2`.

1. Every regular tournament has an even number of 5-cycles.
2. Let `R` be doubly regular of order `n = 4t + 3`.
   - `n = 7 (mod 8)`: `rank_2 A = (n+1)/2`; every element of `K_0` has weight `0 (mod 4)`; `K_1 = 1 + K_0`;
     `G_35(R) = 2^n`, i.e. every one-vertex extension of `R` has `c_3 + c_5` even.
   - `n = 3 (mod 8)`: `A A^T = I` over `F_2`; `G_35(R) = 2^((n+1)/2)`, the least value possible at `m = 3 (mod 4)`.
3. `G_35(P_q) = 2^q` if 2 is a square in `GF(q)` and `2^((q+1)/2)` if not; `G_35(CD_k) = 2^(2^k - 1)` for `k >= 3`.

Proof. (1) THM-5132 (1): `c_5 = C(q_1, 2) + C(q_2, 2) (mod 2)`; in a regular tournament `s = A1` is `0` or `1` and
`t = As = s`, so `q_1 = q_2`. (2) `A A^T = (t+1) I + t J`, `A 1 = A^T 1 = (2t+1) 1`. For `t` odd: `A A^T = J`; `u` in
`K_0` gives `|u| 1 = 0` and `C(|u|, 2) = u^T A u = 0`, so `|u| = 0 (mod 4)`; `ker A` and `ker A^T` meet in `0` (their
intersection lies in `ker (J + I) = {0, 1}` and `A 1 = 1`), so `rank A >= (n+1)/2`, and `rank(A A^T) = 1 >= 2 rank A - n`;
THM-5132 (3) for `m = 3 (mod 4)` gives `2^((n+1)/2) (n_00 + n_10) = 2^((n+1)/2) 2^((n-1)/2)`. For `t` even: `A A^T = I`,
`K_0 = 0`, `K_1 = {1}` of weight `3 (mod 4)`, and the formula gives `2^((n+1)/2)`.

## (R) Records

`U_T = sum_S 2^|S| prod_(C in S) p_len(C) p_1^(rest)` (sets `S` of disjoint directed odd cycles; Grinberg-Stanley).
For a composition `alpha` let `s_j` be its partial sums and `pi_u(alpha) = prod_j s_j`; `Psi` denotes the
noncommutative power sums of the first kind.

1. `U` is multiplicative over ordered sums, and over uniformly random labelled tournaments on `n` vertices

       E[U_T] = n! sum_(lambda |- n, odd parts) 2^(l(lambda) - n) p_lambda / z_lambda = n! 2^(-n) q_n
              = (n!/2^(n-1)) sum_(k=0)^(n-1) s_(n-k, 1^k),      sum_n q_n t^n = prod_i (1 + x_i t)/(1 - x_i t).

2. With `U_T^nc = sum_S 2^|S| Psi^(tau(S))`, `tau(S)` the block sizes of `S` (cycles and singletons) listed by
   increasing largest vertex: `E[U_T^nc] = n! sum_(alpha |= n, odd parts) (2^(l(alpha) - n)/pi_u(alpha)) Psi^alpha`, the
   canonical lift of `E[U_T]` in the sense of Chen-Mogielnicki-Ono. In particular
   `sum_(alpha |= n, odd parts) 2^l(alpha)/pi_u(alpha) = 2`.
3. Let `l_n` be the number of labelled tournaments on `[n]` in which the vertex `n` reaches every vertex. Then
   `sum 2^C(n,2) x^n/n! = exp(sum l_n x^n/n!)`, so `l_n = A001187(n) = 1, 1, 4, 38, 728, 26704, 1866256, ...` (connected
   labelled graphs); `l_n = sum_j C(n-1, j-1) S_j 2^C(n-j, 2)` with `S_j` the number of strong labelled tournaments; and
   the number of labelled tournaments whose pieces (cut the sequence of strong components before every component
   whose largest vertex exceeds all earlier ones), listed by increasing largest vertex, have sizes `alpha` is
   `prod_j C(s_j - 1, alpha_j - 1) l_(alpha_j)`, the number of labelled graphs whose components have that type.

Proofs: source note, section 6.

## Remark (elementary): the monoid of path counts has elasticity 5/2

With the atoms of `M = {odd} \ {7, 21}` (KNOWN, THM-4094: the primes other than 7, the numbers `7p` with `p` an odd
prime other than 3, 63, 343): the largest ratio of two factorisation lengths of one element is `5/2`, attained
exactly at the even powers of 63 (`63^2 = 3^4 * 49`); the smallest element with factorisations of different lengths
is `315`. For the actual set of path counts the lower bound is unconditional and the upper bound holds if the
Chan-Pak conjecture does (source note, section 7).

## Falsifier

A doubly regular tournament violating (D) 2; a labelled tournament count disagreeing with (R) 2 or 3.
