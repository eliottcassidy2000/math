---
id: THM-4805
title: "Aliquot sources of 2, 5, 7 and 21, the twin-prime collision on the aliquot table, the powers-and-doublings chain structure, and abc hits on the row of 1. (A) s(n) = sigma(n) - n: s^-1(2) = s^-1(5) = {} (classical), s^-1(7) = {8} (so 7 = 1 + 2 + 4 is forced), s^-1(21) = {18, 51, 91}; every value <= 14 has only sources of the shapes p^2, pq, p^3, and the Mersenne value 2^k - 1 = s(2^k) has the single source 2^k exactly for k = 2, 3. (B) For twin primes p - 2, p: s(p^2) = s(2(p - 2)) = p + 1; the row-type sources 2^a q of p + 1 are exactly the twin source; the collision value is a fixed point only at p = 5 (6 = s(6) = s(25)); 6 is the only even perfect number with a prime-square source; FINITE-EXACT: s^-1(28) = {28}, s^-1(496) = {496, 652}, s^-1(8128) = {8128, 10682, 10828}. (C) With columns P(n) = {n^k} and rows D(n) = {2^j n}: a number is old (on a ray of a smaller number) iff even or a perfect power, the rays of new numbers are disjoint from all earlier rays, n = 2^a r (r odd) has n^2 on the row of r^2, a later number iff a >= 1 and r > 2^a (first instance 6: 36 = 4 x 9), the column of 2^a r is the column of r shifted along the rows, the row of an odd n meets an earlier column iff n is a perfect power, and the row of 1 is the only row that is a column. (D) (1, s(n), n) is a coprime additive triple iff n is almost perfect; (1, 2^k - 1, 2^k) is an abc hit iff 2^k - 1 is not squarefree, iff ord_(p^2)(2) | k for some odd p; the hits with k <= 64 are k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63."
status: >
  PROVED: (A) Lemma 1 (shapes by the number of proper divisors) and case analysis; (B) from the inherited formulas
  s(p^2) = p + 1 and s(2^a q) = (2^a - 1) q + 2^(a+1) - 1 (the guarded prime-extension formula with C = 2^a, CITED);
  (C) unique factorisation; (D) rad(2^k - 1) <= (2^k - 1)/3 when not squarefree. Proofs in the note.
  FINITE-EXACT: aliquot_chains_abc_20261026.py (numpy/sympy, 26 s): the source table to v = 60 (exact, n < v^2), the
  89 untouchables below 1000 (n < 10^6), 35 twin pairs, the perfect numbers to 8128 by the parity lemma, the chain
  bookkeeping to n = 40, the 121 abc hits with c <= 10^4, the Mersenne hits to k = 64 (sympy factorisation), H(P_7) = 189.
  CITED: S57 (s57_untouchables_root27_20261026: census to 10^10, the parity lemma, the Goldbach dependence of odd
  values, the single-source odd values below 10^8, the F_1 reading), aliquot_guarded_cycles_20261014e and the Lehmer
  five note (the affine and dyadic formulas), perfect_core_twin_edges_20261019 (a different twin-prime statement),
  THM-1370 and THM-4790 (7 and 21 never Hamiltonian-path or presentation counts), THM-126 (H(P_7) = 189), the S18 F_21
  proposition, THM-3833 (abc consequences). The untouchability of 2 and 5 and "5 the only odd untouchable under
  distinct Goldbach" are classical (Erdos 1973; S57 section 3).
  ANALOGY (not a theorem): the crevasse dictionary between untouchables and forbidden path counts (note section 5).
  NUMEROLOGY: 7 = 2 + 5, 21 = 4 x 5 + 1, 63 = 7 x 9 = 21 x 3 (note section 5).
  AUDITED: pending (see the note, section 9).
  OPEN: HYP-9415 (twin-value source counts, self-only perfect numbers, squarefree Mersenne numbers as abc non-hits).
session: opus-2026-10-26-S60
source: 05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 05-knowledge/results/aliquot_guarded_cycles_20261014e.md
  - 05-knowledge/results/collatz_aliquot_lehmer_five_20260927.md
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-4790-presentation-count-fibres-7-and-21-are-absent-for-every-tournament-on-at-most-11-vertices.md
  - 01-canon/theorems/THM-126-paley-circulant-uniquely-maximizes-z7.md
  - 01-canon/theorems/THM-3833-abc-conditional-cube-radical-growth-and-hyperbolic-power-finiteness.md
  - 05-knowledge/results/perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md
  - 05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md
  - 05-knowledge/hypotheses/HYP-9415-aliquot-table-predictions-twin-value-sources-self-only-perfect-numbers-mersenne-abc-hits.md
---

# THM-4805: the aliquot table, the chain structure, and abc on the row of 1

Let `s(n) = sigma(n) - n` and `t(n) = tau(n) - 1` (number of proper divisors). A composite `n` has `s(n) >= 1 + sqrt(n)`, so every source of a value `v` is below `v^2` and all source sets below are exact.

**Lemma 1.** `t = 2`: `n = p^2`, `s = p + 1`. `t = 3`: `n = p^3` (`s = 1 + p + p^2`) or `pq` (`s = 1 + p + q`). `t = 4`: `p^4`, `s >= 15`. `t = 5`: `p^5` (`s >= 31`) or `p^2 q` (`s = 1 + p + p^2 + q + pq >= 16`). `t >= 6`: `s >= 22` (six distinct proper divisors sum to at least 21, with equality forcing `60 | n`).

## Theorem A (sources of 2, 5, 7, 21)

(i) `s^-1(2) = s^-1(5) = {}`. (ii) `s^-1(7) = {8}`. (iii) `s^-1(21) = {18, 51, 91}`. (iv) Every value `v <= 14` has only sources of the shapes `p^2`, `pq`, `p^3`; the first value with a four-proper-divisor source is `15 = s(16)`, the first with five is `16 = s(12)`. (v) `2^k - 1 = s(2^k)` has the single source `2^k` exactly for `k = 2, 3`.

*Proof.* (iv) is Lemma 1. (i) Composites have `s >= 3`; `5 = p + 1` needs `p = 4`, `5 = 1 + p + q` needs `2 + 2`, `5 = 1 + p + p^2` has no solution. (ii) `p = 6` fails, `p + q = 6` needs `3 + 3`, `1 + p + p^2 = 7` gives `p = 2`. (iii) `p = 20` fails; `p + q = 20` gives `{3, 17}`, `{7, 13}`; `p^2 + p = 20` fails; `t = 4`: `s(16) = 15`, `s(81) = 40`; `t = 5`: `s(p^2 q) = (1 + p + p^2) + q(1 + p) = 21` only for `p = 3`, `q = 2`; `t >= 6` gives at least 22. (v) The dyadic formula and the exact table (`15 <- {16, 33}`, `31 <- {32, 125, 161, 209, 221}`, `63 <- {64, 177, 817}`); for `k >= 4` a distinct-prime partition of `2^k - 2` gives a semiprime source, verified through `k = 26` by S57's list of single-source odd values below `10^8` (`3, 7, 9, 11, 39`). QED

## Theorem B (the twin collision)

(i) For primes `p - 2, p` (`p >= 5`): `s(p^2) = s(2(p - 2)) = p + 1`. (ii) For a prime `p`, the sources of `p + 1` of the form `2^a q` (`a >= 1`, `q` an odd prime) are exactly `2(p - 2)` when `p - 2` is prime: `(2^a - 1) q = p - 2(2^a - 1)` forces `2^a - 1 | p`. (iii) The value is a source of itself only at `p = 5`: `s(2q) = 2q` iff `q = 3`. (iv) 6 is the only even perfect number with a prime-square source: an even perfect `N > 6` is `1 mod 3`, so `N - 1` is composite. (v) FINITE-EXACT (parity lemma: even sources `m <= 2N - 2`, odd sources are odd squares): `s^-1(6) = {6, 25}`, `s^-1(28) = {28}`, `s^-1(496) = {496, 652}`, `s^-1(8128) = {8128, 10682, 10828}`; the row source `4q` of an even perfect `N` needs `q = (N - 7)/3` prime, true for `N = 28, 496, 8128` and false for the next four; 115 of the 205 twin values `p + 1 <= 10^4` have exactly the two twin sources.

## Theorem C (the chain structure)

Columns `P(n) = {n^k : k >= 1}`, rows `D(n) = {2^j n : j >= 0}`, numbers added in order; a number is old when it lies on a ray of a smaller number. Write `n = 2^a r`, `r` odd.
1. `n` is old iff `n` is even or a perfect power; the rays of a new number meet no ray of a smaller number.
2. `n^2` lies on `D(r^2)`; it already occurs iff `r^2 < n` or `n` is a perfect power, and it is a later doubling (of `r^2 > n`, with `a >= 1`) iff `r > 2^a`. First instance 6; instances to 40: `6, 10, 14, 18, 20, 22, 26, 28, 30, 34, 36, 38`; `X/2 - O(sqrt X)` instances up to `X`.
3. `n^k` lies on `D(r^k)`, a later number iff `r^(k-1) > 2^a`; the column of `2^a r` is the column of `r` shifted along the rows (`(2^a r)^k = 2^(ak) r^k`).
4. The row of an odd `n` meets an earlier column iff `n` is a perfect power (`9 -> 36`, `25 -> 100, 400`, `27 -> 216, 1728, 13824`).
5. `D(1) = P(2)` is the only row that is a column; its segment `1, 2, 4` is the divisor set of 8 and the 1-2-4 decomposition of `7 = s(8)`.

*Proof.* Unique factorisation; the four ray intersections for an odd non-perfect-power `n` and `m < n` each force `m >= n` (note, section 3). QED

## Theorem D (abc on the row of 1)

(i) `(1, s(n), n)` is a pairwise coprime additive triple iff `sigma(n) = 2n - 1`. (ii) `(1, 2^k - 1, 2^k)` has `rad(abc) < c` iff `2^k - 1` is not squarefree, iff `ord_(p^2)(2) | k` for some odd prime `p`; for `k <= 64` exactly `k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63` (multiples of 6, 20, 21). (iii) FINITE-EXACT: 121 abc hits with `c <= 10^4`, 94 containing an odd perfect power, 48 with a power of two as the even term, 30 with `a = 1`; the first hit `1 + 8 = 9` is the row of 1 meeting the column of 3.

*Proof.* (ii) If `2^k - 1` is squarefree, `rad = 2(2^k - 1) >= 2^k`; otherwise `rad(2^k - 1) <= (2^k - 1)/3`. QED

## Relation to {7, 21}

THM-1370 (7 and 21 are never Hamiltonian-path counts) and THM-4790 (never presentation counts for `n <= 11`) are cited, not extended. The proofs on the two sides share a shape (a small value between consecutive terms of every finite family of attained-value sequences, with no rescue by factorisation); that is an ANALOGY. `7 = 2 + 5`, `21 = 4 x 5 + 1`, `63 = 7 x 9 = 21 x 3` (the first non-squarefree Mersenne number, the first abc hit of the row of 1 and the first attained `H` of the form `7 x 3^k`) are NUMEROLOGY. The exact content of "the 1-2-4 decomposition of 7" is Theorem A (ii) with Theorem C (5), together with the cited `{1, 2, 4} = <2> <= (Z/7)^x` as the multiplier group of `F_21 = Aut(P_7)` (THM-126, S18).
