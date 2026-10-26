---
id: THM-4805
title: "The aliquot table after session69 and aliquot124. (A') Every value 2 <= v <= 14 of s(n) = sigma(n) - n has only sources of the shapes p^2, pq, p^3 (the t = tau - 1 <= 3 regime; the first four-proper-divisor source is 15 = s(16), the first five-divisor source 16 = s(12)); CITED FINITE-EXACT for 1 <= k <= 26, and CONDITIONAL on distinct-prime Goldbach for every k >= 1: the Mersenne value 2^k - 1 = s(2^k) has the single source 2^k exactly for k = 2, 3 (3 <- {4}, 7 <- {8}; 15 <- {16, 33}, 31 and 63 have five and three sources). (B') For a prime P the sources of P + 1 of the row form 2^a q (q an odd prime) are exactly the twin source 2(P - 2); FINITE-EXACT: s^-1(28) = {28} (28 is its own only source), s^-1(496) = {496, 652}, s^-1(8128) = {8128, 10682, 10828}, s^-1(33550336) = {33550336, 44655764}. (C') In the two-ray carrier (rows 2^j n, columns n^k): old = even or perfect power, the column of 2^a r is the shifted column of r, the row of an odd n meets an earlier column iff n is a perfect power, and for j >= 0 and k >= 1, D(1) = {1} union P(2), while D(2) = P(2) is the only exact row-column equality (the segment 1, 2, 4 of D(1) is the 1-2-4 decomposition of 7 = s(8)). (D') For n >= 2, (1, s(n), n) is a positive coprime additive triple iff n is almost perfect; for k >= 1, (1, 2^k - 1, 2^k) is an abc hit iff 2^k - 1 is not squarefree, iff ord_(p^2)(2) | k for some odd prime p; the hits with k <= 64 are k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63."
status: >
  PROVED: (A') (i), for 2 <= v <= 14, from Lemma 1 (shapes by the number of proper divisors); the unit and prime fibres are separate; (B') (i) from the guarded
  prime-extension formula s(2^a q) = (2^a - 1) q + 2^(a+1) - 1 (CITED); (C') corollaries of session69's earliest-seed
  and row-intersection theorems and aliquot124's address theorems P1-P2 (CITED); (D') rad(2^k - 1) <= (2^k - 1)/3 when
  not squarefree. Proofs in the note.
  CITED FINITE-EXACT: (A') (ii) for 1 <= k <= 26, from S57's odd single-source list below 10^8; not independently recomputed in this bounded repair.
  CONDITIONAL: the all-k assertion in (A') (ii), assuming every even integer at least 8 is a sum of two distinct primes.
  Unconditionally, that global single-source classification is unproved here.
  FINITE-EXACT: aliquot_chains_abc_20261026.py (numpy/sympy, 26 s: the source table to 60, the 89 untouchables below
  1000, the perfect numbers to 8128 by the parity lemma, the Mersenne hits to k = 64, the abc census to 10^4);
  aliquot_sources_preregistered_20261026.cpp (the fifth perfect number, twin-value counts to 10^5, abc census to 10^5,
  even untouchables below 10^5 = 13,862, consistent with A057978).
  CITED (concurrent, same day, audited): session69_prime_ports_power_aliquot_20261026 (complete small fibres with the
  bound n <= (m-1)^2; 2 and 5 untouchable; the conditional route distinct Goldbach => 5 the only odd untouchable =>
  {2, 5} the only prime untouchables; the twin identity s(2p) = p + 3 = s((p+2)^2) with 2p = p + 3 iff p = 3; the
  unique perfect prime-square target (N, p) = (6, 5); the square-later predicate a >= 1, u > 2^a with first instance 6
  and density 1/2; the earliest-seed theorem; the Paley-7 / F_8 bridge; the tournament firewall for {7, 21}; the
  primitive abc carriers), aliquot_power_124_20261026 (address theorems P1-P3, delayed powers, the perfect-number
  diagram, the S57 Goldbach-antecedent correction). Earlier: S57 (census to 10^10, parity lemma, F_1 reading), the
  owner's ray tableau (collatz_mod6_20260922_lead_reflection2 sec. 3), arithmetic_seams_20260921 (forests),
  aliquot_guarded_cycles_20261014e and the Lehmer five note (formulas), THM-1370, THM-4790, THM-126, THM-3833.
  ANALOGY (not a theorem): the crevasse dictionary between untouchables and forbidden path counts.
  NUMEROLOGY: 7 = 2 + 5, 21 = 4 x 5 + 1, 63 = 7 x 9 = 21 x 3.
  AUDITED: full theorem audit remains pending (note, section 7). A bounded domain/quantifier repair is recorded in twin_pythagorean_124_20261026/INCOMING_S60.md; it does not certify unrelated claims.
  OPEN: HYP-9415.
session: opus-2026-10-26-S60
source: 05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 05-knowledge/results/session69_prime_ports_power_aliquot_20261026/README.md
  - 05-knowledge/results/aliquot_power_124_20261026/README.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 05-knowledge/results/aliquot_guarded_cycles_20261014e.md
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-4790-presentation-count-fibres-7-and-21-are-absent-for-every-tournament-on-at-most-11-vertices.md
  - 01-canon/theorems/THM-126-paley-circulant-uniquely-maximizes-z7.md
  - 01-canon/theorems/THM-3833-abc-conditional-cube-radical-growth-and-hyperbolic-power-finiteness.md
  - 05-knowledge/results/collatz_mod6_20260922_lead_reflection2.md
  - 05-knowledge/results/perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md
  - 05-knowledge/hypotheses/HYP-9415-aliquot-table-predictions-twin-value-sources-self-only-perfect-numbers-mersenne-abc-hits.md
---

# THM-4805: the aliquot table after session69 and aliquot124

`s(n) = sigma(n) - n`, `t(n) = tau(n) - 1`. A composite `n` has `s(n) >= 1 + sqrt n`, so every source of a target v > 1 satisfies n < v^2. The fibres at 0 and 1 are handled separately.

**Lemma 1.** `t = 2`: `p^2`, `s = p + 1`. `t = 3`: `p^3` (`1 + p + p^2`) or `pq` (`1 + p + q`). `t = 4`: `p^4`, `s >= 15`. `t = 5`: `p^5` (`s >= 31`) or `p^2 q` (`s >= 16`). `t >= 6`: `s >= 22`. Least sources with `t = 4, 5`: 16 and 12.

## Theorem A' (the t = tau - 1 <= 3 regime; Mersenne single sources)

(i) Every value `2 <= v <= 14` has only sources of the shapes `p^2`, `pq` (distinct primes), `p^3`; the first four- and five-proper-divisor sources are `16 -> 15` and `12 -> 16`. Separately, `s^-1(0)={1}` and `s^-1(1)` is the set of all primes.

(ii) **CITED FINITE-EXACT:** for `1 <= k <= 26`, the Mersenne value `2^k - 1 = s(2^k)` has the single source `2^k` exactly for `k = 2, 3`, by S57's odd single-source list below `10^8`. This cited range was not independently recomputed in the bounded repair. At `k=1` the target1 has infinitely many prime sources. The small fibres include `3 <- {4}`, `7 <- {8}`, `15 <- {16, 33}`, `31 <- {32, 125, 161, 209, 221}`, `63 <- {64, 177, 817}`; 127 has eleven sources. **CONDITIONAL:** if every even integer at least8 is a sum of two distinct primes, the same iff holds for every `k>=1`. Unconditionally, the all-k classification is unproved here.

*Proof.* (i) Lemma 1 applies once targets0 and1 are separated. The only source of0 is1; precisely the primes have aliquot sum1. (ii) The finite statement is cited from S57. Under distinct-prime Goldbach, for each `k>=4` write `2^k-2=p+q` with distinct primes. They are odd, so `pq` differs from `2^k`, and `s(pq)=1+p+q=2^k-1` supplies a second source. The exact fibres at k2,3 and the infinite fibre at k1 finish the conditional statement. QED

## Theorem B' (row-source uniqueness; the perfect numbers)

(i) For a prime `P` and `a >= 1`, `q` an odd prime: `s(2^a q) = P + 1` iff `a = 1` and `q = P - 2`. Proof: `(2^a - 1) q = P - 2(2^a - 1)` forces `2^a - 1 | P`; `2^a - 1 = P` gives `q < 0`. (ii) FINITE-EXACT by the parity lemma: `s^-1(6) = {6, 25}`, `s^-1(28) = {28}`, `s^-1(496) = {496, 652}`, `s^-1(8128) = {8128, 10682, 10828}`, `s^-1(33550336) = {33550336, 44655764}` (`652 = 4 x 163`, `10828 = 4 x 2707`, `10682 = 2 x 7^2 x 109`, `44655764 = 4 x 1669 x 6689`). The row source `4q` needs `q = (N - 7)/3` prime: true for 28 (giving 28 itself), 496, 8128, false for the next four even perfect numbers. (iii) 115 of the 205 twin values `P + 1 <= 10^4` have exactly the two twin sources; in `(10^4, 10^5]` 494 of 1019, mean 2.816, maximum 11.

## Theorem C' (the two-ray carrier)

Rows `D(n) = {2^j n : j>=0}`, columns `P(n) = {n^k : k>=1}`, with positive bases n (so `P(1)={1}`); `n = 2^a r`, `r` odd. (1) `n` is old iff even or a perfect power; the rays of an odd non-perfect-power meet no ray of a smaller number. (2) `n^2` lies on `D(r^2)`, already occurs iff `r^2 < n` or `n` is a perfect power, and is a later doubling iff `a >= 1`, `r > 2^a` (session69's predicate; first instance 6). (3) `(2^a r)^k = 2^(ak) r^k`: the column of an even number is its odd part's column shifted along the rows. (4) The row of an odd `n` meets an earlier column iff `n` is a perfect power. (5) `D(1)={1} union P(2)`, while `D(2)=P(2)` is the unique exact equality between an unrestricted doubling row and a positive-exponent power column. The segment `1,2,4` of `D(1)` is the 1-2-4 decomposition of `7=s(8)`; 2 is the only prime on `D(1)`.

*Proof.* Unique factorisation, as in session69's row-intersection law and aliquot124's P1-P2; (4): `2^j n = m^k` with `n` odd forces `n` to be a `k`-th power, and `n = b^k` gives `(2b)^k`. (5) If `D(n)=P(m)`, then `m>=2`, since `P(1)` is finite. Comparing least elements gives `n=m`, and comparing second elements gives `2n=m^2`, hence `m=n=2`. Conversely `D(2)=P(2)`, and adjoining1 gives `D(1)`. QED

## Theorem D' (abc on the row of 1)

(i) For `n>=2`, `(1,s(n),n)` is a positive pairwise-coprime additive triple iff `sigma(n) = 2n - 1`. At n1, `sigma(1)=1=2*1-1` but the triple `(1,0,1)` is not positive and is excluded. (ii) For `k>=1`, `(1, 2^k - 1, 2^k)` has `rad(abc) < c` iff `2^k - 1` is not squarefree, iff `ord_(p^2)(2) | k` for some odd prime `p`; for `k <= 64` exactly `k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63`. (iii) 121 abc hits with `c <= 10^4` (94 with an odd perfect power, 48 with a power-of-two even term, 30 with `a = 1`), 419 with `c <= 10^5` (323, 119, 63).

*Proof.* (ii) If `2^k - 1` is squarefree, `rad = 2(2^k - 1) >= 2^k`; otherwise `rad(2^k - 1) <= (2^k - 1)/3`. QED

## Relation to {7, 21}

Cited, not extended: THM-1370, THM-4790, and session69's tournament firewall. The shared proof shape (a small value between consecutive terms of every finite family of attained-value sequences, with no rescue by factorisation) is an ANALOGY; `7 = 2 + 5`, `21 = 4 x 5 + 1`, `63 = 7 x 9 = 21 x 3` are NUMEROLOGY.
