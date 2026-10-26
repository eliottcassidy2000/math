---
id: THM-4806
title: "Named legality languages reconstruct tournaments (no single arc can hide; exact through order 8), Goldbach-free shape lemmas for odd aliquot values, and the row companions of even perfect numbers. (R) For every tournament T on n >= 3 vertices and every pair {u, v}, the arc between u and v lies on a Hamiltonian path of T, or the reversed arc lies on a Hamiltonian path of T^(uv); hence L(T) != L(T^(uv)) for every single flip. L(T) determines the ordered strong components and their languages, so all-order injectivity of T -> L(T) is equivalent to injectivity on strong tournaments; T -> L(T) is injective on all labelled tournaments of order <= 8 (exhaustive). (G) An odd v has a source p^2 q (p, q distinct odd primes) iff (p + 1) | (v - 1) and (v - 1)/(p + 1) - p is an odd prime; a source p^k (k >= 3) iff v - 1 = p + ... + p^(k-1); if v - 1 = 2q' with q' > 3 prime, no such source exists. Census: 4,117,840 Goldbach-dependent odd values below 10^8 (S57's count plus the untouchable 5), decade shares 0.1525, 0.1067, 0.0903, 0.0812 for 10^4..10^7. (P) For N = 2^(p-1)(2^p - 1): (2^a - 1) | (N - 1) iff a | (p - 1); N has a source 2^a q (q an odd prime) iff a | (p - 1) and q_a = (N + 1 - 2^(a+1))/(2^a - 1) is an odd prime (a = p - 1 is N itself; a = 1 gives N - 3, divisible by 5 when p = 3 mod 4); the row companions of the first eight even perfect numbers are 652, 10828, 2^8 x 33685759, 2^6 x 36600682668888127 (p = 5, 7, 17, 31) and none for p = 2, 3, 13, 19. (A) In the two-ray carrier the number of addresses of x >= 2 is v_2(x) + tau(gcd of the exponents of x)."
status: >
  PROVED: (R) (i)-(ii) (Theorem R of the note; Camion's theorem is the only input), (G) (i)-(ii), (P) (i)-(ii), (A).
  FINITE-EXACT: (R) strong form, legality_strong_lemma_20261026.cpp (zero violations among strong tournaments, n <= 7 and X_8
  exhaustive, samples at 9 and 10; the implication strong form => injectivity is PROVED); (R) (iii) legality_language_reconstruction_20261026.cpp (all 2^C(n,2) labelled tournaments for n <= 7
  with the single-pair lemma tested on every dark arc, and all 2^21 members of X_8; zero language collisions; dark-arc
  census); (G) goldbach_dependent_odd_values_20261026.cpp to 10^8 (reproduces S57's 8,250 / 56,266 / 462,731 /
  4,117,839 plus the value 5; five pre-registered predictions filed at 0946a65ace: four confirmed, the 2 x pq pattern
  share 0.1990 missed its range [0.200, 0.240]); (P) perfect_companions_shapes_20261026.cpp (row theorem checked for
  all a <= 2p + 2, p <= 31; shape census with at most four odd prime factors within stated bounds; exact for p <= 13
  through THM-4805); (A) two_ray_addresses_planarity_20261026.py (brute force to 10^4).
  CITED: session70_attachment_guard_alignment_20261026 (A1-A4, the scope statement), strong_atom_continuation_vitali
  (C1-C5; the guard-signature injectivity through n = 7 that (R) (iii) reproduces and extends), HYP-9416 (fresh
  unnamed dynamics, a different operation menu), S57 and aliquot124 (the (S)/(O) split and the distinct-Goldbach
  antecedent), aliquot_guarded_cycles (the row formula), aliquot124 P1-P2 (addresses), THM-4805.
  OPEN: all-order injectivity of T -> L(T) (HYP-9417 (1)); the limit of the Goldbach-dependent share (HYP-9417 (2));
  whether the seventh perfect number has a companion (HYP-9417 (3)).
  AUDITED: pending (note, section 10).
session: opus-2026-10-26-S61
source: 05-knowledge/results/legality_language_goldbach_free_perfect_companions_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md
  - 05-knowledge/results/session70_attachment_guard_alignment_20261026/README.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md
related:
  - 05-knowledge/hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md
  - 05-knowledge/hypotheses/HYP-9417-legality-language-injectivity-goldbach-dependent-decay-and-perfect-hermits.md
  - 05-knowledge/results/aliquot_power_124_20261026/README.md
  - 01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
---

# THM-4806: legality languages, Goldbach-free shapes, row companions

## Theorem R (named legality language)

`L(T) = {p in S_n : p is a directed Hamiltonian path of T}`; an arc is *dark* when it lies on no Hamiltonian path; `T^(uv)` reverses one arc.

(i) For every tournament `T` (`n >= 3`) and every pair `{u, v}`: the arc of `T` between `u` and `v` lies on a Hamiltonian path of `T`, or the reversed arc lies on a Hamiltonian path of `T^(uv)`. (ii) `L(T)` determines the ordered strong components `C_1 -> ... -> C_m` and `L(T) = L(C_1) ... L(C_m)` (concatenation). (iii) `T -> L(T)` is injective on all labelled tournaments of order `<= 8`.

*Proof of (i).* Let `W = T - {u, v}` have strong components `W_1 -> ... -> W_m` and suppose `u -> v` in `T` with both orientations dark. Using Hamiltonian paths of `W` that start at any prescribed vertex of `W_1`, end at any prescribed vertex of `W_m`, or split between `W_1 u ... u W_j` and `W_(j+1) u ... u W_m` at prescribed vertices (Camion), darkness of `u -> v` in `T` gives `W_1 -> v`, `u -> W_m`, and for each `j < m`: no vertex of `W_j` beats `u` or no vertex of `W_(j+1)` is beaten by `v`; darkness of `v -> u` in `T^(uv)` gives the same with `u, v` exchanged. From `W_1 -> u, v` induction yields `W_i -> u, v` for every `i`, so `W_m -> u`, against `u -> W_m`. QED

*Proof of (ii).* The first vertices of the words of `L(T)` form `C_1`; deleting `C_1` from every word gives `L(T - C_1)`. QED

*Proof of (iii).* Exhaustive hashing of the languages (note, section 1); every collision pair at order 8 can be relabelled into `X_8`. QED

The gap to all orders: two tournaments differing on several pairs could both avoid every differing pair although every single flip is visible (the flip path of (i) may use another differing pair). Dark arcs exist in strong tournaments from order 5 (120 of the 544 strong labelled 5-tournaments). **Strong form.** If, for every strong tournament `C` and every dark arc `u -> v` of `C`, `C^(uv)` has a Hamiltonian path through `v -> u` avoiding every other dark arc of `C`, then `T -> L(T)` is injective for all orders (PROVED implication, note section 1); the strong form fails for non-strong tournaments (transitive order 4) and holds for every strong tournament of order `<= 7`, every strong member of `X_8`, and 30,000 + 3,000 random members of `X_9`, `X_10` (FINITE-EXACT / NUMERICAL, legality_strong_lemma_20261026.cpp).

## Theorem G (shape lemmas for odd values)

(i) `s(p^2 q) = v` iff `v - 1 = (p + 1)(p + q)`, i.e. `(p + 1) | (v - 1)` and `q = (v - 1)/(p + 1) - p`; `s(p^k) = v` iff `v - 1 = p + ... + p^(k-1)`. (ii) For `v - 1 = 2q'`, `q' > 3` prime, no source of shape `p^2 q` or `p^k` (`k >= 3`) exists. Census and the reflection-matching reading: note, section 3.

## Theorem P (row companions of even perfect numbers)

`N = 2^(p-1) M`, `M = 2^p - 1` prime. (i) `(2^a - 1) | (N - 1)` iff `a | (p - 1)` (reduce `2^(2p-1) - 2^(p-1) - 1` modulo `2^a - 1`: the residue `2^x - 2^r - 1` with `r = (p - 1) mod a` vanishes only for `r = 0`). (ii) `s(2^a q) = N` with `q` an odd prime iff `a | (p - 1)` and `q_a = (N + 1 - 2^(a+1))/(2^a - 1)` is an odd prime; `a = p - 1` gives `q = M`; `a = 1` gives `q = N - 3`, divisible by 5 when `p = 3 mod 4`. (iii) FINITE-EXACT: row companions `652` (p = 5, a = 2), `10828` (p = 7, a = 2), `2^8 x 33685759` (p = 17, a = 8), `2^6 x 36600682668888127` (p = 31, a = 6); none for `p = 2, 3, 13, 19`; other companions found: `10682 = 2 x 7^2 x 109` (p = 7), `44655764 = 4 x 1669 x 6689` (p = 13); the census is exact for `p <= 13` (THM-4805) and covers at most four odd prime factors otherwise.

## Proposition A (two-ray addresses)

For `x >= 2` the two-ray carrier (rows `2^j m`, columns `m^k`, `k >= 1`) holds `x` at `v_2(x) + tau(g(x))` addresses, `g(x)` the gcd of the exponents of `x`; aliquot124's combined table adds the mixed cells (`36`: four against five).
