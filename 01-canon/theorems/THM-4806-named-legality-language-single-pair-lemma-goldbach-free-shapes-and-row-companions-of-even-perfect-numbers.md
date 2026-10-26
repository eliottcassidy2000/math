---
id: THM-4806
title: "The named legality language after session70 and THM-4807, Goldbach-free shape lemmas for odd aliquot values, and the row companions of even perfect numbers. (R) For every tournament T on n >= 3 vertices, strong or not, and every pair {u, v}, the arc between u and v lies on a Hamiltonian path of T, or the reversed arc lies on a Hamiltonian path of T^(uv) (single-pair lemma); L(T) determines the ordered strong components and their languages; (R') for a strong tournament and a dark arc u -> v, the arc is the chord of a Hamiltonian path u ... v all of whose arcs are traced, so T^(uv) has a Hamiltonian path through v -> u avoiding every dark arc of T (the strong form; THM-4807 (b)), and the strong form with the reduction gives a second self-contained proof that T -> L(T) is injective at every order (first proofs: session70's guards lane, audited, relative to Busch-Jacobson-Reid and Camion; THM-4807, audited, Redei and Camion); exhaustive labelled census: injective on every labelled tournament of order <= 8. (G) An odd v has a source p^2 q (p, q distinct odd primes) iff (p + 1) | (v - 1) and (v - 1)/(p + 1) - p is an odd prime; a source p^k (k >= 3) iff v - 1 = p + ... + p^(k-1); if v - 1 = 2q' with q' > 3 prime, no such source exists. Census: 4,117,840 Goldbach-dependent odd values below 10^8 (S57's count plus the untouchable 5), decade shares 0.1525, 0.1067, 0.0903, 0.0812 for 10^4..10^7. (P) For N = 2^(p-1)(2^p - 1): (2^a - 1) | (N - 1) iff a | (p - 1); N has a source 2^a q (q an odd prime) iff a | (p - 1) and q_a = (N + 1 - 2^(a+1))/(2^a - 1) is an odd prime (a = p - 1 is N itself; a = 1 gives N - 3, divisible by 5 when p = 3 mod 4); the row companions of the first eight even perfect numbers are 652, 10828, 2^8 x 33685759, 2^6 x 36600682668888127 (p = 5, 7, 17, 31) and none for p = 2, 3, 13, 19. (A) In the two-ray carrier the number of addresses of x >= 2 is v_2(x) + tau(gcd of the exponents of x)."
status: >
  PROVED: (R) (i)-(ii) and (R') (Theorems R and R' of the note; Camion's theorem is the only input; the strong form's
  proof was supplied by the independent audit from the facts in the proof of (i), and is THM-4807 (b) reached
  independently); (G) (i)-(ii); (P) (i)-(ii); (A). All-order injectivity is a second proof: session70's guards lane
  (audited, relative to Busch-Jacobson-Reid and Camion) and THM-4807 (audited, Redei and Camion) proved it first the
  same day; the pre-audit version typed it OPEN (MISTAKE-677).
  FINITE-EXACT: (R) legality_language_reconstruction_20261026.cpp (all 2^C(n,2) labelled tournaments for n <= 7 with
  the single-pair lemma tested on every dark arc, and all 2^21 members of X_8; zero language collisions; dark-arc
  census) and legality_strong_lemma_20261026.cpp (strong form: zero violations among strong tournaments, n <= 7 and
  X_8 exhaustive, 50,280 + 2,195 dark arcs; samples at 9 and 10, seed 0x9E3779B97F4A7C15); (G)
  goldbach_dependent_odd_values_20261026.cpp to 10^8 (reproduces S57's 8,250 / 56,266 / 462,731 / 4,117,839 plus the
  value 5; five pre-registered predictions filed at 0946a65ace: four confirmed, the 2 x pq pattern share 0.1990
  missed its range [0.200, 0.240]); (P) perfect_companions_shapes_20261026.cpp (row theorem checked for all
  a <= 2p + 2, p <= 31, and by the audit for all primes p <= 61; row companions complete; shape census exact for
  p <= 13 through THM-4805, complete for at most three odd prime factors at p = 17, for at most two at p = 19, and
  for the row, q^2 r, q^3 r shapes at p = 31; odd-square companions not enumerated for p >= 17); (A)
  two_ray_addresses_planarity_20261026.py (brute force to 10^4; audit to 2000).
  CITED: session70_attachment_guard_alignment_20261026 (final README: the all-order closure, the alignment-memory
  mechanism, the graph lane), THM-4807 (untraced arcs are cycle-closing chords; the elementary all-order proof; the
  concurrent lanes guard_reconstruction, growing_sieve_path, sieve_path_merge named there), session71's reconciliation
  note, strong_atom_continuation_vitali (C1-C5; guard-signature injectivity through n = 7), HYP-9416 (fresh unnamed
  dynamics, a different operation menu), Busch-Jacobson-Reid 2006 and Busch-Jacobson 2007 (nontraceable arcs =
  dark arcs), S57 and aliquot124 (the (S)/(O) split and the distinct-Goldbach antecedent), aliquot_guarded_cycles
  (the row formula), aliquot124 P1-P2 (addresses), THM-4805, THM-4470 with Lee-Shee.
  OPEN: the limit of the Goldbach-dependent share (HYP-9417 (1)); whether the seventh perfect number has a
  companion (HYP-9417 (2)); the leakage-rank test (HYP-9417 (3)).
  AUDITED: PASS WITH CORRECTIONS (every proof checked and every number reproduced; four majors of scope and
  inheritance and thirteen minors applied; MISTAKE-677; report beside the note).
session: opus-2026-10-26-S61
source: 05-knowledge/results/legality_language_goldbach_free_perfect_companions_20261026.md
depends_on:
  - 01-canon/definitions.md
  - 05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md
  - 05-knowledge/results/session70_attachment_guard_alignment_20261026/README.md
  - 01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md
  - 05-knowledge/results/s57_untouchables_root27_20261026/README.md
  - 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md
related:
  - 05-knowledge/hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md
  - 05-knowledge/hypotheses/HYP-9417-legality-language-injectivity-goldbach-dependent-decay-and-perfect-hermits.md
  - 05-knowledge/results/session71_dissection_norm_alignment_20261026/INCOMING_RECONCILIATION.md
  - 05-knowledge/results/aliquot_power_124_20261026/README.md
  - 01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
---

# THM-4806: legality languages, Goldbach-free shapes, row companions

## Theorem R (named legality language)

`L(T) = {p in S_n : p is a directed Hamiltonian path of T}`; an arc is *dark* (nontraceable) when it lies on no Hamiltonian path; `T^(uv)` reverses one arc.

(i) For every tournament `T` (`n >= 3`) and every pair `{u, v}`: the arc of `T` between `u` and `v` lies on a Hamiltonian path of `T`, or the reversed arc lies on a Hamiltonian path of `T^(uv)`. (ii) `L(T)` determines the ordered strong components `C_1 -> ... -> C_m` and `L(T) = L(C_1) ... L(C_m)` (concatenation). (iii) `T -> L(T)` is injective on all labelled tournaments of order `<= 8` (exhaustive). (iv) `T -> L(T)` is injective at every order (session70; THM-4807; and Theorem R' below with (ii)).

*Proof of (i).* Let `W = T - {u, v}` have strong components `W_1 -> ... -> W_m` and suppose `u -> v` in `T` with both orientations dark. Using Hamiltonian paths of `W` that start at any prescribed vertex of `W_1`, end at any prescribed vertex of `W_m`, or split between `W_1 u ... u W_j` and `W_(j+1) u ... u W_m` at prescribed vertices (Camion), darkness of `u -> v` in `T` gives `W_1 -> v`, `u -> W_m`, and for each `j < m`: no vertex of `W_j` beats `u` or no vertex of `W_(j+1)` is beaten by `v`; darkness of `v -> u` in `T^(uv)` gives the same with `u, v` exchanged. From `W_1 -> u, v` induction yields `W_i -> u, v` for every `i`, so `W_m -> u`, against `u -> W_m`. QED

*Proof of (ii).* The first vertices of the words of `L(T)` form `C_1`; deleting `C_1` from every word gives `L(T - C_1)`. QED

*Proof of (iii).* Exhaustive hashing of the languages (note, section 1); every collision pair at order 8 can be relabelled into `X_8`. QED

## Theorem R' (the strong form, and the second proof of (iv))

(a) Let `C` be strong and `u -> v` dark. With `W = C - {u, v}` and components `W_1 -> ... -> W_m`, darkness gives `W_1 -> v` and `u -> W_m`; `m >= 2` (else `v` is a sink); `W_1` needs an in-neighbour outside `W`, and `v` beats nothing in `W_1`, so `u` beats some `b in W_1`; dually some `a in W_m` beats `v`. Then `P = u, b, ..., a, v` (a Hamiltonian path of `W` from `b` to `a` in the middle) is a Hamiltonian path of `C`, all its arcs are traced, and `b, ..., a, v, u` is a Hamiltonian path of `C^(uv)` through `v -> u` using no dark arc of `C`. (This is THM-4807 (b).) (b) If `L(T) = L(T')` with `T != T'`, by (ii) they have the same strong components; the differing pairs lie inside components and are dark on both sides (an arc of a component is dark in `T` iff dark in the component); for a differing `u -> v` inside `C`, the path of (a) in `C^(uv)` uses `v -> u` and otherwise arcs of `C` outside the differing set, so it is a Hamiltonian path of `C'` through a differing pair, which extends to `T'`: contradiction. QED

The strong form fails for non-strong tournaments (transitive order 4, dark arc `2 -> 0`: both flip paths use a dark arc), which is why (i) alone does not give (iv) and the reduction (ii) is needed. Dark arcs exist in strong tournaments from order 5 (120 of the 544 strong labelled 5-tournaments, one each).

## Theorem G (shape lemmas for odd values)

(i) `s(p^2 q) = v` iff `v - 1 = (p + 1)(p + q)`, i.e. `(p + 1) | (v - 1)` and `q = (v - 1)/(p + 1) - p`; `s(p^k) = v` iff `v - 1 = p + ... + p^(k-1)`. (ii) For `v - 1 = 2q'`, `q' > 3` prime, no source of shape `p^2 q` or `p^k` (`k >= 3`) exists. Census and the reflection-matching reading: note, section 3.

## Theorem P (row companions of even perfect numbers)

`N = 2^(p-1) M`, `M = 2^p - 1` prime. (i) `(2^a - 1) | (N - 1)` iff `a | (p - 1)` (reduce `2^(2p-1) - 2^(p-1) - 1` modulo `2^a - 1`: the residue `2^x - 2^r - 1` with `r = (p - 1) mod a` vanishes only for `r = 0`). (ii) `s(2^a q) = N` with `q` an odd prime iff `a | (p - 1)` and `q_a = (N + 1 - 2^(a+1))/(2^a - 1)` is an odd prime; `a = p - 1` gives `q = M`; `a = 1` gives `q = N - 3`, divisible by 5 when `p = 3 mod 4` and composite for every odd `p <= 31`. (iii) FINITE-EXACT: row companions `652` (p = 5, a = 2), `10828` (p = 7, a = 2), `2^8 x 33685759` (p = 17, a = 8), `2^6 x 36600682668888127` (p = 31, a = 6); none for `p = 2, 3, 13, 19` (complete: each `q_a` is tested); other companions found: `10682 = 2 x 7^2 x 109` (p = 7), `44655764 = 4 x 1669 x 6689` (p = 13); the census is exact for `p <= 13` (THM-4805), and for larger `p` has the scope stated in the status.

## Proposition A (two-ray addresses)

For `x >= 2` the two-ray carrier (rows `2^j m`, columns `m^k`, `k >= 1`) holds `x` at `v_2(x) + tau(g(x))` addresses, `g(x)` the gcd of the exponents of `x`; aliquot124's combined table adds the mixed cells (`36`: four against five).
