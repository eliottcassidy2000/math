---
id: HYP-9357
title: "Paley extremality stops where Kantor's 2-transitive designs stop: for prime q = 3 mod 4 the Paley tournament P_q maximises the number of Hamiltonian paths among circulants, and conjecturally among all tournaments of order q, iff q is 3, 7 or 11. These are exactly the q = 3 mod 4 for which the design of translates of the quadratic residues is 2-transitive (trivial, the Fano plane under PSL(2,7), the (11,5,2) biplane under PSL(2,11)). For q >= 19 the rotational circulant Z_q{1..(q-1)/2} has more, as in Alon's and Wormald's conjecture that the rotational tournament is the maximiser for large odd n. At n = 13 a rigid regular tournament beats every circulant."
status: >
  OPEN: the iff, and the global form for q >= 19.
  FINITE-EXACT:
  - the circulant maxima for every odd order 3 <= n <= 23 (circ_counts.c, all multiplier classes; reproduced by the
    independent audit); P_q is the circulant maximum for prime q = 3, 7, 11 and is beaten by the rotational circulant
    at q = 19 and 23;
  - at n = 13 a regular tournament with |Aut| = 1 and H = 3719893 (local search, 160 restarts, found 51 times with its
    converse; H, C, strongness and |Aut| re-verified by the census DP) beats every circulant (best 3711175, the
    rotational one), and its 168361 Hamiltonian cycles beat every circulant too (best 166848).
  CITED:
  - global maxima for n <= 11 (OEIS A038375: 3, 15, 189, 3357, 95095 at n = 3, 5, 7, 9, 11, all attained by circulants);
  - W. M. Kantor, Classification of 2-transitive symmetric designs, Graphs Combin. 1 (1985) 165-166, DOI 10.1007/BF02582940;
  - Galois' exceptional actions of PSL(2,p) on p points (p <= 11);
  - N. Alon, The maximum number of Hamiltonian paths in tournaments, Combinatorica 10 (1990) 319-324;
  - N. C. Wormald, Tournaments with many Hamilton cycles (preprint), who conjectures that the rotational tournament
    maximises for large odd n and whose blow-ups give the constant about 2.8559;
  - Adler-Alon-Ross, Random Structures Algorithms 18 (2001), for the lower bound (e - o(1)) n!/2^(n-1).
  PROVED (elementary, here, given Kantor): the Paley design is 2-transitive iff q = 3, 7, 11 (argument below,
  repaired after the audit, MISTAKE-651; nauty confirms that it is not 2-transitive for q = 19, 23, 27, 31, 43, 47).
  AUDITED 2026-10-21: citation corrected; the design-isomorphism step repaired; prior conjecture credited;
  "circulant" restricted to prime q.
  CONTINUATION (2026-10-21, mac-mini-2026-10-21-fanopsp with opus-2026-10-21-S51; two independent methods): the
  decisive tests at q = 27 (prime power) and q = 31 PASS prediction (3): H(R_27) = 426990256001863046739 >
  H(P_27) = 418425945688463432031 (ratio 1.0205) and H(R_31) = 20351875726012821643802839 > H(P_31) =
  19876729650849384218015865 (ratio 1.0239). Method A (mac-mini, paley_extremal_20261021/ie_ham.c): inclusion-exclusion
  over vertex subsets in Z/2^128 with necklace (circulants) or F_27-translation-orbit reduction, validated on P_7, R_7,
  P_11, R_11, P_19, P_23. Method B (S51, repunits_hyp9357_20261021/hyp9357_rotational_vs_paley.py): Karp inclusion-exclusion
  modulo three 58-bit primes. H(P_27) and H(P_31) also equal |Aut| times S50's orbit counts (1053 * 397365570454381227,
  465 * 42745655163116955307561).
  MECHANISM (2026-10-22, THM-4751 (2), (6), (7); HYP-9362):
  - P_q is B_q = ((q-1)/2 - (2-(2/q))h(-q))/2 flips from R_q, and B_q <= 1 iff q in {3, 7, 11}.
  - The circulant H-maximum is the additive-energy minimum for n <= 11 and R_n for 13 <= n <= 25. The complete n = 25
    table is new, and R_n is a strict local maximum to n = 31.
  - In H = sum 2^k alpha_k the flat class wins alpha_1, alpha_2 and R_n wins every alpha_k with k >= 3. At n = 13 R_n
    overtakes only through alpha_4, i.e. once four disjoint triangles fit. This is a tournament reason for the
    cut-off at 11, with the small-numbers caveat.
session: mac-mini-2026-10-21-fanopsp
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
related:
  - 01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md
  - 01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md (P_q converse-Gluck exceptional iff q = 3, 7)
  - 01-canon/theorems/THM-211-paley-design.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (top structure of the spectrum)
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/circ_counts.c (+ circ_*.out, paley19.out, paley23.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/hmax_search.c (+ hmax_search_13.out, hmax_search_13b.out, hmax_search_15.out, hmax13_best.txt)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/audit/a18_thm4747_checks.py, a20_designs.py (independent audit)
---

# HYP-9357 — Paley extremality and the exceptional primes 3, 7, 11

## Data (FINITE-EXACT)

The maximum of `H` over circulant tournaments `Cay(Z_n, S)` (all connection sets) is:

| `n` | circulant maximiser | `H` | `H·2^(n−1)/n!` | `H(P_n)` if Paley exists |
|---|---|---|---|---|
| 3 | C_3 | 3 | 2.00 | 3 |
| 5 | Z_5{1,2} | 15 | 2.00 | — |
| 7 | **P_7** | 189 | 2.40 | 189 |
| 9 | Z_9{2,3,4,8} | 3357 | 2.37 | — |
| 11 | **P_11** | 95095 | 2.44 | 95095 |
| 13 | rotational | 3711175 | 2.44 | — |
| 15 | rotational | 198464295 | 2.49 | — |
| 17 | rotational | 13689269499 | 2.52 | — |
| 19 | rotational | 1184212824763 | 2.55 | 1172695746915 |
| 21 | rotational | 125547534942879 | 2.58 | — |
| 23 | rotational | 16011537490557279 | 2.60 | 15760206976379349 |
| 27 | rotational only (full circulant maximum not computed) | 426990256001863046739 | 2.63 | 418425945688463432031 (`P_27` on `F_27`) |
| 31 | rotational only (full circulant maximum not computed) | 20351875726012821643802839 | 2.66 | 19876729650849384218015865 |

* *Rotational* means `Z_n{1, …, (n−1)/2}`.
* For orientation, Wormald's blow-ups of rotational tournaments give the best known asymptotic constant for the maximum, about 2.8559. The Adler–Alon–Ross bound is `e`.
* At `n = 3, 5, 7, 9, 11` the circulant maximum equals the global maximum (A038375).
  * Their automorphism orders are 3, (5 and 3), 21, 9, 55.
  * The class is unique at `n = 7, 9`; at `n = 5` two classes reach 15.
* **At `n = 13` the circulant maximum is not global.**
  * A regular tournament with **trivial automorphism group** has `H = 3719893`, against `3711175` for the best circulant. Its converse has the same `H`; both were found 51 times in 160 restarts of a steepest ascent over cyclic-triangle reversals.
  * It also has more Hamiltonian cycles than any circulant: 168361 against 166848.
  * This gives the lower bound `max H(13) ≥ 3719893`. OEIS A038375 stops at `n = 11`.
* At `n = 15`, 16 restarts found nothing above the rotational circulant. That search is weak; this is not evidence either way.
* The Paley tournament wins among circulants at 7 and 11 and loses at 19, 23 and 31. At the prime power 27 it also loses to the rotational circulant.
  * Rotational/Paley ratios: `1.0098, 1.0159, 1.0205, 1.0239` at `q = 19, 23, 27, 31`; the gap grows.
  * At `n = 27, 31` only the rotational circulant was computed, not the full circulant maximum.

## Conjecture

1. For prime `q ≡ 3 (mod 4)`, `P_q` is `H`-maximal among all tournaments of order `q` iff `q ∈ {3, 7, 11}`.
2. Equivalently, by Kantor's classification: iff the Paley design of translates of `QR_q` is 2-transitive.
3. The rotational circulant beats `P_q` for every prime `q ≥ 19` with `q ≡ 3 (mod 4)`.

**Credit.** (3) and the "only if" half of (1) for large `q` are special cases of the conjecture of Alon (1990) and Wormald: the rotational tournament maximises Hamilton cycles, and analogously paths, for large odd `n`.

* The naive all-`n` form fails at small `n`: `C(P_7) = 24 > C(R_7) = 17`, and the `n = 9` maximum is `Z_9{2,3,4,8}`, not `R_9`.
* The `n = 13` record shows that "large" must exceed 13, for paths and for cycles.

## The design statement (PROVED given Kantor; repaired after the audit)

The Paley design has parameters `(q, (q−1)/2, (q−3)/4)`. Kantor's 2-transitive symmetric designs, up to complementation, are:

* the point-hyperplane designs of `PG(m−1, r)`;
* the unique `(11,5,2)` biplane;
* the `(176,50,14)` Higman design;
* the symplectic designs on `2^(2m)` points.

The Higman and symplectic designs have even order. A point-hyperplane design (or its complement) with `k = (v−1)/2` forces `r = 2`, so `q = 2^m − 1`. This is prime, because `2^m − 1` is never a proper prime power.

**Matching against `PG(m−1, 2)` for `q = 2^m − 1 ≥ 31`.**

* `q ∥ |PGL(m,2)|`, so every regular cyclic subgroup is a Sylow `q`-subgroup, and they are all conjugate.
* A design isomorphism can therefore be chosen to intertwine the translation groups. It is then affine, so it is an equivalence of difference sets.
* Equivalent difference sets have equal multiplier groups.
  * For the Singer set, the multiplier group is `⟨2⟩`, of order `m`. The automorphism group of `PG(m−1,2)` is `PGL(m,2)` for `m ≥ 3`, and the normaliser of a Singer cycle there is the cycle extended by the Frobenius, of order `qm`.
  * For `QR_q`, the multiplier group is `QR_q`, of order `(q−1)/2`.
* `(q−1)/2 = m` forces `q = 7`.

**Conclusion.** The Paley design is 2-transitive exactly for `q = 3` (trivial), `7` (Fano) and `11` (biplane). The audit confirmed with nauty that the designs for `q = 31` and `127` are not isomorphic to `PG(4,2)` and `PG(6,2)`.

## Why this is a live alignment

**Three Paley lists nest.**

* `{3, 7}`: converse-Gluck exceptional (THM-4743).
* `{3, 7, 11}`: `H`-maximal and 2-transitive design.
* Galois' exceptional primes are those `p` for which `PSL(2,p)` has a subgroup of index `p`: `p ∈ {2, 3, 5, 7, 11}` (via `A_4 < A_5`, `S_4 < PSL(2,7)`, `A_5 < PSL(2,11)`). Restricted to `q ≡ 3 (mod 4)`, the list is exactly `{3, 7, 11}`.

**The index-`q` actions are design actions.**

* For `q = 7` and `q = 11` the index-`q` action is the action on the out-neighbourhood design of `P_q`: the Fano plane, resp. the biplane.
* The octahedral stabiliser `S_4` of a Fano point is the full symmetry group of the tetrahedron, whose six edges plus one formal point give the Fano plane (session23's tetrahedral bridge).

**Symmetric maximisers stop at 11.** The rigid `n = 13` record fits a picture in which design symmetry (`PSL(2,7)`, `PSL(2,11)`) balances the 5- and 6-vertex statistics and keeps `P_7`, `P_11` on top. Beyond 11 no such symmetry is available and the maximisers lose all symmetry. This explanation is a HEURISTIC; no inequality is proved.

## Cheapest decisive tests

* DONE: `H(P_27)` (Paley over `F_27`; prime power) and `H(P_31)` against the rotational circulant. Both lose (table above), by two independent exact methods (`paley_extremal_20261021/ie_ham.c`, Z/2^128 inclusion-exclusion with necklace or `F_27`-translation reduction; S51's `repunits_hyp9357_20261021/hyp9357_rotational_vs_paley.py`, Karp inclusion-exclusion modulo three 58-bit primes).
* Still open: `H(P_43)`. The same inclusion-exclusion needs about `2^43/43 · 43^3 ≈ 1.6·10^16` operations, out of reach here; a meet-in-the-middle or transfer-matrix count would be needed.
* The full circulant maximum at `n = 25, 27, 29, 31` (all multiplier classes), to see whether the rotational circulant is still the circulant maximum there.
* A stronger search at `n = 15, 17` (simulated annealing, larger move sets) for non-circulant tournaments above the rotational one.
* NEGATIVE CONTROL (2026-10-21, spark from the wide-net survey). Does the global maximum on `q + 1` vertices lie in the switching class of "`P_q` + sink"? That class is the Legendre chirotope class of `P¹(F_q)` (THM-4602), on which `PSL(2,q)` acts 2-transitively.
  * REFUTED: the class maxima are 3, 333 and 256023 at `q = 3, 7, 11`, against the global maxima `H_max(4) = 5` and `H_max(8) = 661`.
  * The class spectra are `{3}`, `{189, 333}` and `{95095, 196815, 236447, 256023}`. The minimum is `H(P_q)`, attained at "`P_q` + sink" itself.
  * Script: `paley_extremal_20261021/paley_sink_switching.c` (+ `.out`). So 2-transitivity on `q + 1` points does not produce `H`-maximisers.
