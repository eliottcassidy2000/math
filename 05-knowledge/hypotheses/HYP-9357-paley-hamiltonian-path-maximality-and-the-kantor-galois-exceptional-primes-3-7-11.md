---
id: HYP-9357
title: "Paley extremality stops where Kantor's 2-transitive designs stop: the Paley tournament P_q maximises the number of Hamiltonian paths (among circulants, and conjecturally among all tournaments of order q) iff q is 3, 7 or 11, which are exactly the q = 3 mod 4 for which the design of translates of the quadratic residues is 2-transitive (trivial, the Fano plane under PSL(2,7), the (11,5,2) biplane under PSL(2,11)); for q >= 19 the rotational circulant Z_q{1..(q-1)/2} has more"
status: >
  OPEN (the iff, and the global form for q >= 19). FINITE-EXACT: the circulant maxima for every odd order 3 <= n <= 23
  (circ_counts.c, all multiplier classes); P_q is the circulant maximum for q = 3, 7, 11 and is beaten by the
  rotational circulant at q = 19 and 23. CITED: the global maxima for n <= 11 (OEIS A038375: 3, 15, 189, 3357, 95095
  at n = 3, 5, 7, 9, 11, all attained by circulants found here); Kantor's classification of 2-transitive symmetric
  designs (W. M. Kantor, J. Combin. Theory A 38 (1985) 66-74); Galois' exceptional actions of PSL(2,p) on p points
  (p <= 11). PROVED (elementary, here, given Kantor): the Paley design (q, (q-1)/2, (q-3)/4) can match only the
  point-hyperplane design of PG(m-1, 2) (q = 2^m - 1) or the (11,5,2) biplane, the Higman and symplectic designs
  having even order; for prime q = 2^m - 1 >= 31 the quadratic-residue set is not equivalent to a Singer set,
  because equivalent difference sets have equal multiplier groups and these have orders (q-1)/2 and m. So the
  Paley design is 2-transitive iff q = 3, 7, 11.
session: mac-mini-2026-10-21-fanopsp
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
related:
  - 01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md
  - 01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md (P_q converse-Gluck exceptional iff q = 3, 7)
  - 01-canon/theorems/THM-211-paley-design.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (top structure of the spectrum)
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/circ_counts.c (+ circ_*.out, paley19.out, paley23.out)
---

# HYP-9357 — Paley extremality and the exceptional primes 3, 7, 11

## Data (FINITE-EXACT)

The maximum of `H` over circulant tournaments `Cay(Z_n, S)` (all connection sets) is:

| `n` | maximiser | `H` | `H·2^(n−1)/n!` | `H(P_n)` if Paley exists |
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

* *Rotational* means `Z_n{1, …, (n−1)/2}`.
* At `n = 3, 5, 7, 9, 11` the circulant maximum equals the global maximum (A038375). Those are the only orders where the global value is known.
* The Paley tournament wins at 7 and 11 and loses at 19 and 23.

## Conjecture

1. `P_q` is `H`-maximal among all tournaments of order `q` iff `q ∈ {3, 7, 11}`.
2. Equivalently, by Kantor's classification: iff the design of translates of `QR_q` is 2-transitive.
3. The rotational circulant beats `P_q` for every prime `q ≥ 19` with `q ≡ 3 (mod 4)`.

## Why this is a live alignment and not only a coincidence

**Three Paley lists nest.**

* `{3, 7}`: converse-Gluck exceptional (THM-4743).
* `{3, 7, 11}`: `H`-maximal and 2-transitive design.
* Galois' exceptional primes are those `p` for which `PSL(2,p)` has a subgroup of index `p`: `p ∈ {2, 3, 5, 7, 11}` (via `A_4 < A_5`, `S_4 < PSL(2,7)`, `A_5 < PSL(2,11)`). Restricted to `q ≡ 3 (mod 4)`, the list is exactly `{3, 7, 11}`.
* For `q = 7` and `q = 11` the index-`q` action is the action on the out-neighbourhood design of `P_q`: the Fano plane, resp. the biplane. The octahedral stabiliser `S_4` of a Fano point is the full symmetry group of the tetrahedron, whose six edges plus one formal point give the Fano plane (session23's tetrahedral bridge).

**One structural explanation.** Extra design symmetry should mean extra balance in the 5- and 6-vertex subtournament statistics that drive `H`.

* All doubly regular tournaments of a given order share their 3- and 4-vertex statistics.
* The Fano and biplane symmetry (`PSL(2,7)`, `PSL(2,11)` acting on the out-neighbourhood design) is a candidate reason why `P_7` and `P_11` sit at the top. Its absence for `q ≥ 19` is a candidate reason why the locally transitive rotational circulant overtakes.
* This explanation is a HEURISTIC; no inequality is proved.

## Cheapest decisive tests

* `H(P_27)` (Paley over `F_27`; prime power), `H(P_31)` and `H(P_43)` against the rotational circulant. These need an orbit-reduced or meet-in-the-middle path count; the present subset DP stops at `n = 23`.
* A local search at `n = 13, 15` for non-circulant tournaments beating the rotational circulant. This would refine (1) to "Paley is maximal iff …" without any circulant restriction.
