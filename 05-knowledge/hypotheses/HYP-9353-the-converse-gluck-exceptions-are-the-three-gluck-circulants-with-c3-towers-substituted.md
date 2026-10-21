---
id: HYP-9353
title: "The converse-Gluck exceptions: a tournament T has no vertex subset with trivial stabilizer in its group of automorphisms and anti-automorphisms (equivalently it needs three colours to be distinguished up to converse) if and only if T = X[C_3^[k]] with X one of C_3, the circulant Z_5{1,2}, the Paley tournament P_7, and k >= 0: the three prime circulants singled out by Gluck's Lemma 6, with towers of directed triangles substituted"
status: >
  OPEN (the 'only if' direction). PROVED around it (THM-4743): every X[C_3^[k]] is exceptional; the primitive case
  (exactly C_3, Z_5{1,2}, P_7); Paley tournaments (exactly q = 3, 7); the orbit-pairing lemma; the substitution rule for
  A prime, B strong. FINITE-EXACT: all tournaments with n <= 9 (exceptions C_3, Z_5{1,2}, P_7, C_3[C_3]); all circulant
  tournaments on Z_n, n odd <= 15 (exceptions are exactly the family members C_3, Z_5{1,2}, P_7, C_3[C_3], Z_5{1,2}[C_3]);
  Z_5{1,2}[C_3] exceptional and C_3[Z_5{1,2}] not, at n = 15.
source: opus-2026-10-21-S49, 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md (section 3)
related:
  - 01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md
  - 05-knowledge/results/signed_tournament_palette_rigidity_20261021/README.md (concurrent: three colours always suffice)
  - 01-canon/theorems/THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers.md
scripts:
  - 04-computation/experiments/seven_point_spheres_20261021/converse_gluck_tournaments.py (+ .out)
---

# HYP-9353 — the converse-Gluck exceptions

## Statement

* Call a tournament *exceptional* if no vertex subset has trivial stabilizer in `Aut±(T)`, the automorphisms together with the anti-automorphisms.
* By the concurrent session23b result, three colours always suffice. So exceptional means "needs exactly three colours up to converse".
* **Conjecture.** `T` is exceptional iff `T ≅ X[C_3^[k]]` with `X ∈ {C_3, Z_5{1,2}, P_7}` and `k ≥ 0`.
* Here `C_3^[k] = C_3[C_3[…]]` is the `k`-fold tower and `X[C_3^[0]] = X`.
* The orders are `n = 3^(k+1)`, `5·3^k`, `7·3^k`. The smallest members are `C_3 (3), Z_5{1,2} (5), P_7 (7), C_3[C_3] (9), Z_5{1,2}[C_3] (15), P_7[C_3] (21), C_3^[3] (27)`.

## Why

* **Reduction.** By THM-4743 (1), `T` is exceptional iff its anti-automorphisms fix every regular `Aut(T)`-orbit on subsets. This is a rigidity condition, so large symmetric pieces are needed.
* **The primitive case is closed.** By Gluck's Theorem 1, the only primitive solvable groups without regular subsets and of order `≡ 2 (mod 4)` are `S_3`, `D_10` and `F_42`.
* **Imprimitivity through modules is closed for lexicographic products** `A[B]` with `A` prime and `B` strong. Such a product is exceptional iff `D(Aut±(A)) > r(B)`.
  * Exceptional `B` have `r(B) ∈ {2, 4, 6, …}`.
  * `r(B) = 2` occurs for the triangle towers, and only for them in the data.
  * Exceptional prime `A` need three colours. So only towers can be substituted into `X`.
* **Remaining routes to a counterexample:**
  * substitutions with unequal blocks `A[B_1, …, B_m]`;
  * transitive or stacked (non-strong) composites;
  * imprimitive `Aut±` whose blocks are not modules, e.g. Cayley tournaments on non-cyclic groups.

## What would refute it

An exceptional tournament outside the family. Concrete searches:

1. All self-converse tournaments on 10 and 11 vertices (A002785). The prediction is no exceptions, since `10` and `11` are not of the form `3^(k+1)`, `5·3^k` or `7·3^k`.
2. Vertex-transitive tournaments on 21 vertices (`04-computation/circtourn21.d6`, `othervttourn21.d6`). The prediction is that `P_7[C_3]` and its isomorphic copies are the only exceptional ones.
3. Cayley tournaments on `Z_3^3` and `Z_9 × Z_3` (`n = 27`). The prediction is exceptional iff isomorphic to `C_3^[3]`. (`Z_15` is cyclic and already covered by the circulant census.)
4. Unequal substitutions `C_3[B_1, B_2, B_3]` with exceptional `B_i` of different types, e.g. `C_3[C_3, C_3, Z_5{1,2}]`. The prediction is not exceptional, since the blocks are distinguishable.

## Notes

* The family is the maximally symmetric corner of THM-2450's grammar (`stack`, `C_3[T1, T2, T3]`) for rigid self-converse tournaments, whose leaves are the opposite extreme.
* `n = 7` and `n = 21` occur here as orders of exceptional tournaments. They are not linked to the Hamiltonian-path gaps 7 and 21 of THM-1370: no mechanism is claimed (NUMEROLOGY guard).
