---
id: HYP-9353
title: "The converse-Gluck exceptions: a tournament T has no vertex subset with trivial stabilizer in its group of automorphisms and anti-automorphisms (equivalently it needs three colours to be distinguished up to converse) if and only if T = X[C_3^[k]] with X one of C_3, the circulant Z_5{1,2}, the Paley tournament P_7, and k >= 0: the three prime circulants singled out by Gluck's Lemma 6, with towers of directed triangles substituted"
status: >
  PROVED (Session24, cba114930c, independently audited), relative to the cited Gluck and Feit-Thompson theorems:
  the structural argument in 05-knowledge/results/converse_exception_classification_20261021/session24-structure/STRUCTURE.md
  closes the 'only if' direction. It was OPEN at the S49 checkpoint; the S49 note's pre-audit sentence asserting
  completeness was unsupported (MISTAKE-650).
  Inherited (THM-4743, PROVED):
  - every X[C_3^[k]] is exceptional;
  - the primitive case (exactly C_3, Z_5{1,2}, P_7) and Paley tournaments (exactly q = 3, 7);
  - the orbit-pairing lemma, and that exceptional tournaments are strong (1b, from the S49 audit);
  - the substitution rule for A prime and B exceptional.
  FINITE-EXACT (S49 and its audit, the pre-proof evidence):
  - all tournaments with n <= 9: exceptions C_3, Z_5{1,2}, P_7, C_3[C_3];
  - all self-converse tournaments with n <= 11: no others;
  - all circulant tournaments on Z_n, n odd <= 15: exactly the family members;
  - the 110 vertex-transitive 21-tournaments: only P_7[C_3];
  - all Cayley tournaments on Z_3^3, Z_9 x Z_3 and Z_27: only C_3^[3];
  - Z_5{1,2}[C_3] is exceptional and C_3[Z_5{1,2}] is not, at n = 15;
  - unequal substitutions and stacks tried are not exceptional.
  Predictions 1-4 below were made before these tests and all passed.
source: opus-2026-10-21-S49, 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md (section 3); proof by Session24
related:
  - 01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md
  - 05-knowledge/results/converse_exception_classification_20261021/session24-structure/STRUCTURE.md (the proof)
  - 05-knowledge/results/signed_tournament_palette_rigidity_20261021/README.md (concurrent: three colours always suffice)
  - 01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md (path-rigid tournaments; a different classification)
scripts:
  - 04-computation/experiments/seven_point_spheres_20261021/converse_gluck_tournaments.py (+ .out)
  - 04-computation/experiments/seven_point_spheres_20261021/audit_S49/tourn_sc_census.py, tourn_vt.py, tourn_checks.py (+ .out)
  - 05-knowledge/results/converse_exception_classification_20261021/session24-structure/verify_structure.py (+ receipts)
---

# HYP-9353 — the converse-Gluck exceptions

**Closed by Session24.** The [independently audited structural classification](../results/converse_exception_classification_20261021/session24-structure/STRUCTURE.md) proves the statement for every finite tournament of order at least two.
* Ordinary vertex-transitivity excludes unequal canonical blocks and non-strong composites.
* Minimal signed blocks are triangles, and the full ternary-kernel argument proves that they are modules.
* Induction then gives exactly the stated family.

The S49 rationale and searches below are kept for provenance. They are not open obligations.

## Statement

* Call a tournament *exceptional* if no vertex subset has trivial stabilizer in `Aut±(T)`, the automorphisms together with the anti-automorphisms.
* By the concurrent session23b result, three colours always suffice. So exceptional means "needs exactly three colours up to converse".
* **Theorem (formerly the conjecture).** `T` is exceptional iff `T ≅ X[C_3^[k]]` with `X ∈ {C_3, Z_5{1,2}, P_7}` and `k ≥ 0`.
* Here `C_3^[k] = C_3[C_3[…]]` is the `k`-fold tower and `X[C_3^[0]] = X`.
* The orders are `n = 3^(k+1)`, `5·3^k`, `7·3^k`. The smallest members are `C_3 (3), Z_5{1,2} (5), P_7 (7), C_3[C_3] (9), Z_5{1,2}[C_3] (15), P_7[C_3] (21), C_3^[3] (27)`.

## Rationale at the S49 checkpoint

* **Reduction.** By THM-4743 (1), `T` is exceptional iff its anti-automorphisms fix every regular `Aut(T)`-orbit on subsets. This is a rigidity condition, so large symmetric pieces are needed.
* **Strongness is forced.** By THM-4743 (1b), exceptional tournaments are strong. Session24 proves more: they are vertex-transitive.
* **The primitive case is closed.** By Gluck's Theorem 1, the only primitive solvable groups without regular subsets and of order `≡ 2 (mod 4)` are `S_3`, `D_10` and `F_42`.
* **Equal substitutions into a prime tournament are closed.** `A[B]` with `A` prime and `B` exceptional is exceptional iff `D(Aut±(A)) > r(B)`.
  * In the data, exceptional `B` have `r(B) ∈ {2, 4, 6}`, and `r(B) = 2` occurs only for the triangle towers.
  * Exceptional prime `A` need three colours. So only towers can be substituted into `X`.
* **Routes that were open, both closed by Session24:**
  * substitutions with unequal blocks `A[B_1, …, B_m]` into a prime `A`;
  * imprimitive `Aut±` whose blocks are not modules.

## Pre-registered searches

All four were named before testing, and all passed. They are now controls for the theorem.

1. All self-converse tournaments on 10 and 11 vertices: no exceptions (audit census, `2^25` and `2^30` labelled cases).
2. Vertex-transitive tournaments on 21 vertices: only `P_7[C_3]` (88 + 22 tournaments).
3. Cayley tournaments on `Z_3^3`, `Z_9 × Z_3` and `Z_27`: only `C_3^[3]`.
4. Unequal substitutions such as `C_3[C_3, C_3, Z_5{1,2}]`, and stacks: not exceptional.

Session24's palette recursion gives the substitution counts in closed form. For `B` strong and self-converse with `q` regular types, `f` of them fixed by the converse:
* `C_3[B]` has `q' = (q³ − q)/3` and `f' = f(q − 1)`;
* if `B` is exceptional, `C_3[B]` has `C(q, 3)` signed-regular orbits;
* so the triangle is the only exceptional block for `C_3[·]`, which explains the 15,000 regular subsets of `C_3[Z_5{1,2}]`.

## Notes

* THM-2450's finite-scope grammar is historical.
  * [THM-2477](../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md) and the [all-order bypass proof](../results/merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified) classify the path-rigid tournaments (`H = |Aut|`) as the ordered stacks of `{1, C_3}`.
  * THM-2454 gives the strong pure-blue atoms `{1, C_3, T_5}`.
  * Neither is this classification: `C_3[C_3]` is exceptional with `H/|Aut| = 39`.
* `n = 7` and `n = 21` occur here as orders of exceptional tournaments. They are not linked to the Hamiltonian-path gaps 7 and 21 of THM-1370: no mechanism is claimed (NUMEROLOGY guard). S50 examines the spectrum side.
