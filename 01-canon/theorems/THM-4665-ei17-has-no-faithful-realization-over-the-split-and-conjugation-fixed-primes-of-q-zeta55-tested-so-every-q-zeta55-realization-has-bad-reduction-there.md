---
id: THM-4665
title: "EI17 over the primes of Q(zeta55): the 17-vertex Exoo-Ismailescu graph has no faithful realization in the split model over F_p for 26 of the 32 primes p = 1 (mod 55) below 12100, and no embedding at all in the finite unit-distance graph G_p = Cay(F_(p^2), mu_(p+1)) for each of the 20 primes p = 54 (mod 55) below 6000; hence every realization of EI17 with all directions in Q(zeta55) has bad reduction at every prime of Q(zeta55) above each of these 46 primes"
status: >
  PROVED: the reduction lemma (statement 3), the exhaustiveness of the elimination (parallelogram forcing over any field
  of characteristic not 2 by the two-circle bound, of which THM-4663 (2) is the G_q case; each degenerate branch forces a
  vertex collision in every model, enumerated by the split script and skipped by the twisted one), the gauge reduction.
  Independently audited (04-computation/experiments/s39_audit_ei17_20261017/REPORT.md): every count reproduced by a
  lemma-free embedding search with zero disagreements; one interpretive sentence corrected (MISTAKE-625). FINITE-EXACT: the counts (exact arithmetic modulo p; every faithful solution
  re-checked on all 31 edges), including the controls (faithful points do occur at other primes). NOT CLAIMED: an
  exclusion of EI17 from Q(zeta55) (a realization of large height could have bad reduction at all these primes; S24's LLL
  bound excludes heights below about 1e63). HYP-9247 is updated, not settled.
session: opus-2026-10-17-S39
source: 05-knowledge/results/ei17_finite_fields_20261017.md
scripts:
  - 04-computation/experiments/ei17_finite_field_realizations_20261017.py (split model; + .out)
  - 04-computation/experiments/ei17_twisted_unit_graph_embeddings_20261017.py (twisted model = embeddings in G_p; + .out)
related:
  - 05-knowledge/results/ei17_placement_q_zeta55_20261012.md (S24: direction structure, 3000-digit realization, LLL to height 1e63, the refuted degree-20 certificate)
  - HYP-9247 (every finite unit-distance graph in Q(zeta55) is three-colourable), THM-4661, THM-4662 (EI19, fish, EI21 excluded)
  - THM-4628 (conjugation-fixed primes of Q(zeta55) have residue degree 2 or 10; conjugation = Frobenius on the residue field), THM-4663 (the graphs G_p, codegree lemma)
---

# THM-4665 — EI17 over the split and conjugation-fixed primes of `Q(zeta55)`

## Setting

EI17: 17 vertices, 31 edges, triangle-free, four-chromatic, minimally rigid (Exoo-Ismailescu 2016; House of Graphs
51375). `L = Q(zeta55)`. In every faithful realization (17 distinct points, unit edges) over any field of characteristic
not 2, the nine unit 4-cycles are parallelograms (in a 4-cycle `a b c d` with distinct vertices, `w = c - a` is nonzero,
and two unit circles with distinct centres meet in at most two points for any nondegenerate binary form in characteristic
not 2, since subtracting the norm equations leaves a line; THM-4663 (2) is the case `G_q`. The set `{b - a, d - a}` is
closed under `x -> w - x` and a fixed point would force `b = d`, so `c - b = d - a`: opposite steps are negatives), so the 31 edges carry 13 direction classes with six
integer relations (S24). In the gauge `v17 = 0`, direction of `(11, 17)` equal to 1, a realization is the six free
directions `U0, U1, U2, U3, U4, U6` with

    U7 = U1 + U2 - U3 - 1,      U8 = -U0 - U2 + U4 + 1,      U9 = -U0 + U1 - U2 + U3 + U6 + 1,
    U10 = 2U2 - U3 - U6 - 1,    U11 = U1 + U2 + U3 - U4 + U6 - 1,    U12 = -U1 - U2 - U6 + 2

all of norm one, and the 17 vertices are the integer combinations listed in the scripts.

* **Split model over `F_p`:** directions `x ∈ F_p^×` with conjugate `1/x`, points `(z, z̄) ∈ F_p^2`, norm `z z̄`.
* **Twisted model over `F_p`:** directions in `μ_(p+1) ⊂ F_(p^2)`, points in `F_(p^2)`: embeddings of EI17 in
  `G_p = Cay(F_(p^2), μ_(p+1))`; every embedding is moved into the gauge by the arc-transitive affine automorphisms.

## Statements

1. **(Split primes, FINITE-EXACT.)** For the 26 primes `p ≡ 1 (mod 55)` below 12100 in the list
   `331, 661, 881, 991, 1321, 1871, 2311, 2531, 2971, 3301, 3631, 3851, 5171, 5281, 5501, 6271, 6491, 7151, 7591, 8581,
   9241, 9461, 9791, 9901, 10891, 11551`, the split model has no faithful realization. At the other six primes
   `p ≡ 1 (mod 55)` below 12100 (`2861, 3191, 4621, 4951, 7481, 10781`) it has 8 or 16.
2. **(Conjugation-fixed primes, FINITE-EXACT.)** For each of the 20 primes `p ≡ 54 (mod 55)` below 6000
   (`109, 439, 659, 769, 1319, 1429, 1759, 1979, 2089, 2309, 2749, 2969, 3079, 3299, 3739, 4289, 4729, 5059, 5279, 5939`),
   EI17 is not a subgraph of `G_p`.
3. **(Reduction, PROVED.)** If all directions of a realization of EI17 lie in `L`, then for every prime `P` of `L` above
   a prime `p ≡ 1 (mod 55)` (residue field `F_p`; complex conjugation is `σ_(-1)`, so the conjugate of a unit direction
   `u` is `1/u` in `L`) the reduction modulo `P` is a point of the split model, and for every prime `P` above
   `p ≡ 54 (mod 55)` (residue field `F_(p^2)`, on which conjugation acts as Frobenius, THM-4628) it is a point of the
   twisted model; the point is faithful unless some direction is not a `P`-unit or two vertices are congruent modulo
   `P` and modulo `σ_(-1) P`. Hence, by 1 and 2, every such realization has bad reduction at every prime of `L` above each
   of the listed 46 primes. At the conjugation-fixed primes every unit direction is automatically a `P`-unit
   (`v_P(u) = -v_P(σ_(-1) u) = -v_P(u)`), so there bad reduction means two vertices congruent modulo `P`.
4. **(Controls, FINITE-EXACT.)** Faithful points do occur: in the split model at `p = 103, 139, 281, 313, 401, 499, 521,
   541, 631, 641, 677` (below 700) and at the six primes of statement 1; in the twisted model (embeddings of EI17 in
   `G_p`) at `p = 127, 149, 227, 277, 283, 311, 317, 349, 383, 397, 443, 463, 547, 587` (below 600) and, of the 38 primes `p ≡ 10 (mod 11)` below 2600, exactly at `2441` (16) and `2551` (8) (independent
   audit run, `s39_audit_ei17_20261017/a2_twisted_10mod11_table.out`; the split model has no point at any of these 38).
   All counts are multiples of 4; all but one (`p = 383`, twisted, 4) are multiples of 8.

## Proof

The elimination is exhaustive: choose `U1, U2`; `|U7| = 1` and `|U12| = 1` are quadratics in `U3` and `U6` (the codegree
lemma: `c̄ x^2 - N(c) x + c = 0` for `N(c - x) = 1`); check `|U10| = 1`; `|U11| = 1` and `|U8| = 1` are quadratics in `U4`
and `U0`; check `|U9| = 1`; keep solutions with 17 distinct vertices and re-check all 31 edges. A quadratic degenerates
only when its constant vector `c` vanishes: `U1 + U2 = 1`, `U1 + U2 = 2`, `U1 + U2 + U3 + U6 = 1` or `U4 = U2 - 1`, which
force `v5 = v11`, `v7 = v10`, `v2 = v11` or `v1 = v6` respectively (sympy identities), so these branches contain no
faithful realization (the split script enumerates them anyway). Statement 3: reduction of the algebraic identities
`u σ_(-1)(u) = 1`, `U_k = L_k(U)`, `L_k(U) σ_(-1)(L_k(U)) = 1` modulo `P`.

## Reading

S24 left EI17 the only one of the four smallest triangle-free four-chromatic unit graphs not excluded from `L`, with a
numerical exclusion to height `1e63` and an elimination that did not finish. The finite-field count replaces the
elimination by about `16 p^2` operations per prime and turns the question into one about reduction: a realization in
`L` would have to be degenerate (a non-unit direction or a vertex collision) at all forty primes above each of 26
split primes and at the primes above 20 conjugation-fixed primes. A zero count is the typical outcome at any prime (108 of 119
split and 91 of 105 twisted control primes); the class `1 (mod 55)` is not depleted (6 of 32, binomial tail 0.07) and the
class `54 (mod 55)` only mildly (0 of 20, probability 0.06). What is specific to the classes `±1 (mod 55)` is the
inference: only there does a zero count force bad reduction of a hypothetical realization in `L`. EI17 does embed in many
`G_p` (`G_127`, ...). (Corrected after the audit; the earlier sentence called the obstruction itself specific to these
classes, MISTAKE-625.)

## Open

An unconditional exclusion (a bound on the bad primes of a hypothetical realization, e.g. from the discriminant of the
realization scheme over `Z`); the residue-degree-ten primes of `L` (twisted model over `F_(p^5)`); the Galois structure
of the faithful realization scheme. The third involution is identified by the audit: the block re-placement
`ι: (U0, U4) -> (U9, U11)`, `v4 -> 2m - v4`, `v6 -> 2m - v6`, with `m` the common midpoint of `v1 v14` and `v2 v11`
(`v1 + v14 = v2 + v11` identically). `ι` is fixed-point free and the mirror is fixed-point free for `p > 17`, so faithful
counts are even for `p > 17` (PROVED); `<mirror, automorphism, ι> = (Z/2)^3` acts freely at every tested prime with points
except `p = 383` (twisted), where mirror∘automorphism fixes all four. That `ι` always preserves faithfulness is open.
The reduction mechanism is also THM-4558 Lemma R (KNOWN: Madore, arXiv:1509.07023, Prop. 3.2).
