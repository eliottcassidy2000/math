---
id: THM-4752
title: "Tournament scores in II_{25,1}: for a 25-vertex tournament T the lift v_T = (s | 70) has norm -2 c3(T) and is null exactly for the transitive tournaments, the Leech cusps (Conway's Weyl vector rho = (0,1,...,24 | 70), cannonball 70^2); the backward length of every ranking pi is the pairing -v_T . rho_pi, minimised by the score rankings, with the minimising cusp also the nearest (squared distance 2(U - c3)); and the score set lies on the level -25 slice of the A_24 cusp nu = (1/2^25 | 5/2) of Borcherds' table (KNOWN: nu^perp/nu = N(A_24), nu . rho = -h), so every tournament has the class s - 12 in N(A_24)"
status: >
  PROVED, elementary: the tournament dictionary (a), (b), (d), (e) and (f).
  KNOWN, reproduced FINITE-EXACT: the lattice facts in (c).
  - Borcherds' table of norm-0 vectors of II_{25,1} (https://math.berkeley.edu/~reb/lattices/norm0) lists 2nu for A_24
    with height 25 and 25 / 5 / 5 simple roots at inner products 0 / -1 / -2.
  - Borcherds' thesis (arXiv math/9911195), Lemmas 2.2.4, 2.4.2 and 2.5.1, gives rho_N^2 = 2h(h+1), z . rho = -h and
    rho = (rho_N, h, h+1).
  - Conway 1983 (Leech roots); Conway-Parker-Sloane 1982 (deep holes).
  New here is only the tournament reading of these lattice facts.
  Classical tournament ingredients:
  - c3 = C(n,3) - sum C(s_i,2) (Kendall-Babington Smith 1940); the score variance form is THM-785 (1) and THM-4739;
  - Landau (1953) = the permutohedron (THM-1825-q-redei-census, THM-207);
  - tournaments as root choices in A_{n-1}, with sum x^score = s_delta (camion_busch §6.1; Gessel 1979);
  - score rankings minimise backward length (Kano-Sakamoto, Networks 13 (1983) 143-151);
  - the cannonball (Watson 1918).
  The closed form U = c3 + |s_sorted - iota|^2 / 2 is elementary and possibly folklore. Monsuur, EJOR 164 (2005), is
  unread (paywalled).
  Prior repo sightings:
  - ogg_triangular_chowla_20260927.md named (0,...,24; 70) "the Leech vector in II_(25,1)", typed NUMEROLOGY;
  - statement (a) was announced in MSG-1913 (a08f25ecbf) and independently verified by Session27, which also gives the
    hostile collision (LORENTZIAN_SCORE_BRIDGE.md).
  Independently audited (audit_S52/REPORT.md; MISTAKE-661). Session opus-2026-10-21-S52.
session: opus-2026-10-21-S52
source: 05-knowledge/results/lattices_modular_repunits_20261021.md
scripts:
  - 04-computation/experiments/lattices_repunits_20261021/lorentzian_scores.py (+ .out)
  - 04-computation/experiments/lattices_repunits_20261021/lattice_cusps.py (+ .out) (the II_{9,1}, II_{17,1} contrast)
  - 04-computation/experiments/lattices_repunits_20261021/audit_S52/a_lorentzian_audit.py (+ .out) (independent)
related:
  - THM-868 (n = 8: centred scores in E8; n = 9: the A_8 slice)
  - THM-869 / THM-872 (n = 16: D16+; n = 24: D24+; their named-open n = 24 Niemeier question stays open)
  - THM-785 / THM-4739 (score variance n(n^2-1)/12 - 2 c3, the norm in (a) and (d))
  - Session27 LORENTZIAN_SCORE_BRIDGE.md (statement (a); the padded collision: equal v_T, H = 657 vs 643)
  - HYP-9363 (Leech-root counts at the A_24 cusp)
---

# THM-4752 — tournament scores in II_{25,1}

## Setting

II_{25,1} is the even unimodular Lorentzian lattice. Its vectors are x = (x_0, ..., x_24 | x_25) with all coordinates
in Z or all in Z + 1/2 and x_0 + ... + x_24 - x_25 even, under the form x . y = sum_{i<25} x_i y_i - x_25 y_25.

Fix

    rho = (0, 1, ..., 24 | 70),   nu = (1/2, ..., 1/2 | 5/2),   r_i = e_i - e_(i+1) (0 <= i <= 23),
    r* = (-1/2, 1/2, ..., 1/2, 3/2 | 5/2)   (23 middle coordinates equal to 1/2).

For a tournament T on vertices 0..24, write s_v for the scores, c3(T) for the number of cyclic triangles, and
v_T = (s | 70).

A ranking is a bijection pi from the vertices to 0..24, with 24 the top; set rho_pi = (pi | 70). The backward length
BL_T(pi) is the sum of pi(w) - pi(u) over the arcs u -> w with pi(u) < pi(w) (Kano-Sakamoto).

## Statement

**(a) Norm.** v_T lies in II_{25,1} and v_T . v_T = -2 c3(T). It is null exactly when T is transitive, and then
v_T = sigma(rho) for a coordinate permutation sigma. These 25! vectors are cusps whose null quotient is the Leech
lattice (Conway; BCQS 1984).

**(b) Backward length is the Weyl-vector pairing.** For every ranking pi,

    BL_T(pi) = -v_T . rho_pi.

The rankings sorted like the scores, and only those, minimise it. The minimum is

    U(T) = sum_k (s_(0) + ... + s_(k-1) - C(k,2)) = c3(T) + |s_sorted - (0, ..., 24)|^2 / 2,

the total Landau slack. So U(T) = min over sigma of (-v_T . sigma(rho)): the horospherical height of v_T over the
Leech cusps.

The minimising cusp is also the nearest one. The difference v_T - rho_pi = (s - pi | 0) is spacelike with
|v_T - rho_pi|^2 = 2(BL_T(pi) - c3), so the squared distance to the nearest cusp is |s_sorted - iota|^2 = 2(U - c3).

**(c) The A_24 cusp (lattice side KNOWN).** r_0, ..., r_23 and r* are Leech roots (r . r = 2, r . rho = -1). Their
Gram matrix is the affine Cartan matrix of A~_24, a 25-cycle with all marks 1, and

    r_0 + ... + r_23 + r* = nu.

Consequently:
- nu is a primitive, future-pointing null vector in the closure of the Leech chamber;
- the Leech roots orthogonal to nu are exactly these 25;
- nu . rho = -25 = -h(A_24);
- nu^perp / nu is the Niemeier lattice N(A_24).

This is the entry 2nu of Borcherds' norm-0 table, given in Conway's standard coordinates, which are read here as
tournament coordinates. The 25 vertices of the A_24 deep hole are the 24 adjacent transpositions of the score ranks
together with r*.

nu is canonical: the S_25-fixed future null rays are spanned by nu and by nu^- = (-1^25 | 5), and only nu lies in the
closed Leech chamber (nu^- . (-e_0 - e_1) = 2 > 0, while -e_0 - e_1 is a Leech root).

**(d) Every tournament has an A_24 class.** v_T . nu = -25 for every T, so the score vectors lie on the level -25
slice of the A_24 cusp, intersected with x_25 = 70.

The vector u_T = v_T - 10 (0^25 | 1) = (s | 60) lies in nu^perp, and u_T - 24 nu = (s - 12 | 0). So the class of T in
N(A_24) is its centred score vector d_T = s - 12. This vector:
- lies in the root lattice A_24 itself (glue class 0);
- has |d_T|^2 = 1300 - 2 c3(T);
- is zero exactly when T is regular.

Any S_25-invariant w with w . nu != 0 used in place of (0^25 | 1) gives the same class. The transitive class
(12, 11, ..., -12) is the Weyl vector rho_{A_24}, with rho_{A_24}^2 = 1300 = 2h(h+1) (Venkov). This is Borcherds'
decomposition rho = (rho_N, h, h+1) read in tournament coordinates.

*Glue classes.* N(A_24) = { x - (sum x / 25) 1 : x in Z^25, 5 | sum x }. So the centred out-degree vector of any
25-vertex digraph with m arcs, 5 | m, lies in N(A_24), in the glue class [-m mod 25]. Tournaments (m = 300) fill
class 0 only.

**(e) The light cone.** The reflections in the roots +-(e_i - e_j) of nu^perp generate W(A_24) = S_25, the vertex
relabellings. By Landau's theorem the score vectors are exactly the lattice points of the permutohedron conv(S_25 rho).
Its 25! vertices are the null cusps sigma(rho); every other point is timelike.

**(f) Why 25.** 0^2 + ... + (n-1)^2 is a nonzero square only for n = 2 and 25 (Watson). II_{n,1} exists only for
n = 1 (mod 8). So n = 25 is the only order at which the transitive score vector is a nonzero null vector of an even
unimodular Lorentzian lattice.

## Proof

**(a)** The coordinates are integers, and 300 - 70 is even. Every transitive triple has exactly one source, so
c3 = C(25,3) - sum_v C(s_v, 2). With sum s_v = 300 this gives sum s_v^2 = 4900 - 2 c3. A tournament is transitive
iff c3 = 0, and then its scores are 0..24 in some order.

**(b)** Each pair {u, w} with pi(u) < pi(w) contributes pi(w) to sum_v pi(v) s_v if w -> u, and
pi(u) = pi(w) - (pi(w) - pi(u)) if u -> w. Summing over pairs,

    sum_v pi(v) s_v = sum_j j^2 - BL_T(pi) = 4900 - BL_T(pi).

Since v_T . rho_pi = sum_v pi(v) s_v - 4900, this is BL_T(pi) = -v_T . rho_pi.

By the rearrangement inequality, sum_v pi(v) s_v is maximal exactly when pi is sorted like s; this is Kano-Sakamoto's
theorem. At such a pi, using sum_j (s_(j) - j) = 0,

    -v_T . rho_pi = sum_j j (j - s_(j)) = sum_j (s_(j) - j)(24 - j) = sum_k Landau slack_k.

Expanding |s_sorted - iota|^2 and inserting (a) gives U = c3 + |s_sorted - iota|^2 / 2. For the distance, note
|(s - pi | 0)|^2 = sum s^2 - 2 sum pi(v) s_v + sum pi(v)^2 = 2(BL - c3).

These identities hold at every order n, with 4900 replaced by sum_{j<n} j^2. Only the Lorentzian reading needs n = 25.

**(c)** Norms, pairings, the Gram matrix and the sum are direct computations. nu/2 has coordinates 1/4, so nu is
primitive, and nu_25 > 0, so nu is future-pointing.

By Conway's theorem, the simple roots for the chamber at rho are the roots r_lambda with r . rho = -1, indexed by the
Leech lattice. They satisfy r_lambda . r_mu = 2 - |lambda - mu|^2 / 2 <= 0 for lambda != mu. So for a Leech root r
outside the 25,

    r . nu = sum (r . r_i) + r . r* <= 0,

and nu lies in the closed chamber.

If r . nu = 0, every term vanishes and r is orthogonal to the 25 roots. These span nu^perp over R: the Gram matrix has
rank 24, and its kernel vector sums the roots to the nonzero null vector nu. Then r would lie in R nu, which is
impossible for a root. Over Z the 25 roots span an index-5 sublattice of nu^perp modulo nu, which is the glue.

nu^perp / nu is even, unimodular, positive definite of rank 24, and contains the 600 roots +-(e_i - e_j), an A_24
system. By Niemeier-Venkov its root system has a single Coxeter number and full rank, so it is A_24 and
nu^perp / nu = N(A_24). The audit also computes this directly: A_24 plus glue [10], with glue minima 4, 6, 6, 4.

nu . rho = -25 = -h is Borcherds' Lemma 2.4.2.

**(d)** v_T . nu = 150 - 175 = -25. u_T . nu = 150 - 150 = 0, and 300 - 60 is even. Since 24 nu = (12^25 | 60),
u_T - 24 nu = (s - 12 | 0), and |s - 12|^2 = sum s^2 - 3600 = 1300 - 2 c3.

For a different invariant w = (a 1 | b) with w . nu != 0, the projected vector differs from u_T by a vector of the
form (c 1 | c'), and that vector lies in nu^perp only along nu. So the class is unchanged.

The glue description holds because A_24^* / A_24 = Z/25, and N(A_24) is the preimage of the order-5 subgroup, i.e. the
vectors all of whose coordinates are congruent to i/25 mod 1 with 5 | i.

**(e)** Landau's theorem in majorization form is membership of the permutohedron of (0, ..., 24) plus integrality.
Timelikeness follows from (a).

**(f)** Watson (1918); and II_{p,q} exists iff p - q = 0 (mod 8).

## Scope

- The Lorentzian point does not determine H. Session27 pads two 8-vertex tournaments to the identical v_T, with
  H = 657 and 643.
- The Leech quotient rho^perp / rho receives only the transitive tournaments. For cyclic T, v_T^perp / v_T is not a
  null quotient (Session27). The A_24 cusp of (d) is the canonical quotient defined for every T, and it remembers only
  the scores.
- The reflection in r* changes the time coordinate, so it is not a tournament operation.
- Nothing here concerns the gaps 7 and 21 or repunit H-values.

## Relation to the score-lattice ladder

THM-868, THM-869 and THM-872 place centred scores as follows:
- n = 8: E8;
- odd n: A_{n-1} (THM-868 §4 at n = 9);
- n = 16: D16+;
- n = 24: D24+.

At n = 25 the centred score vector lies in A_24, the odd-n case. The new step is the canonical identification of the
ambient Niemeier lattice: it is nu^perp / nu for the unique S_25-fixed cusp in the closed Leech chamber, with the Leech
lattice entering as the cusp of the transitive tournaments.

The n = 24 question named in THM-869 and THM-872 (which Niemeier lattice, if any, the scores select) is not answered
here.

In lower rank there is no tournament analogue: 0^2 + ... + 8^2 and 0^2 + ... + 16^2 are not squares. Both lattices
are Vinberg-reflective (lattice_cusps.py):
- II_{9,1}: 10 simple roots (E10), one cusp E~8 -> E8, with rho . z = -30 = -h(E8) and rho^2 = -1240;
- II_{17,1}: 19 simple roots, two cusps E~8 + E~8 and D~16, the two heterotic lattices.
