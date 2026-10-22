---
id: THM-4752
title: "Tournament scores in II_{25,1}: for a 25-vertex tournament T the lift v_T = (s | 70) has norm -2 c3(T) and is null exactly for the transitive tournaments, which are the Leech cusps (Conway's Weyl vector rho = (0,1,...,24 | 70), the cannonball 0^2 + ... + 24^2 = 70^2); the backward length of every ranking pi is the Weyl-vector pairing -v_T . rho_pi, minimised by the score rankings; and the score hyperplane is the level -25 of the A_24 deep-hole cusp nu = (1/2^25 | 5/2) = (sum of the 25 Leech roots e_i - e_(i+1), r*), so every tournament has a class in the Niemeier lattice N(A_24) = nu^perp/nu, namely its centred score vector s - 12 in A_24"
status: >
  PROVED (elementary; statement (c) uses Conway's theorem that the Leech roots of II_{25,1} are indexed by the Leech
  lattice with r_lambda . r_mu = 2 - |lambda - mu|^2/2, and the Niemeier-Venkov theorem that a Niemeier root system has
  one Coxeter number and rank 0 or 24). FINITE-EXACT checks: lorentzian_scores.py. KNOWN ingredients: the c3 formula
  (Kendall-Babington Smith 1940), Landau's theorem (1953), score rankings minimise backward length (Kano-Sakamoto 1983),
  the cannonball (Watson 1918), the Weyl vector and Leech roots (Conway-Sloane 1982, Conway 1983, Borcherds-Conway-
  Queen-Sloane 1984), cusps and deep holes (Conway-Parker-Sloane 1982, Borcherds 1985), forest counts (Stanley 1980).
  Statement (a) was announced in MSG-1913 (a08f25ecbf) and independently verified, with a hostile collision, by
  Session27 (05-knowledge/results/session27_repunits_lattices_20261022/LORENTZIAN_SCORE_BRIDGE.md). No prior instance
  of the dictionary itself was found (repo canon, web and MathOverflow searches, 2026-10-21/22). Session opus-2026-10-21-S52.
session: opus-2026-10-21-S52
source: 05-knowledge/results/lattices_modular_repunits_20261021.md
scripts:
  - 04-computation/experiments/lattices_repunits_20261021/lorentzian_scores.py (+ .out)
  - 04-computation/experiments/lattices_repunits_20261021/lattice_cusps.py (+ .out) (the II_{9,1}, II_{17,1} contrast)
related:
  - THM-868 (n = 8: centred scores in E8; n = 9: the A_8 slice)
  - THM-869 / THM-872 (n = 16: D16+, n = 24: D24+; the named-open question which Niemeier lattice the scores select)
  - THM-484 (Remark 1: the cannonball recorded as context, its meeting with tournaments left open)
  - Session27 LORENTZIAN_SCORE_BRIDGE.md (statement (a), the padded collision: equal v_T, H = 657 vs 643)
---

# THM-4752 — tournament scores in II_{25,1}

## Setting

II_{25,1} is the even unimodular Lorentzian lattice of vectors x = (x_0, ..., x_24 | x_25) with all coordinates in Z or
all in Z + 1/2 and x_0 + ... + x_24 - x_25 even, with the form x . y = sum_{i<25} x_i y_i - x_25 y_25. Put

    rho = (0, 1, ..., 24 | 70),   nu = (1/2, ..., 1/2 | 5/2),   r_i = e_i - e_(i+1) (0 <= i <= 23),
    r* = (-1/2, 1/2, ..., 1/2, 3/2 | 5/2)   (23 middle coordinates 1/2).

For a tournament T on vertices 0..24 with scores s_v and c3(T) cyclic triangles, let v_T = (s | 70). A ranking is a
bijection pi from the vertices to 0..24 (24 = top); rho_pi = (pi | 70); the backward length is
BL_T(pi) = sum of pi(w) - pi(u) over the arcs u -> w with pi(u) < pi(w) (Kano-Sakamoto).

## Statement

**(a) Norm.** v_T lies in II_{25,1} and v_T . v_T = -2 c3(T). It is null exactly when T is transitive, and then
v_T = sigma(rho) for a coordinate permutation sigma.

**(b) Backward length is the Weyl-vector pairing.** For every ranking pi, BL_T(pi) = -v_T . rho_pi. It is minimised
exactly by the rankings sorted like the scores, and the minimum is

    U(T) = sum_k (s_(0) + ... + s_(k-1) - C(k,2)) = c3(T) + |s_sorted - (0, 1, ..., 24)|^2 / 2,

the total Landau slack. So U(T) is the Lorentzian distance from v_T to the nearest of the 25! Leech cusps sigma(rho).

**(c) The A_24 deep hole.** r_0, ..., r_23 and r* are Leech roots (r . r = 2, r . rho = -1). Their Gram matrix is the
affine Cartan matrix of A~_24 (a 25-cycle with all marks 1), and r_0 + ... + r_23 + r* = nu. Hence:
- nu is a primitive null vector in the closure of the Leech chamber;
- the Leech roots orthogonal to nu are exactly these 25;
- nu . rho = -25 = -h(A_24);
- nu^perp / nu is the Niemeier lattice N(A_24).

In Conway's coordinates nu corresponds to a deep hole of the Leech lattice of type A_24, whose 25 vertices are the
24 adjacent transpositions of the score ranks together with r*.

**(d) Every tournament has an A_24 class.** v_T . nu = -25 for every T: the score hyperplane is the level -25 slice of
the A_24 cusp. The vector u_T = v_T - 10 (0^25 | 1) = (s | 60) lies in nu^perp, and u_T - 24 nu = (s - 12 | 0). So
the class of T in N(A_24) is its centred score vector d_T = s - 12. This lies in the root lattice A_24 itself (glue
class 0), with |d_T|^2 = 1300 - 2 c3(T), and d_T = 0 exactly when T is regular.

**(e) The light cone.** The reflections in the roots +-(e_i - e_j) of nu^perp generate W(A_24) = S_25, the vertex
relabellings. By Landau's theorem the score vectors are exactly the lattice points of the permutohedron
conv(S_25 rho). Its 25! vertices are the null vectors sigma(rho), and every other point is timelike.

**(f) Why 25.** 0^2 + ... + (n-1)^2 is a square only for n = 1, 2, 25 (Watson), and II_{n,1} exists only for
n = 1 (mod 8). So n = 25 is the only order at which the transitive score vector is a null vector of an even unimodular
Lorentzian lattice.

## Proof

(a) The coordinates are integers and 300 - 70 = 230 is even. Every transitive triple has one source, so
c3 = C(25,3) - sum_v C(s_v, 2). With sum s_v = 300 this gives sum s_v^2 = 4900 - 2 c3 = 70^2 - 2 c3. A tournament is
transitive iff c3 = 0, and then its scores are 0..24 in some order.

(b) Count each pair {u, w} with pi(u) < pi(w) once. Its arc adds pi(w) to sum_v pi(v) s_v when w -> u, and
pi(u) = pi(w) - (pi(w) - pi(u)) when u -> w. Summing over all pairs,

    sum_v pi(v) s_v = sum_{j=0}^{24} j * j - BL_T(pi) = 4900 - BL_T(pi),

and 4900 - sum_v pi(v) s_v = -v_T . rho_pi. By the rearrangement inequality, sum_v pi(v) s_v is maximal exactly when
pi is sorted like s; this is the Kano-Sakamoto theorem for tournaments. At such pi,

    -v_T . rho_pi = sum_j j (j - s_(j)) = sum_j (s_(j) - j)(24 - j) = sum_k Landau slack_k,

using sum_j (s_(j) - j) = 0. Expanding |s_sorted - iota|^2 = sum s^2 - 2 sum j s_(j) + sum j^2 and inserting (a)
gives U = c3 + |s_sorted - iota|^2/2. (This identity holds at every order n, with 4900 replaced by sum_{j<n} j^2; the
Lorentzian reading needs n = 25.)

(c) The norms, the pairings with rho and nu, the Gram matrix and the sum are direct computations. nu is primitive
because nu/2 has coordinates 1/4.

By Conway's theorem, the simple roots of II_{25,1} for the chamber at rho are the roots with r . rho = -1. They are
r_lambda (lambda in the Leech lattice), with r_lambda . r_mu = 2 - |lambda - mu|^2/2 <= 0 for lambda != mu, since
the Leech minimum is 4. For any Leech root r outside the 25,

    r . nu = sum (r . r_i) + r . r* <= 0,

so nu lies in the closure of the chamber. If r . nu = 0, then every term vanishes and r is orthogonal to all 25
roots. The Ã24 Gram matrix has rank 24 and its kernel vector gives the null vector nu != 0, so the 25 roots span
nu^perp. Then r would be a multiple of nu, which is impossible for a root. So the Leech roots orthogonal to nu are
exactly the 25.

nu^perp / nu is positive definite, even, unimodular, of rank 24, and contains the images of +-(e_i - e_j): 600
roots forming A_24. By Niemeier-Venkov its root system has one Coxeter number and full rank. A_24 is the only
rank-24 system containing A_24 with that property, so nu^perp / nu = N(A_24).

In coordinates Lambda + U with rho = (0, 0, 1), a primitive null vector with z . rho = -h is z = (h c, h, h|c|^2/2)
and z . r_lambda = h(1 - |lambda - c|^2/2). For nu this gives h = 25, and c is a point at distance >= sqrt 2 from
the Leech lattice, at distance exactly sqrt 2 from 25 points forming A~_24: the deep hole of type A_24
(Conway-Parker-Sloane, Borcherds).

(d) v_T . nu = 300/2 - 70 * 5/2 = -25. u_T . nu = 150 - 150 = 0, and u_T is in the lattice since 300 - 60 is even.
24 nu = (12^25 | 60) gives u_T - 24 nu = (s - 12 | 0). Then |s - 12|^2 = sum s^2 - 24 * 300 + 25 * 144 =
sum s^2 - 3600 = 1300 - 2 c3.

(e) Landau's theorem in majorization form says s is a score vector iff its sorted partial sums dominate
C(k, 2), with equality at k = 25. That is membership of the permutohedron of (0, ..., 24) together with
integrality. By (a), v_T is null exactly at the vertices.

(f) Watson (1918); even unimodular II_{p,q} exists iff p - q = 0 (mod 8).

## Scope

- The Lorentzian point does not determine H: Session27 pads two 8-vertex tournaments to identical v_T with
  H = 657 and 643.
- The Leech quotient rho^perp / rho receives only the transitive tournaments. For cyclic T, v_T^perp / v_T is not a
  null quotient (Session27). The A_24 cusp of (d) is the canonical quotient defined for every T, but it remembers only
  the scores.
- The extra root r* changes the time coordinate, so its reflection is not a tournament operation. Nothing here
  concerns the gaps 7 and 21 or repunit H-values.

## Relation to the score-lattice ladder

THM-868 places centred scores at n = 8 in E8 and at n = 9 in A_8. THM-869 and THM-872 place them in D16+ at n = 16 and
in D24+ at n = 24, and leave open which Niemeier lattice the score slice selects. At n = 25 the question has a
canonical answer: the scores are the level -h slice of the A_24 cusp nu, every tournament's class is its centred
score vector in N(A_24), and the Leech lattice enters as the cusp of the transitive tournaments.

The lower ranks are Vinberg-reflective with no such tournament reading (lattice_cusps.py):
- II_{9,1} has 10 simple roots (E10) and one cusp, E~8 -> E8, with rho . z = -30 = -h(E8) and rho^2 = -1240;
- II_{17,1} has 19 simple roots and two cusps, E~8 + E~8 and D~16, which are the two heterotic lattices;
- 0^2 + ... + 8^2 and 0^2 + ... + 16^2 are not squares, so there is no tournament analogue at n = 9 or 17.
