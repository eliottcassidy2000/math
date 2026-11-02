---
id: THM-4662
title: "The genus-one fish curve C3 (11 z0^2 = D(1-D), 11 z1^2 = D(4-D), 11 z2^2 = (1-D)(D+3)) has no point over the real cyclic quartic field F_4 = Q(tau), tau^2 = (55 + 11 sqrt5)/2: C3 is the two-covering of E: y^2 = x(x-9)(x+7) with class (33, 33), E has rank 0 over Q(sqrt5) and its twist by tau^2 has rank 2 there, so E(F_4) has rank 2 and its two-descent image is spanned by four explicit classes not containing (33, 33). Consequently the Hochberg-O'Donnell fish and the Exoo-Ismailescu graph EI21 have no faithful unit-distance placement in Q(zeta55), nor in any conjugation-stable odd-degree extension of Q(zeta5, sqrt(-11)), the cubic carriers Q(zeta5, sqrt(-11), theta) and Q(zeta55, theta) included"
status: >
  PROVED (computer-assisted): the two-Selmer bounds rest on exact rational arithmetic and local images verified to have
  their known orders (4 at odd places, 16 at the prime 2 of Q(sqrt5), 2 at real places), reproduced with two sampling
  seeds and checked on controls with known ranks (congruent-number curves n = 1, 2, 3, 5, 6, 7 over Q; y^2 = x^3 - x
  and y^2 = x^3 - 121x over Q(sqrt5)); rank E(Q(sqrt5)) = 0 also agrees with PARI's ellrank over Q for E and its
  quintic twist (fish_marked_elliptic_cover_20261013). The two generators of the twist are exact points.
  The consequences use THM-4661 statements 4(b), 4(c) and the fish note's equations (1) and reduction (5).
  Independently audited 2026-10-16 (04-computation/experiments/s38_audit_six_colour_20261016/REPORT.md): an independent
  two-descent reproduces both Selmer groups, every local image is maximal isotropic for the local Tate pairing, and a
  residue-symbol test at 59 split primes confirms independence and non-membership of (33, 33); statements 1-5
  CONFIRMED. clock_duality_oct15 found and repaired a false-negative square oracle in the original script (5 = (sqrt5)^2
  was rejected) and replayed the descent with the same conclusion (cubic_descent_audit_20261016e); MISTAKE-622.
session: opus-2026-10-16-S38
source: 05-knowledge/results/six_colour_odd_carrier_twenty_vertices_20261016.md (section 4)
scripts:
  - 04-computation/experiments/two_descent_q_sqrt5_fish_curve_20261016.py (+ .out, .json)
related:
  - 05-knowledge/results/fish_marked_elliptic_cover_20261013.md (C3, its elliptic quotient E, rank E(Q(sqrt5)) = 0, the reduction C3(L^+) != empty <=> C3(F_4) != empty; left C3(F_4) OPEN)
  - THM-4661 (statements 4(b), 4(c): fish and EI21 placements in odd-degree carriers over K give points of C3 over F_4)
  - HYP-9247 (the fish and EI21 were its two OPEN placement tests)
  - 05-knowledge/results/ei19_fish_placement_q_zeta55_20261012.md, 05-knowledge/results/ei21_kappa_q_20261012.md
---

# THM-4662 — the fish curve has no point over the real quartic field

**2026-10-16 focused computational repair (clock_duality_oct15).** The original
`is_square_K` omitted purely irrational roots of rational inputs, incorrectly
rejecting `5 = (sqrt5)^2`. The implementation now handles that branch and uses
exact rational comparisons at real places. The focused replay
`04-computation/experiments/cubic_descent_audit_20261016e.py` recovers the same
Selmer sizes 4 and 16, checks the two displayed points, and verifies the
four-dimensional image and target exclusion by both corrected field square
tests and a separate F2 calculation. Normal and optimized executions agree.
The original receipts preserve provenance; the repaired replay is
`05-knowledge/results/cubic_descent_audit_20261016e.json`. This is not a new
independent implementation of the entire local descent. See also B6--B8 of
`05-knowledge/results/cubic_adelic_20261016e.md` for the completed local-solubility
argument, index two, and an explicit even-degree reroute for Gamma_9.

## Setting

`K0 = Q(sqrt5)`, `phi = (1 + sqrt5)/2`, `delta = (55 + 11 sqrt5)/2 = 11 sqrt5 phi`, `F_4 = K0(tau)`, `tau^2 = delta`: the real
cyclic quartic subfield of `K = Q(zeta5, sqrt(-11))` (`tau = sqrt(-11)(zeta5 - zeta5^-1)`). The fish curve `C3` is the
smooth projective genus-one curve defined by

    11 z0^2 = D(1-D),   11 z1^2 = D(4-D),   11 z2^2 = (1-D)(D+3),

`E: y^2 = x(x-9)(x+7)`, and `E^delta: Y^2 = X(X - 9 delta)(X + 7 delta)` its quadratic twist by `delta` over `K0`.
The descent map is `P -> (x - 0, x - 9)` modulo squares (with the usual values at the two-torsion points).

## Statement

1. `C3` is isomorphic over `Q` to the two-covering of `E` with class `(33, 33)`.
2. `rank E(K0) = 0` and `rank E^delta(K0) = 2`; the points
   `P1 = (484, 1573 - 3751 sqrt5)` and `P2 = ((-65 + 19 sqrt5)/2, 160 + 400 sqrt5)` of `E^delta(K0)` have descent classes
   `(1, 11)` and `(-sqrt5 phi, -(4 - sqrt5))`, independent modulo the torsion classes. Hence `rank E(F_4) = 2`.
3. The image of `E(F_4)/2E(F_4)` in `(F_4^x / F_4^x2)^2` is the span of `(-7, -1)`, `(-3, -3)` (the rational torsion),
   `(1, 11)` and `(-sqrt5 phi, -(4 - sqrt5))`, and `(33, 33)` is not in it. Therefore **`C3(F_4)` is empty**.
4. Consequently the Hochberg-O'Donnell fish (23 vertices) and the Exoo-Ismailescu graph EI21 (21 vertices) have no
   injective unit-edge placement in `L = Q(zeta55)`, and none in any conjugation-stable field `M ⊇ K` with `[M:K]` odd,
   in particular in `K' = K(theta)` and `L' = L(theta)`, `2 theta^3 + 2 theta^2 - 2 theta - 1 = 0`.
5. **(A forbidden nine-vertex configuration.)** Let `Gamma_9` be the graph on `h, a, b, c, d, w_ab, w_ac, w_bd, w_cd` with
   the 14 edges `h-a, h-b, h-c, h-d`, `a-w_ab, b-w_ab, a-w_ac, c-w_ac, b-w_bd, d-w_bd, c-w_cd, d-w_cd`, `w_ab-w_cd`,
   `w_ac-w_bd` (four rhombi at the hub `h` and two closing edges). No unit-distance graph containing `Gamma_9` (as an
   injective, not necessarily induced subgraph) has a placement in `Q(zeta55)` or in any field `M` as in 4. `Gamma_9`
   occurs in EI21 (hub 15) and EI19 (hub 12); it is realizable in the plane (a one-parameter family).

## Proof

**1.** For `D != 0`, `x = 12/D - 3` gives `x = 3 D(4-D)/D^2 = 33 (z1/D)^2` and `x - 9 = 12 D(1-D)/D^2 = 33 (2 z0/D)^2`,
while `x + 7 = 4 D(D+3)/D^2` and `D(D+3) = (11 z0 z2/(1-D))^2`. The assignments `w1 = z1/D`, `w2 = 2 z0/D`,
`w3 = 22 z0 z2/((1-D) D)` map `C3` birationally onto the covering `x = 33 w1^2`, `x - 9 = 33 w2^2`, `x + 7 = w3^2`, with
inverse `D = 12/(x + 3)`, `z1 = w1 D`, `z0 = w2 D/2`. A two-covering has a point over a field `F` exactly when its class
lies in the descent image of `E(F)`.

**2.** `F_4 = K0(sqrt delta)` gives `E(F_4) ⊗ Q ≅ E(K0) ⊗ Q ⊕ E^delta(K0) ⊗ Q`. The full two-torsion descent over `K0`
(class number one; `S` = primes over 2, 3, 7 for `E`, and over 2, 3, 7, 5, 11 for `E^delta`; S-unit generators
`-1, phi, 2, 3, 7, sqrt5, 4 + sqrt5, 4 - sqrt5`) gives two-Selmer groups of orders 4 and 16. So `rank E(K0) = 0` and
`rank E^delta(K0) <= 2`. The script finds exact points in Selmer classes that, together with the torsion classes, span
all 16 elements; `P1` and `P2` are two of them (their coordinates satisfy the equation exactly), so the rank is 2 and
`Sha(E^delta/K0)[2] = 0`. Hence `rank E(F_4) = 2`.

**3.** `E(F_4)[2]` is all of `E[2]`, so `dim E(F_4)/2E(F_4) = 2 + 2 = 4`. A point `(X, Y)` of `E^delta(K0)` gives the point
`(X/delta, Y/(delta tau))` of `E(F_4)` with descent class `(X, X - 9 delta)`, because `delta` is a square in `F_4`. An
element of `K0` is a square in `F_4` iff it or its quotient by `delta` is a square in `K0`. Write classes in the basis
`-1, phi, 2, 3, 7, sqrt5, pi = 4 + sqrt5, pi' = 4 - sqrt5` of `K0(S,2)`, modulo `delta ≡ phi sqrt5 pi pi'`. The four
classes `(-7, -1)`, `(-3, -3)`, `(1, 11)`, `(-sqrt5 phi, -pi')` are independent there (the `7`- and `3`-coordinates and the
sign force the coefficients of the first, second and fourth to vanish, then the `pi`-coordinate of the second component
forces the third). So they span the image. If `(33, 33)` were in the span, the first component would force the
coefficients `(0, 1, *, 1)` and one factor `delta`. The second component then needs one `delta` factor for its
`pi`-coordinate and none for its `pi'`-coordinate, a contradiction. The script confirms this by exact square tests.

**4.** A fish placement in `M` gives a point of `C3` over `M^+` (fish note, equations (1)), and so does an EI21 placement
(THM-4661 4(c)). Since `[M^+ : F_4] = [M : K]` is odd (respectively `[L^+ : F_4] = 5`), the degree-eight-divisor argument of
THM-4661 4(b) gives a point of `C3` over `F_4`, which does not exist. **5.** The proof of THM-4661 4(c) uses only the
nine vertices and fourteen edges of `Gamma_9`: the four rhombi at the hub (non-degenerate because the vertices are
distinct) give the four sums, the closing edges give `|P|^2 + |Q|^2 = 1` and `P ⊥ Q`, and `0 < X = |a - d|^2 < 1` because
`a != d` and `b != c`. So a placement of `Gamma_9` gives a point of `C3` over `M^+`, and 4 applies. ∎

## Remarks

- This settles the two OPEN placement tests of HYP-9247 (the fish, EI21): of the four smallest known triangle-free
  four-chromatic unit-distance graphs, EI19, the fish and EI21 are now excluded from `Q(zeta55)` and from every
  conjugation-stable odd-degree carrier over `K` (EI19 from every odd-degree extension of `K` or `L`). Only EI17
  remains, excluded from `L` only numerically to coordinate height `1e63`.
- The obstruction is global: `E(F_4)` has rank 2, but `(33, 33)` is not a descent class of any point.
  The initial hand-check suggestion of everywhere local solubility is now proved in B6 of
  `05-knowledge/results/cubic_adelic_20261016e.md`: good odd places use the full unramified local
  two-descent image, and the places over 2, 3 and 7 use points at infinity. Thus `C3` represents a
  nonzero class of `Sha(E/F_4)[2]`; its index is exactly two. This strengthens the remark, not the
  dependencies of the global nonexistence proof above.
- The fish note's bounded searches (36,518 quadratic and 1,728 quartic parameters) were consistent with this; the
  proof replaces them.
