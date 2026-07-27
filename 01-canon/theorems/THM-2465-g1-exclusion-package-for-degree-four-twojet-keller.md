---
id: THM-2465
title: "The G1 exclusion package: degree-four 2-jet Keller maps are pinned"
status: >
  PROVED (universal layer: no 2-jet Keller map has A != 0 with
  A x B == 0; B != 0 always; universal z-rationality -- for A != 0
  the generic fiber injects into the (x,y)-plane, so field degree =
  degree of a plane 0-cycle and the z-quadratic never doubles it; at
  degree 4 a wild cover can have monodromy D4, A4, or S4 (regular
  C4/V4 would be Galois and hence invertible), and at least one of
  x,y is still primitive by the corrected subgroup argument; every
  quadratic Keller map is injective, so any G1
  witness has total degree >= 3) + PROVED (stratum exclusions,
  QQ-exact Groebner with cross-prime adversarial re-derivation:
  point-cap b == 0 face empty; line-cap doubly-degenerate face
  empty; conic-cap empty of Keller maps of ANY field degree in four
  boxes through (A=2,B<=2,C<=3); box (1,2,3) r2a stratum forces
  field degree exactly 1; the (1,2,3) cusp engine needs deg B >= 3)
  + PROVED (transfer/hardness floor: field degrees 1 and 3 are
  realized on both point- and line-cap strata by explicit elementary
  extensions of the THM-1310 wild map, so an unconditional degree-4
  exclusion on those strata implies the degree-4 case of the
  order-{1,3} conjecture) + PROVED (tower constraints: purity forces
  a nonempty Jelonek set and the resolvent cubic cover to ramify
  over it; S4 witnesses need odd Jelonek valuation of Delta_4,
  A4 witnesses a square Delta_4 with cyclic cubic layer; D4 witnesses
  have a 1+2 resolvent algebra and a proper quadratic intermediate; a
  THM-1310-as-resolvent witness must satisfy the five conditions
  (N1)-(N5)). VERDICT: no witness found; G1 remains OPEN with
  detection floors recorded per box. Nothing here resolves G1,
  the order-{1,3} conjecture, or any JC.
source: kind-pasteur-2026-07-26-S134
depends_on:
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
related:
  - THM-1310-conic-pair-fibers-and-design-equations
  - THM-1340-engine-trichotomy-zaffine-keller
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/jacobian_g1_stratum_exclusions_kps_S134.py
output: 05-knowledge/results/jacobian_g1_stratum_exclusions_kps_S134.out
companions: >
  jacobian_g1_degree_arithmetic_kps_S134, jacobian_g1_strata_kps_S134,
  jacobian_g1_line_strata_kps_S134,
  jacobian_g1_staircase_replication_kps_S134,
  jacobian_g1_conic_cap_hunt_kps_S134 (+ master out),
  jacobian_g1_conic_kernel_structure_kps_S134,
  jacobian_g1_gb_wildcert_kps_S134 (GB wildcert boxes: PENDING,
  no verdict is claimed from them)
hash_basis: working-tree bytes (LF); per-file hashes in INDEX entry
---

# THM-2465 -- what a degree-four 2-jet Keller map must be

**PROVED + VERDICT OPEN** as itemized in the status. Produced by a
six-agent workflow (theory strata, degree arithmetic, staircase and
conic-cap hunts, adversarial verification, synthesis); every
load-bearing computational claim was independently re-derived over
QQ and fresh primes by the adversarial pass, with zero refutations.

## 1. Universal layer (PROVED)

For `F = A z^2 + B z + C` Keller: (U0) `B != 0` (else `D0 == 0`);
(U1) `A != 0` with `A x B == 0` is impossible (`B = rho A` makes
`det J` vanish along `z = -rho/2`); (U2) with `A != 0`, the frame
identity `z (A x B) = A x (P - C)` (THM-2455) makes `z` rational
over `C(P)(x,y)` on the generic fiber, which therefore injects into
the `(x,y)`-plane: **field degree = degree of the plane 0-cycle
`{g1 = 0} intersect {lambda = mu^2}`** -- the z-quadratic shape
never doubles the degree; (U3) a wild degree-four cover has monodromy
`D4`, `A4`, or `S4`.  Indeed the regular transitive groups `C4,V4`
would make the original extension Galois and hence invertible by Campbell.
In the `S4/A4` lanes the point stabilizer is maximal.  In the `D4` lane
the quartic root field has exactly one proper intermediate `K_2`; since the
frame identity gives `K_4=K(x,y)`, if neither coordinate were primitive then
both `K(x)` and `K(y)` would be contained in `K_2` (with `K` allowed), forcing
`K(x,y) subseteq K_2`, a contradiction.  Thus in all three surviving lanes
**at least one of
`x,y` is a primitive element with a quartic minimal polynomial**, but the
old no-intermediate-field rationale was false.

The elementary quartic

```text
f(T)=T^4-2,                 R(W)=W(W^2+8)                 (U4)
```

is the sharp field-theory hostile.  Its Galois group is `D4`;
`Q(2^(1/4))` contains the proper intermediate `Q(sqrt(2))`, while
the nontrivial resolvent factor gives `Q(sqrt(-2))`.  Thus even the
two visible quadratic sidecars need not coincide.  This is not a
Keller example; it refutes only the discarded Smith/maximality shortcut.
Every quadratic Keller map is
injective (midpoint argument), so **any G1 witness has total degree
at least 3**; the staircase boxes have even BKK ceilings
`(d+2)^3 - d^3 = 8, 26, 56, 98` with the excess over the realized
degree dropping by exactly 2 per rank-<=1 point of the fiber system.

## 2. Stratum exclusions (PROVED, QQ-exact)

- Point-cap (`A = f v0`) with `b == 0`: empty of Keller maps
  (`det J = jac(C1, C2) (2 f z + B3)`).
- Line-cap doubly-degenerate face (`B3 == 0` and `j(Abar, Bbar) == 0`):
  empty (collapses into U1).
- Conic-cap `A`: empty of Keller maps of ANY field degree in the
  boxes `(A=C1, B<=1, C<=2)`, `(A=C2, B<=1, C<=2)`,
  `(A=C1, B<=2, C<=2)`, `(A=C1, B<=2, C<=3)` (normal forms
  `C1 = (x^2, xy, y^2)`, `C2 = (x^2, x, 1)` exhaust the cap modulo
  gauge); plus all 47 enumerated B-rays of the minimal quartic box
  and their gauge orbits, and the (2,4,4) kernel directions. The
  exact kernel law `ker D4|_{C1} = {hA + fA_x + gA_y} + {linear B}`
  is part of the record.
- Box (1,2,3): the `r2a` stratum forces field degree exactly 1 (the
  Euler-eigenvalue kill `D4 = 0 <=> B3 = 0`, then `C3` constant,
  then a polynomially z-solvable shape which is an automorphism);
  the cusp engine needs `deg B >= 3`.

## 3. Realized degrees and the hardness floor (PROVED)

`W1 = (F1, F2, F3 + F1^2)` (point-cap) and
`W2 = (F1 + F2 F3, F2 + F3^2, F3)` (line-cap), both built on the
THM-1310 wild map, are 2-jet Keller maps of field degree 3; degree 2
is impossible because a quadratic extension is Galois and Campbell's
Galois-Keller theorem forces invertibility. Hence field degrees `{1, 3}` are realized on
both residual strata, and **an unconditional degree-4 exclusion
there would prove the degree-4 case of the z-affine order-{1,3}
conjecture** -- no cheap all-strata resolution exists. Conversely
this makes the G1 question a sharp, well-posed door: the residual
arena is the point/line-cap z-affine quartic-engine problem plus the
named pockets with floors (r2b, r1, staircase branches, conic-cap
off-ray kernel families, deg >= 5 tails, and the pending GB boxes).

## 4. Tower constraints on any witness (PROVED)

Purity (Zariski-Nagata) plus `pi_1^et(A^3) = 1` re-proves that a
wild Keller map has nonempty Jelonek set `A_F`, and forces the
quartic field, every nontrivial normalized resolvent factor, and the
Galois closure to be unramified in codimension one off `A_F`; every
nontrivial connected factor must ramify over some Jelonek component.
`S4` witnesses have odd `v_D(Delta_4)` along a Jelonek component (the
HYP-9027 shape); `A4` witnesses have square `Delta_4` and a cyclic
cubic layer with order-3 inertia.  In a `D4` lane the cubic resolvent
has a rational section plus a quadratic factor, and the quartic root
field has a proper quadratic intermediate which need not be that
resolvent quadratic, as (U4) shows.  This lane remains open. A witness
carrying THM-1310's map as its resolvent
layer must satisfy: (N1) `G = S4`; (N2) `Delta_4 == -L o iota` mod
squares; (N3) its Jelonek set contains a copy of `{L = 0}`;
(N4) `v(Delta_4)` odd along it; (N5) its Galois closure contains
`K3(sqrt(-L))` as the degree-6 layer.

## 5. Detection floors and pendings (honest negatives)

Sampled-only (floors, not proofs): `r2b` (3/3 inconsistent), `r1`
(3/3 degenerate), staircase branches (8/8), conic-cap B off the 47
rays inside the 13/15-dimensional kernel families, `deg C >= 5`,
`deg B >= 5`, `deg A >= 3` caps, and the `(C2, B<=2, C<=2)` joint
box (watchdog kill, NO VERDICT). The GB wildcert boxes
(`(1,1,1)+fiber4`, `(1,2,2)+fiber4`, `(2,2,2)+alpha3`) are built
but unfinished: **no verdict is claimed from them**; the decisive
next computation is finishing `(1,2,2)+fiber4` -- `GB = [1]` there
would be the first stratum-free G1 exclusion in a box not killed by
the quadratic-injectivity argument (BKK ceiling 19 >= 4).

## 6. Reproduction

```bash
python 04-computation/jacobian_g1_stratum_exclusions_kps_S134.py   # 28/28 asserts
python 04-computation/jacobian_g1_degree_arithmetic_kps_S134.py
python 04-computation/jacobian_g1_conic_cap_hunt_kps_S134.py       # long
python 04-computation/jacobian_g1_gb_wildcert_kps_S134.py box122   # pending box
```

Outputs in 05-knowledge/results with matching names; the stored
outs are the workflow transcripts whose load-bearing claims were
adversarially re-derived cross-prime in-session.
