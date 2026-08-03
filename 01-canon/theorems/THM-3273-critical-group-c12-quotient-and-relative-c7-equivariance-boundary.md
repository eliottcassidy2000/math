---
id: THM-3273
title: "Critical-group C12 quotient and relative C7 equivariance boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  THM-3260's twelve-vertex core has cyclic critical group C_74748 and
  therefore a canonical Hall-{2,3} quotient C12.  Its Abel--Jacobi vertex
  image occupies only six classes, with multiplicities 3,3,3,1,1,1, so it
  is not a C12 torsor; all twelve pair differences nevertheless occur.  Core
  directed edges miss exactly phases 4 and 8, and the delayed pair
  (11,17),(11,21) supplies exactly 8 and 4, raising the augmentation sampler
  from rank ten to eleven.
  The delayed relative H1 is saturated Z^6, but its intrinsic C2 signature
  (5,1) is incompatible with both trivial and affine-reflection actions on
  Aug(Q[C7]).  The full critical group has no 7- or 13-primary torsion.
source: root/2026-08-03
audit: >
  The assertion-independent exact companion pins THM-3254 and THM-3260,
  parses the literal covering-pair banks, reconstructs all three reduced
  Laplacians, computes their exact Smith forms, verifies a primitive adjugate
  coordinate and every C12 vertex/pair/edge label, computes the saturated
  relative boundary and explicit basis, checks the core and repaired
  augmentation samplers, checks both C2 signatures, and separates absence of
  C13 torsion from the order-13 automorphisms of the C131 primary factor.
  Normal, optimized and stored output agree byte-for-byte. An independent
  reconstruction recovered all three Smith groups from determinantal
  divisors, the primitive adjugate chart, every vertex/pair/edge difference,
  the automorphism multiplier, the full saturated relative kernel and both
  possible affine-C7 involution signatures. It also closed the companion's
  nonblocking lattice-index audit gap directly from the kernel equations.
  A second immutable audit rederived the delayed increments 8 and 4, the
  rank-ten kernel and rank-eleven repaired sampler, the C2 anti-line, the
  C131 null chart and the exhaustive order-13 affine nonpreservation scan.
depends_on:
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
related:
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
script: 04-computation/gmc_critical_group_phase_carrier_thm3273.py
output: 05-knowledge/results/gmc_critical_group_phase_carrier_thm3273.out
script_sha256: ab50ab5f2c5ea0fd9f7c25504f38bc3982fbfa9cc8fa6892b9ee3dfdf13b738b
output_sha256: da1d102a73ae19b3d28b54dbe6a8243b1ef587c0782bda7172fb17e6f019d7d0
hash_basis: LF-normalized bytes
---

# THM-3273 -- critical-group C12 quotient and relative C7 equivariance boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3260 shows that the first-link core has Betti number eleven and the delayed
relative layer has rank six, matching the dimensions of the `C_12` and `C_7`
augmentation modules without supplying either cyclic action.  The graph
Jacobian is the next intrinsic carrier to test: unlike Betti rank, it remembers
the integral incidence lattice and can contain canonical torsion quotients.

It produces a real but incomplete repair.  A canonical `C_12` quotient exists,
and its vertex differences see every phase, but its vertex image is not a
torsor.  On the delayed side even this torsion mechanism is absent, and the
intrinsic involution has the wrong rational signature for an affine `C_7`
carrier.

## 1. Exact critical groups

For a connected graph `X`, write

```text
Jac(X)=coker(L_X'),                                     (1)
```

where `L_X'` is any reduced integral Laplacian.  Let `G_0`, `G` and `H` be
THM-3260's leafless first-link core, its one-leaf reset graph, and the full
31-edge covering graph.  Exact Smith reduction gives

```text
Jac(G_0) = C_74748,
Jac(G)   = C_74748,
Jac(H)   = C_4 direct_sum C_12 direct_sum C_113184.     (2)
```

The leaf bridge explains the equality of the first two groups.  Their orders
factor as

```text
74748  = 2^2 * 3 * 6229 = 12 * 6229,
5432832= |Jac(H)| = 2^9 * 3^4 * 131.                  (3)
```

In particular none of the three critical groups has 7- or 13-primary torsion.
Since `gcd(12,6229)=1`, the Hall-`{2,3}` part of `Jac(G_0)` is canonically
cyclic of order twelve.  Equivalently,

```text
J_12 := Jac(G_0)/12 Jac(G_0) is isomorphic to C_12.    (4)
```

The quotient `(4)` is intrinsic as an abstract group.  An isomorphism with a
numbered copy of `Z/12`, and hence a generator, is not intrinsic.

## 2. Abel--Jacobi landing: six absolute classes, twelve differences

Root `G_0` at row 22 and map a vertex `v` to the divisor class
`[v-22]`.  A primitive row of the reduced-Laplacian adjugate is

```text
(12951,7527,2466,8883,7593,6369,14988,34229,4932,5961,9313), (5)
```

in nonroot order

```text
(2,3,7,10,11,13,16,17,18,19,21).
```

Its gcd with 74748 is one, and multiplying `(5)` by the reduced Laplacian
gives 74748 times the row-17 coordinate functional.  It therefore identifies
`Jac(G_0)` with `Z/74748`.  Reduction modulo twelve gives the following chart
of `(4)`:

```text
row:    2  3  7 10 11 13 16 17 18 19 21 22
phase:  3  3  6  3  9  9  0  5  0  9  1  0.          (6)
```

Thus the twelve vertices occupy only

```text
{0,1,3,5,6,9}, with multiplicities {3,1,3,1,1,3}.     (7)
```

Changing the root translates `(6)`; changing the isomorphism in `(4)`
multiplies it by a unit.  Both operations preserve the six-class image and
the multiplicity multiset `3,3,3,1,1,1`.  Consequently no such choice turns
the vertices into a `C_12` torsor.

The strongest survivor is genuinely useful:

```text
{phase(v)-phase(w): v,w in V(G_0)} = Z/12.             (8)
```

So the intrinsic quotient contains every relative phase.  Local first-link
edges do not quite suffice: after allowing both orientations their difference
set is

```text
(Z/12) \ {4,8}.                                        (9)
```

This gap is repaired exactly by two delayed edges.  Let
`A=Aug(Q[J_12])`, regarded as the zero-sum rational functions on `J_12`, and
for an oriented edge set `E` define the phase sampler

```text
R_E:A --> Q^E,
(R_E f)(u,v)=f(phase(v)-phase(u)).                      (9a)
```

For the directed core edges, `(9)` gives

```text
rank R_core=10,
ker R_core=Q(delta_4-delta_8).                          (9b)
```

The delayed edges

```text
(11,17) has increment 8,
(11,21) has increment 4.                               (9c)
```

Adjoining their two orientations makes `(9a)` injective, of rank eleven.
Thus the delayed layer contains precisely the missing one-dimensional sampler,
not merely another six-dimensional coincidence.

The nontrivial graph automorphism `(17 21)` acts on `J_12` as multiplication
by 5.  It swaps phases 4 and 8 and swaps the two edges in `(9c)`.  Hence the
anti-invariant line `delta_4-delta_8` is compatible with the graph involution.
This is a functorial action on the critical quotient, not a faithful `C_12`
action on the graph.

## 3. Saturated delayed layer and the wrong C2 signature

Order the delayed relative edges by

```text
e1=(2,7), e2=(3,9), e3=(7,14), e4=(11,17),
e5=(11,21), e6=(12,13), e7=(12,19), e8=(14,19).       (10)
```

Relative to `G`, the only new vertices are 9 and 12, and the relative boundary
matrix is

```text
             e1 e2 e3 e4 e5 e6 e7 e8
[9]           0  1  0  0  0  0  0  0
[12]          0  0  0  0  0 -1 -1  0.                (11)
```

Its Smith form is `diag(1,1)`.  Hence the relative homology is saturated and
has the explicit integral basis

```text
H_1(H,G;Z)=<e1,e3,e4,e5,e8,e6-e7> = Z^6.              (12)
```

The intrinsic edge-coloured involution swaps `e4,e5` and fixes the other four
basis vectors.  Over `Q`, its `(+,−)` signature on `(12)` is therefore

```text
(5,1).                                                  (13)
```

Its anti-invariant line is `Q(e4-e5)`, exactly the delayed detector that
separates the two missing phase atoms in `(9b)--(9c)`.

An affine involution of a `C_7` torsor is either trivial or a reflection.  On
`Aug(Q[C_7])` these give signatures `(6,0)` and `(3,3)`, respectively.  Thus
the relative layer is not equivariantly identifiable with the `C_7`
augmentation module even after one weakens the missing order-seven action to
the graph's actual `C_2` symmetry.  The absence of 7-primary torsion in `(2)`
also rules out repairing it through a critical-group `C_7` quotient.

## 4. Carrier theorem and boundary

Equations `(4)`, `(8)` and `(9a)--(9c)` are the positive result: the response
graph itself contains a canonical order-twelve *relative-phase* carrier, and
the core plus two delayed edges faithfully samples its full eleven-dimensional
augmentation space.  This is stronger than the bare dimension match in
THM-3260 and requires no external finite group object.  It can type pairwise
phase differences and every rational augmentation word.

Equations `(7)` and `(9)` locate the remaining defect.  The carrier does not
label the twelve vertices bijectively, does not choose an absolute origin or
generator, and misses two phases on individual first-link edges.  Therefore it
cannot replace THM-3260's external cyclic vertex label or THM-3255's absolute
phase marker.  The sampler is an injective coefficient-space map, not an edge
action, a cycle projection or a physical response.  Pairwise difference
completeness is not absolute reference.

The delayed rank-six space is still further away: it is saturated free
homology with neither `C_7` critical torsion nor the correct involution type.
The same computation excludes a graph-Jacobian `C_13` torsion subgroup or
quotient as a mod-13 response label.  It does **not** exclude all order-13
automorphisms.  The characteristic `C_131` primary factor of `Jac(H)` has
automorphism group `C_130`, hence a unique subgroup of order thirteen.
The rooted full-graph Abel--Jacobi image modulo 131 occupies eleven classes,
and exact enumeration shows that no nonidentity element of this order-13
subgroup, even followed by a translation, preserves that image.  Thus an
abstract `C_13` automorphism survives, but it does not act on the vertex-image
chart and is not a `C_13` phase label.

No map from response rows to owner phases, positive Gaussian-moment response,
NC2 theorem, row exclusion, or `LRC(14)` decrement follows.  No claim is made
that a critical-group character is an existing physical observable.

## 5. Exact companion

Run

```text
python3 04-computation/gmc_critical_group_phase_carrier_thm3273.py
python3 -O 04-computation/gmc_critical_group_phase_carrier_thm3273.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
exact integer Smith reduction and adjugate arithmetic, contains no floating
point or randomness, and raises explicitly on every failed certificate.

QED.
