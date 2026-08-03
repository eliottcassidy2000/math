---
id: THM-3245
title: "Pointed divisor median cubes, saturation-band no-go, and the z219 sorted-position boundary"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/modular-divisor-synthesis/2026-08-03
depends_on:
  - THM-3174-projected-k3-z223-terminal-descent-and-cap222
  - THM-3179-projected-k3-z222-composite-divisor-square-terminal-descent-and-cap221
  - THM-3207-projected-k3-z221-coprime-terminal-descent-and-cap220
  - THM-3218-projected-k3-z220-valuation-product-terminal-descent-and-cap219
  - THM-3230-projected-k3-z219-common-gcd-three-terminal-descent-and-cap218
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-3157-pointed-resolvent-c3-lift-edge-hexagon-and-relative-tournament
  - THM-3242-projected-k3-z217-exact-status-annihilation
script: 04-computation/lrc14_pointed_divisor_median_cubes_saturation_band_no_go_z219_supplier_support_thm3245.py
output: 05-knowledge/results/lrc14_pointed_divisor_median_cubes_saturation_band_no_go_z219_supplier_support_thm3245.out
hash_basis: LF-normalized bytes
script_sha256: 86ccb2ea0bb13f802a0e2bdbac9d2f2478ce420d5702cb3a710b8ab7ff4a9c81
output_sha256: e2551728f636fa852f036f277d3904f160b15a5ef9453ffce5cd79b96dc4f12f
semantic_sha256: 58491fb73468826a3cab7e689ff66a20e46e409f0d391c1d7881c30beaef9e88
---

# THM-3245 -- pointed divisor median cubes, saturation-band no-go, and the z219 sorted-position boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Pointed divisor intervals are median partial cubes

Let `z,L` be positive integers, put

```text
g=gcd(z,L),                         f=L/g,                 (1)
```

and suppose `f | D | L`.  Then

```text
D |-> q=D/f                                                     (2)
```

is a pointed lattice isomorphism

```text
[f,L]_(divisibility)  ~=  Div(g),       with D/L=q/g.      (3)
```

If `g=product p^(a_p)`, valuation coordinates identify `Div(g)` with

```text
product_(p|g) {0,1,...,a_p}.                              (4)
```

Meet and join are coordinatewise minimum and maximum.  Replacing an exponent
`e` by its threshold block

```text
(1_(e>=1),...,1_(e>=a_p))                                 (5)
```

embeds the Hasse graph isometrically into a cube: Hamming distance equals
valuation `l1` distance.  Coordinatewise median is the unique graph median of
each triple and remains in the monotone-block image.  Hence every pointed
divisor interval is a median partial cube.

The word *partial* is essential.  For example, `Div(4)` maps to
`{00,10,11}`; an ambient cube geodesic may leave that set.

## 2. Internal restoration is an idempotent band, not a free-factor action

For `a|g`, define the natural internal restoration

```text
J_a(q)=lcm(a,q).                                             (6)
```

These maps satisfy

```text
J_a J_b = J_lcm(a,b),          J_a^2=J_a,                   (7)
```

and commute.  Moreover, `J_a` is bijective exactly when `a=1`: for `a>1`,
`J_a(1)=J_a(a)=a`.  Thus the internal restoration monoid is a commuting
idempotent band, and every homomorphism from it to a group is trivial.

In particular, the order-two and order-three generators in

```text
PSL_2(Z)=C2*C3                                                (8)
```

cannot be obtained from nontrivial internal restorations.  The primes `2`
and `3` are special in `(8)` because of invertible external motions, not
because `Div(2)` or `Div(3)` has two or three states: both have two states,
while the three-state example `Div(4)` is a path and not a `C3` orbit.

## 3. External actions exist but require a new carrier map

The companion retains both sharp controls.

- The unpointed square underlying `Div(6)` admits regular XOR translations by
  `V4`; every nonidentity translation moves the pointed bottom and fails to
  preserve join.
- The Boolean lattice `Div(30)=B3` admits external coordinate permutations.
  A transposition and a three-cycle give the `S3` quotient of `(8)` and do
  preserve meet, join, top, and bottom.

The second control is not an LRC action.  It is an automorphism of an abstract
divisor lattice, with no supplied map on a common physical carrier, no
coefficient phase, and no owner or ancestry preservation.  THM-3157 provides
a genuine pointed `C3` elsewhere, on six edge states with its own lift data;
that structure cannot be inferred from `(3)`.

## 4. Exact Farey boundary

For a divisor cover `q -> pq` in `Div(g)`, direct reduction gives

```text
det_F(q/g,pq/g) = (p-1)g/(pq).                              (9)
```

It equals one exactly on the top `p=2` edge `g/2 -> g`.  Therefore occupied
vertices `{3,6}` at `z1=222` induce an ambient Hasse edge whose fractions
`1/2 -> 1` happen to be Farey adjacent.  No transition between those vertices
is asserted.  The phenomenon is not general: `2/20 -> 4/20` has determinant
five.  Even the incomparable `Div(6)` atoms `1/3,1/2` are Farey neighbors,
but their mediant is `2/5` while their lattice join is `1`.

Consequently `(9)` is scalar arithmetic, not THM-2056's physical flank
certificate.  That theorem still requires an ordered primitive flank, a
common owner, acute Gram data, and its defect inequality.

## 5. Exact projected-layer specialization

The companion reconstructs the pinned atlas strata at levels `222` through
`215`.  The successive gcd supports are

```text
z219: {3},
z218: {2},
z217: {7},
z216: {8,18,24,36,72},
z215: {1,5}.                                                (10)
```

Thus the adjacent `3,2` is not a structural alternating free-factor clock.
It is followed immediately by `7`, a five-stratum composite layer, and
`{1,5}`.  At `z1=218`, `g=2` on all 119 occupied rows is ambient atlas
metadata only; this theorem asserts no residual quotient or action there.

At `z1=219`, THM-3230 proves that every one of the 424 residual masks has

```text
g=3,               q=D/(L/g)=3,              q/g=D/L=1.   (11)
```

The apparent prime-three coordinate has therefore saturated to one occupied
top vertex.  It has no three-cycle to rotate.

## 6. The z219 sorted-position graph is not a supplier graph

Each stored residual denominator tuple is sorted,

```text
ds=(d_0<=d_1<=d_2<=d_3).
```

Record the *positions* at which `v_3(d_i)=v_3(L)`.  The exact census is

```text
{0}:92, {1}:136, {2}:66, {3}:5,
{0,1}:45, {0,2}:24, {1,2}:52, {1,3}:4.                    (12)
```

No set of size three occurs.  As a conditional combinatorial object, the
position-support graph has edges

```text
01, 02, 12, 13:       a 3-cycle on 0,1,2 with a leaf at 1. (13)
```

Its unweighted graph automorphism group is `C2`, generated by `0<->2`; its
weighted graph automorphism group is trivial.  Position 3 occurs with mass
`5+4=9`.  Deleting that leaf manufactures an unweighted `S3` symmetry of the
remaining 3-cycle, but discards nine actual denominator patterns.  Even after
deletion, the unequal weights leave only the identity.

This is **not** a physical supplier-label graph.  The exact screen first keeps
one labelled representative in `states[ds]["labels"]`, then stores only the
sorted multiset `ds` in the residual checkpoint.  Consequently positions
`0,1,2,3` in `(12)` are numerical order statistics, not role labels.  The
checkpoint has already quotiented away the label map needed to interpret the
graph as a carrier action.  Reconstructing genuine supplier support requires
either retaining `states[ds]["labels"]` or rerunning that labelled stage.

The scalar quotient `q`, together with sorted-position data, records top
nonemptiness and, for higher exponents, attained truncated valuation height.
It forgets supplier identities, intersections, and restoration order.  Those
are precisely the coordinates needed before a divisor shadow can be promoted
to an ordered carry or modular action.

## 7. Lawful positive target and stopping boundary

A nontrivial modular action requires at least:

1. one common physical carrier `X`;
2. bijections `S,R:X->X` preserving the carrier predicate and the actual
   coefficient/phase, owner, and ancestry data used downstream;
3. exact nontrivial relations `S^2=id` and `R^3=id`; and
4. for nonabelianity, `SR!=RS` on that same carrier.

These conditions yield a quotient of `C2*C3`; faithfulness still needs a
reduced-word or Bass--Serre normal-form sidecar.  THM-2596 identifies the
appropriate Farey transport state as an ordered basis with its Gram-owner
data.  The six-state regular `S3` action in the companion is a positive typing
control for such an external action, not evidence that the LRC quotient has
one.

## 8. Evidence and scope

The companion exhaustively verifies `(1)--(9)` on all relevant divisor
lattices; checks every triple median and every restoration composition;
reconstructs the eight exact layer censuses and row hashes; verifies the
`V4`, `S3`, two-state-prime, supplier-coupling, and cross-level hostiles; and
enumerates all conditional position-graph permutations in `(12)--(13)`.  It pins the
promoted THM-3174/3179/3207/3218/3230 evidence and checks THM-3230's exact
quotient census.  The tracked residual-mask data file is extracted from
THM-3230's exact checkpoint; this companion parses it, matches all ten banks
against THM-3230's promoted per-body counts and hashes, and only then derives
the sorted-position census `(12)`.  It does not rerun the heavy labelled
screen and therefore makes no physical supplier-label claim.

This theorem proves no new row closure or ledger decrement.  It proves no
physical `C2`, `C3`, `V4`, `S3`, or `PSL_2(Z)` action, no `z1=218` residual
support, and no THM-2056 flank.  It identifies the exact algebraic survivor,
the first failed implications, and the minimal missing carrier sidecars.
