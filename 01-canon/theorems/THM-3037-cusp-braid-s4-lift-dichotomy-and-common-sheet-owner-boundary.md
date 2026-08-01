---
id: THM-3037
title: "Cusp-braid S4 lift dichotomy and common-sheet owner boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Two distinct transposition meridians of the matching S3 quotient have
  exactly eight braided lifts to S4.  Four are transposition pairs generating
  a point-stabilizer S3 with one common fixed quartic sheet and quartic
  parabolic order two; four are 4-cycle pairs generating S4 with no common
  sheet and parabolic order four.  Their squares are respectively zero or two distinct nonzero
  V4 elements generating the kernel.  Keller inertia or a simple/index-free raw
  discriminant order eliminates the width-four branch, but the surviving
  section belongs first to the finite normalization: an affine-source owner
  and a point over the cusp require separate sidecars.  For a connected
  global degree-four cover, the local point-stabilizer instead forces global
  monodromy S4.  No quartic Keller, G1, JC(2), DC(2), or LRC exclusion follows.
source: codex-quartic-cusp-braid-lift-2026-08-01
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary
  - THM-3034-ordered-quartic-cross-wall-x1-14-and-diamond-quotient
  - THM-3035-level-two-farey-anharmonic-quartic-s3-orbit-diamond
script: 04-computation/quartic_cusp_braid_s4_lift_dichotomy_thm3037.py
output: 05-knowledge/results/quartic_cusp_braid_s4_lift_dichotomy_thm3037.out
script_sha256: ec28cfcf1056562f4bc537f8f60bea0d3ff3041a04e00e563c1428661a54fa69
output_sha256: c2c4a237be1dda2c310ed281436c16e1090fae1504feb634f2ec0ddae79cdc5c
hash_basis: LF-normalized bytes
---

# THM-3037 -- the quartic braid lift has one binary square defect

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and exact statement

[THM-2598, the quartic `V4` resolvent torsor](THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary.md)
identifies the canonical quotient

```text
pi:S4 -> S4/V4 = S3                                      (1)
```

as the action on the three pair matchings.  [THM-2646, the braid-three
central pullback](THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md)
shows why the order-two and order-three modular factors meet on a braid and
why an affine `V4` lift with surjective `S3` shadow kills its centre.  The new
question here is finer: after the centre dies, which quartic lift of the two
geometric braid meridians remains?

Let `tau_1,tau_2` be distinct transpositions in `S3`.  They obey

```text
tau_1 tau_2 tau_1 = tau_2 tau_1 tau_2.                    (2)
```

Among the sixteen pairs

```text
(X,Y) in pi^-1(tau_1) x pi^-1(tau_2),                    (3)
```

exactly eight obey the same braid relation.  They form two free simultaneous
`V4`-conjugacy orbits of size four:

| branch | cycle types of `X,Y` | `<X,Y>` | common fixed sheets | quartic parabolic order |
|---|---:|---:|---:|---:|
| split | transposition, transposition | `S3` | one | `2` |
| full | four-cycle, four-cycle | `S4` | zero | `4` |

There is no mixed braided lift.  The branch bit is intrinsic:

```text
split:  X^2=Y^2=1;
full:   X^2,Y^2 are distinct nonzero elements generating V4. (4)
```

Thus the user's `V4` tree move is literal rather than analogical.  The two
quotient involutions either lift without translation defect, or their two
squares recover the full Klein kernel.

Put

```text
A=XYX,                  B=XY.                              (5)
```

In every braided lift,

```text
(XY)^3=1,       A^2=B^3=1,       AB=X^-1.                 (6)
```

Consequently the representation factors through

```text
PSL2(Z)=C2*C3=<A,B | A^2=B^3=1>.                          (7)
```

The geometric cusp meridian is `X=B^-1 A`, equivalently the modular
parabolic is `AB=X^-1`.  Its order is two in the split branch and four in
the full branch.  This is the maximal cusp width in this four-sheet
permutation action, not the full cusp-width multiset of THM-2975's six-sheet
actions.  It is the same square defect viewed as a Farey cusp-width jump.
It is essential not to apply a sheet-fixing test separately to `A`
or `B`: the inertia element is `X`, and the full branch can have a
sheet-fixing modular involution `A` while `X` is a fixed-point-free
four-cycle.

## 2. Semidirect proof of the eight lifts

Write

```text
S4 = V : H,          V=F2^2,          H=GL2(F2)=S3.       (8)
```

All ordered distinct transposition pairs in `H` are conjugate, so take

```text
s(u,v)=(u+v,v),              t(u,v)=(v,u).                (9)
```

Every lift has a unique form

```text
X=(a,s),       Y=(b,t),       a=(a1,a2), b=(b1,b2).       (10)
```

Expanding the affine multiplication law turns the braid relation into the
single equation

```text
b1=a2+b2.                                                   (11)
```

Hence `a1,a2,b2` are free and there are exactly eight lifts.  Their squares
are

```text
X^2=(a2,0),                    Y^2=(a2,a2).                (12)
```

For `epsilon=a2=0`, both lifts are involutions.  They generate a
point-stabilizer `S3` and fix a unique common quartic sheet.  For `epsilon=1`,
the squares are two distinct nonzero vectors spanning `V`; the quotient of
`<X,Y>` is all `H`, so `<X,Y>=V:H=S4`.  Both lifts then have order four and
no fixed sheet.  Simultaneous conjugation by `V` is free and transitive on
the four solutions of either fixed `epsilon`, proving the orbit statement.

There is a cohomological form of the count.  Modulo simultaneous translation,
affine lifts of the standard modular shadow are classified by

```text
H^1(C2*C3,V).                                               (12a)
```

For the natural two-dimensional `F2` module, `V^(C3)=0` and Maschke
semisimplicity gives `H^1(C3,V)=0`.  The order-two generator is a
transvection: `V^(C2)=im(1-A)` is one-dimensional and
`H^1(C2,V)=0`.  The free-product Mayer--Vietoris sequence therefore gives

```text
H^1(C2*C3,V) = V / V^(C2) = F2.                            (12b)
```

Its zero and nonzero classes are exactly the split and full square-defect
orbits.  Thus neither free factor alone carries the new bit; it is the
gluing coordinate left by their co-occurrence.

This also proves (6): the usual braid centre is
`(XY)^3=(XYX)^2`; it vanishes in both finite branches, and direct braid
algebra gives `X=B^-1 A` and `AB=X^-1`.  The computation companion checks
the result for all six ordered distinct transposition pairs, not just the
chosen affine gauge.

## 3. Conditional `A2` cusp application

The following geometric statement is deliberately conditional on an actual
cusp chart.  Let `B_0` be a small normal complex two-ball, let
`Delta subset B_0` be an irreducible `A2` cusp, and put
`U=B_0\Delta`.  Suppose a normalized degree-four finite cover over `U` has
matching quotient equal to the standard cubic root-cover monodromy.  This
hypothesis is supplied, for example, by the local depressed cubic
`T^3-3uT-2v` after appropriate units and etale coordinates; topological
monodromy alone is not claimed equivalent to that analytic/order normal
form.  Choose a
Wirtinger pair of **genuine geometric meridians** of `Delta`; downstairs
they are distinct transpositions, and upstairs call them `X,Y`.  Then the
cover lies in exactly one of the two branches of Section 1.

Suppose in addition that this is the local normalization of a polynomial
Keller graph and that the selected meridians are the actual inertia of its
affine Jelonek divisor.  [THM-2633, affine inertia normal generation and
fixed branches](THM-2633-derangement-character-obstruction-and-d4-keller-exclusion.md)
forces a geometric divisor meridian to fix a finite affine sheet.  The full
branch has inertia `X` a four-cycle, so it is impossible.  Therefore the
local braid lift is the split branch:

```text
<X,Y>=S3,        Fix(X) intersect Fix(Y)={one sheet},
quartic parabolic order / maximal cusp width=2.             (13)
```

This is a local statement.  It does not put the **common** fixed sheet into
the original affine source.  Under global connectedness it actually forces,
rather than excludes, `S4` monodromy as shown in Section 5.

## 4. Independent simple-order discriminant selector

There is a second selector which does not use THM-2633.  Work in
characteristic zero at the generic point of `Delta`, and write

```text
e_j = raw discriminant valuation of the degree-j order,
d_j = tame discriminant valuation of its normalization,
i_j = index valuation,              e_j=d_j+2i_j.         (14)
```

Assume here that the local quartic order and cubic order are the actual
standard monic-integral resolvent pair, so their **order discriminants**
agree up to a unit.  Equality only of zero sets is not enough.  THM-2598
supplies this exact equality for its paired orders; local monicization by a
unit does not change the valuation.

The matching quotient meridian is a transposition, so `d_3=1`.  The split
quartic meridian has `d_4=1`, whereas the full four-cycle has `d_4=3`.
The raw-order identity gives `e_4=e_3`.  Hence

```text
split:  i4=i3;
full:   i3=i4+1.                                           (15)
```

If the cubic cusp order is index-free/maximal at the generic divisor, its raw
discriminant has simple divisor order `e_3=1`.  Then `e_4=1`, which is
incompatible with `d_4=3`; only the split branch remains.  Without this
order hypothesis the full branch is not excluded: it is paid for by the
one-unit index tax in (15).  This is exactly the local-order boundary in
THM-2598 and must not be erased by comparing only discriminant polynomials.

## 5. What the common sheet really supplies

In the split branch, the local monodromy subgroup has orbits `1+3` on the
four sheets.  Therefore the finite etale cover over `U` decomposes into a
degree-one component and a degree-three component.  The closure of the
degree-one component in the finite normalization over the normal ball
`B_0` is finite and birational over `B_0`, hence is `B_0`.  Thus the finite
normalization has a genuine local section.

This conclusion is stronger than a permutation fixed point and weaker than
an affine owner.  Zariski main factorization permits the original affine
graph source to be an open subset of the finite normalization.  It must
retain all four generic sheets over `U`, but it can omit the boundary trace
of the degree-one sheet over `Delta`; even when that trace is retained
generically, it can still omit the cusp point itself.

The smallest local owner hostile is

```text
(T-1)(T^3-3uT-2v),              Delta: v^2=u^3.            (16)
```

On the normalization `u=r^2,v=r^3`, the cubic factor becomes

```text
(T-2r)(T+r)^2.                                          (17)
```

The sheet `T=1` is fixed by the whole local `S3`, while `T=2r` is the
simple branch that can pay the existence of one affine branch over the
divisor.  In the finite normalization, take the abstract Zariski-main open
which keeps `T=1` off `Delta`, deletes only its boundary trace over `Delta`,
and deletes the cubic ramification locus `u=T^2`.  It retains the simple
branch `T=2r` for `r!=0`.  All four generic sheets over `U` remain, the
chosen local open is etale, and the only retained generic divisor branch is
not the common sheet.  This model is not asserted to be a polynomial Keller map;
it isolates the first invalid implication.  If an actual source has
`k_D=2`, the common sheet is generically present as the second branch, but a
separate properness or closure argument is still needed to land the cusp
point.

The hostile already passes the local conductor test.  If
`g(T)=T^3-3uT-2v`, then

```text
g(1)=1-3u-2v                                                (18)
```

is a unit at the cusp, so the Chinese-remainder idempotent
`g(T)/g(1)` belongs to the raw quartic order and selects the `T=1` factor.
Thus even an order-level fixed-sheet idempotent does not force its divisor
trace to lie in the affine open.

Equivalently, in completed local algebra the split normalization has a
canonical product decomposition into its singleton and cubic factors.  The
general owner test has two stages: whether the singleton idempotent lies in
the completed quartic graph order (or its conductor), and whether its centre
lies in the affine Zariski-main open.  Equation (18) shows that the second
stage remains load-bearing even when the first succeeds.  Both bits are
invisible to the discriminant.

[THM-2570](THM-2570-jelonek-cusp-cylinder-normalization-and-conductor.md)
performs an exact normalization/conductor/boundary analysis for its own
sporadic Jelonek cusp cylinder.  It supplies the right kind of sidecar, but
there is no asserted base-ring identification between that chart and the
quartic order here.

There is an independent permutation/group-level global hostile.  With split representatives

```text
X=(1 2),          Y=(2 3),
```

the local group fixes sheet `4`, but adjoining remote monodromy `(1 4)`
generates all of `S4`.  A local common sheet therefore does not reduce the
global monodromy group.

For a connected global degree-four cover this is forced, not merely a
hostile possibility.  Let `G<=S4` be the transitive global geometric
monodromy and let `H=<X,Y>` be the split local subgroup.  If `x` is its
unique common fixed sheet, then

```text
H = (S4)_x = S3,             H <= G_x <= (S4)_x.          (19)
```

Thus `G_x=H`.  Transitivity and orbit--stabilizer give

```text
|G|=4|G_x|=24,                 hence G=S4.                (20)
```

Changing the base path merely conjugates `H` and `x`.  Connectedness is
load-bearing: a disconnected global `1+3` cover may have only the local
`S3`.  For an irreducible degree-four graph/function-field cover, generic
connectedness supplies the required transitivity.

## 6. Relation to the existing modular and grade-three frontiers

The full branch is conjugate to the width-four `S4` modular quotient behind
[THM-2975, the six-sheet Schreier boundary](THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary.md).
It is now a sharp **non-Keller local hostile**: the free factors and killed
centre are valid, but the geometric parabolic is a derangement.  The split
branch instead supplies a local point-stabilizer `S3` and a section in the
finite normalization.

This also sharpens the original binary/ternary-tree proposal without
reversing the guardrails of
[THM-2596](THM-2596-modular-free-factor-farey-gram-owner-cocycle.md).
The rooted binary Farey children are parabolic **words** in the two torsion
factors, not the `C2` factor by itself.  The literal co-occurrence datum in
the present quartic action is instead the one-bit quotient (12b).  In its
nonzero class,

```text
X^2, Y^2, X^2Y^2
```

are exactly the three nonzero `V4` directions.  These are the three
parity/resolvent channels of
[THM-2606](THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin.md),
now recovered as a square defect of two lifted braid meridians.  This does
not identify the LRC parity and ternary towers with a common modular action;
[THM-2632](THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar.md)
still separates the permutahedral six-state action from the cyclic CRT
coordinate.  A physical LRC intertwiner remains absent.

Neither branch identifies the local cubic with the special affine cubic of
THM-1310.  [THM-2681](THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion.md)
proves that exact field identification impossible on its principal chart.
Consequently the depressed law, cuspidal law, Jelonek lead, or other
grade-three anatomy may transfer only after an explicit base-ring map,
valuation match, and common affine owner are supplied.  Equality of raw
discriminants alone is insufficient.

The exact information ledger is:

| operation | preserved | destroyed | sidecar still needed |
|---|---|---|---|
| `B3 -> C2*C3` | modular braid class after centre kill | integer full-twist height | discriminant winding |
| `S4 -> S3` | matching monodromy | `V4` square defect | lifted meridians |
| split local cover -> finite section | one common normalization sheet | affine ownership | Zariski-main open/closed incidence |
| split local -> connected global | point-stabilizer `S3`; transitivity | its common fixed sheet | group is forced `S4`; owner transport remains missing |
| resolvent cusp -> grade-three model | cusp/discriminant shadow | affine coordinates and valuations | explicit field/base-ring realization |

No quartic Keller-map exclusion, `G1`, Jacobian conjecture, Dixmier
conjecture, or LRC(14) decrement is asserted.

## 7. Exact companion

The companion constructs the matching quotient `S4->S3` from first
principles, checks the homomorphism on all `24^2` products, enumerates all
six ordered downstairs meridian pairs and all `96` candidate lift pairs,
and verifies the `48` braided lifts, two `V4` orbits per downstairs ordered
pair, square defects, modular
relations, tame normalized exponents `d=degree-number of inertia orbits`,
and the remote-monodromy hostile:

```text
python 04-computation/quartic_cusp_braid_s4_lift_dichotomy_thm3037.py
python -O 04-computation/quartic_cusp_braid_s4_lift_dichotomy_thm3037.py
```

Both modes must LF-normalize byte-for-byte to
`05-knowledge/results/quartic_cusp_braid_s4_lift_dichotomy_thm3037.out`.
