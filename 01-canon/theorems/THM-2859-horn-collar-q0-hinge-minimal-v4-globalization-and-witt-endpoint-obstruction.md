---
id: THM-2859
title: "Horn-collar q0 hinge, generator-valued Z8 endpoint germ, minimal V4 globalization, and Witt descent obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The twenty-cell
  THM-2847 E3 horn meets the THM-2825 collar in one physical q0
  live/dead/live hinge.  Twelve long labelled paths realize one directed,
  ancestry-preserving endpoint translation Z^8, whose value generates the
  vertical endpoint C13 and is compatible with a first-carry reference; no
  physical C13 action or endpoint/ancestry intertwiner is proved.  The
  complete 685-path forest realizes only four of the thirteen vertical
  translates, and its 98 rootless paths supply none.  The formal/rootwise
  collar algebra is the punctured-V4 corner pM4p=M3, whose specified
  enveloping V4-action adds one rank-587 object.
  Full endpoint masks are rectangles; the q0/q3/q7/q11 masses are
  81/90/0/81, and the endpoint germ does not descend across this
  allocation/E3 boundary.  The literal q3/q11 attachment of THM-2863's
  four probes is empty before Prony splitting.  The joint
  residue/first-carry quotient is nonsplit C169, not an affine action on
  the split endpoint plane.  No row is excluded and LRC(14) remains OPEN.
source: root/lrc-horn-collar-witt-hinge-2026-07-28
depends_on:
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2818-right-cofiber-positive-copy-stratification-and-alternating-half-step-chains
  - THM-2825-nearest-half-step-common-right-collar-and-semantic-parity-boundary
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
related:
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2852-prime-power-orbit-spectrum-harvest-and-cayley-tournament-nonsingularity
  - THM-2861-endpoint-hermitian-edge-holonomy-and-semilinear-clutch-test
  - THM-2863-endpoint-prony-splitter-and-carry-character-three-intertwiner
  - THM-2870-prime-power-convolution-versus-physical-diagonal-intertwiner-obstruction
scripts:
  - 04-computation/lrc14_horn_collar_hinge_thm2859.py
  - 04-computation/lrc14_horn_collar_endpoint_carry_thm2859.py
  - 04-computation/lrc14_horn_collar_endpoint_orbit_action_thm2859.py
  - 04-computation/lrc14_horn_collar_rootless_endpoint_boundary_thm2859.py
  - 04-computation/lrc14_horn_collar_prony_typed_descent_gate_thm2859.py
  - 04-computation/lrc14_horn_collar_v4_globalization_thm2859.py
  - 04-computation/lrc14_horn_collar_witt_hinge_thm2859.py
  - 04-computation/lrc14_horn_collar_endpoint_coboundary_thm2859.py
outputs:
  - 05-knowledge/results/lrc14_horn_collar_hinge_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_endpoint_carry_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_endpoint_orbit_action_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_rootless_endpoint_boundary_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_prony_typed_descent_gate_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_v4_globalization_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_witt_hinge_thm2859.out
  - 05-knowledge/results/lrc14_horn_collar_endpoint_coboundary_thm2859.out
---

# THM-2859 -- horn-collar hinge, endpoint germ, and Witt descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem identifies a positive physical endpoint germ and the exact
place where it stops.  The THM-2847 horn and THM-2825 collar share a q0
hinge; farther along twelve of its labelled paths, one nonlocal directed
move preserves every audited physical decoration and translates both
endpoint masks by `Z^8`.  Since `8` is a unit modulo thirteen, this arrow has
a generator value in the vertical endpoint `C_13`.  The audit does not
supply a physical identity, inverses, composition, all thirteen powers, or
an intertwiner identifying endpoint `Z` with THM-2851's ancestry carry `T`.
What is absent is precisely that action/intertwiner and its descent to the
q3/q11/q7 allocation horn across the `E3`/complement boundary.

The same hinge also makes two previously formal obstructions concrete.  The
three collar corners form a partial `V_4` translation groupoid whose
specified enveloping action needs one new rank-587 object.  The
endpoint-address plane has the right 169-element cardinality for the joint
residue/first-carry quotient, but its inherited affine law is split and has
exponent thirteen.  The required residue lift is the nonsplit `C_169` law.
These are independent invoices: one generator-valued arrow over an existing
object does not create either a group action or a missing base object.

This is not a row discharge.  It changes the next construction from “find a
candidate carry value” to “extend the endpoint germ to an action, construct
the endpoint/ancestry intertwiner, and descend it while supplying a new
typed object.”  The LRC(14) ledger remains `165`.

## 1. The exact twenty-cell q0 hinge

Let `h=401080680` be the THM-2825 half-step and let

```text
I  =(142004992589460,142005019034340),
R  =I-2h=(142004190428100,142004216872980),
M1 =I-h =(142004591508780,142004617953660),
M2 =I.                                                    (1)
```

The THM-2847 `E3`-only horn consists of the twenty labelled cells

```text
clock=1,
sigma in {0,3,8,9,12},
target in {5,6,9,10}.                                    (2)
```

Reconstructing the complete THM-2825 common/right bank in `(2)` gives

```text
52 labelled right roots,
4,076 labelled common pieces.                            (3)
```

In every cell, root index zero has the literal path

```text
R --(+h)--> M1 --(+h)--> M2=I,                          (4)
```

with semantic pattern `(live,dead,live)`.  The twenty rooted path lengths
have census

```text
144:8, 118:4, 40:4, 14:4.                               (5)
```

All twenty labels forget to the one physical interval triple `(1)`.
Consequently this is one hinge with twenty labels, not twenty independent
physical fibres.

In factor order `(E3,clock,q1,q2,c2,c3)`, the only failed factor name at the
root is `E3`: the four source-native/source-pulled/target-native/target-pulled
frames are `011111/111111/111111/011111`.  The corresponding four frames at
`M1,M2` are all `111111`.  Their carrier and endpoint types are

```text
                    source carrier   target carrier   source endpoint
R                         empty          delta0               0
M1,M2                    delta0          delta0              81.          (6)
```

The target endpoint mask is the same 81-point mask at all three corners.
The source change in `(6)` cannot be repaired by any permutation of the
169 addresses.

For the pulled allocation intervals `I_q=I+q*T_period/13`, where `T_period`
is the physical period, complete interval membership gives

```text
I_0=M2 is in the collar;
I_q is outside common and right support for q=1,...,12.  (7)
```

Thus the q0 hinge is literal, while the q3/q7/q11 horn arrows are not
composable with it merely because their formal degrees add.

## 2. A physical generator-valued endpoint germ

On the twelve cells

```text
clock=1,
sigma in {0,3,12},
target in {5,6,9,10},                                    (8)
```

the distinguished path has at least 118 vertices.  Its index-two vertex is
`M2`; write `K` for index 68.  Then

```text
K=M2+66h.                                                 (9)
```

Direct evaluation of all 169 source and target endpoint addresses gives the
same unique translation on both sides.  Write `Z=T_(0,1)` for the vertical
endpoint-address translation:

```text
E_source(K)=Z^8 E_source(M2),
E_target(K)=Z^8 E_target(M2).                            (10)
```

The step difference `66` is even, so semantic nonvanishing is preserved.
Both vertices are common pieces, retain all six native and pulled factors,
and have source and target carrier `delta0`.  Their literal contributor sets
agree exactly:

```text
|U|=966606,   |V|=28534,
SHA256=15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd,
U symmetric difference=V symmetric difference=0.         (11)
```

The distinguished sheet is supplied at both ends.  Equations `(9)--(11)`
therefore define a fully decorated physical endpoint-address arrow, not
only a Hamming coincidence.  Again, twelve labels collapse to one physical
interval edge.

The positive is a translation plateau with a sharp edge, not an isolated
coincidence.
Relative to `M2`, the source endpoint remains the same `Z^8` translate
at every path index `68,...,90`; the target endpoint does so at
`68,...,89` and fails at index 90 in all twelve labels.  All 276 labelled
vertices in this range retain all factors, carrier, and ancestry chamber.
Semantic truth alternates with path parity, so the first two-frame live
reference is index 68 and the last is index 88.  Index 90 is a
source-frame-only live reference and the exact target-frame boundary.  These
vertices repeat the same `Z^8` value; they are not iterates of a proved
action.

The q0 endpoint mask is the rectangle

```text
H=1_A tensor 1_B,
A={0,1,2,3,4,5,6,7,12},
B={0,1,2,3,4,5,8,9,10}.                                 (12)
```

Its augmentation in `F_13` is `81=3`, so THM-2839 makes `H` a unit of
`F_13[C_13^2]`.  The signed current associated to `(10)` is

```text
D8=(Z^8-1)H.                                             (13)
```

It has eighteen `+1`, eighteen `-1`, and 133 zero coefficients.  Since

```text
Z^8-1=(Z-1)(1+Z+...+Z^7)
```

and the second factor has unit augmentation `8`, convolution by `D8` has
rank `13*12=156` over `F_13`; its vertical projection has rank `12`.
Equivalently, the vertical difference has the one-dimensional constant
kernel and image the twelve-dimensional augmentation ideal.  This is a
local group-algebra statement, not a semisimple character decomposition.
Independently, exact row reduction over `F_1000003` gives ranks
`169,156,156,12` for `H`, the `Z` and `Z^8` currents, and the vertical
projection, respectively.
THM-2870 moreover proves that no invertible same-mask intertwiner can
identify this unit convolution with the proper Boolean physical diagonal.
Thus the rank statement is not an endpoint action.

Because `gcd(8,13)=1`, `Z^8` is a generator of the abstract vertical
endpoint `C_13`.  Thus a demand for the literally normalized label `Z`
would be coordinate-dependent.  What `(10)` proves is one
generator-valued endpoint-side reference compatible with the first-carry
coordinate.  It does not identify `Z` with the ancestry carry `T`, prove a
physical `C_13` action, give a `C_169` lift, or produce an allocation-q
arrow.

## 3. Full endpoint rectangles locate the descent failure

An independent evaluator applies the defining periodic inequalities directly
to all `55,341` vertices in the 587 rooted collar paths.  After the extra
hinge/allocation controls it caches 1,199 physical interval masks, with 73
distinct source masks on the paths.  Every full mask factors as a Cartesian
rectangle

```text
E_x=A_x times B_x in F_13^2.                             (14)
```

This factorization follows from the endpoint representative

```text
REPS(a,b)=(0,-a,-b,0,0,0,a,b,0):
```

the address inequalities separate into an `a` family and a `b` family.
The complete vertex census by `(|A|,|B|,|E|)` is

```text
(0,0,0):          81
(9,9,81):      17,286
(9,10,90):      4,798
(9,11,99):      2,306
(10,9,90):     20,335
(10,10,100):    7,798
(10,11,110):    2,737.                                  (15)
```

On all 54,754 adjacent `+h` edges, and independently on all 54,167 `+2h`
edges, a nonconstant mask difference has augmentation in

```text
{+/-9,+/-10,+81};                                       (16)
```

augmentation zero occurs exactly when the two masks are equal.  Hence every
nonconstant adjacent edge fails both:

- permutation or translation conjugacy of masks; and
- representation as `(T_d-I)phi` over `Z` or `F_13`.

The path graph is a forest, so this failure precedes holonomy: every
nonconstant adjacent edge in the stated censuses lacks a lawful
group-valued label that could be integrated.

The complete search of all `4,208,957` forward pairs `i<j`, summed over all
587 rooted paths and using full source endpoint masks, finds nonzero
translations only in the directions

```text
(0,1), (0,2), (0,3), (0,7), (0,8), (0,9),               (17)
```

with counts `50,712`, `14,812`, `11,109`, `119,560`,
`42,547`, and `9,730`.  It recovers `(10)` but finds no `(0,4)` edge and no
composable physical `8+9=4` triangle, even after swapping coordinates.
Endpoint-address differences in `(17)` are not allocation-q arrows.

A second complete audit includes all 685 components of the THM-2825 forest:
587 cofiber-rooted paths and the 98 common-only paths carrying its signed
rootless `C_13` units.  Across all 63,895 labelled vertices and 1,275 unique
physical intervals, the vertical orbit `{Z^kH:0<=k<13}` meets the source and
target masks only at

```text
k in {0,4,8,9};  unique intervals by k: 12,23,23,23.   (17a)
```

The common-only paths contribute no such hit.  Since `H` has trivial
vertical stabilizer, no endpoint-mask-preserving `C_13` object action
containing `H` can be realized wholly inside the complete inherited forest.
This zero is structural.  The 98 paths de-duplicate to 286 physical
intervals.  Wrong mass rejects 218 source and 220 target masks; among the
mass-81 survivors, a wrong horizontal factor rejects another 11 and 10.
The remaining 57 and 56 masks have one of the three vertical factors

```text
{0,1,4,5,8,9,10,11,12},
{0,1,4,7,8,9,10,11,12},
{0,1,6,7,8,9,10,11,12},
```

none a cyclic translate of the q0 factor `B` in `(12)`.

The directed path relation has only plateaux and the two `Z^8`
label-transition types `0->8` and `9->4`; it has no forward `Z^5` inverse
and no composable `Z^8` square.  There is a useful support reframe:

```text
{0,4,8,9}={0,8}+{0,9},             8+9=4 mod 13.        (17b)
```

Thus all four vertices of an affine `8/9` parallelogram in the ambient
`C_13` Cayley graph occur, but the restricted orbit relation retains only
the parallel `+8` sides `0->8` and `9->4`, not either `+9` bridge
`0->9` or `8->4`.  This is support-dual to the punctured-`V_4` invoice:
here all four mask vertices occur but bridge arrows are absent; there the
fourth base object is absent.  It is an analogy, not a common grading,
`V_4` action, or endpoint/ancestry identification.

On de-duplicated physical intervals the literal `+66h`
displacement realizes the expected `Z^8` value on only 22 of 81 intervals
in each endpoint frame, while `+132h` sends all 81 intervals inside the
forest but outside this thirteen-mask orbit.  Thus THM-2825's rootless
signed units do not provide the missing endpoint action.  This is an
internal-forest obstruction, not a no-go for new objects or external
correspondences.

Finally, direct endpoint evaluation on the four horn allocations gives the
same source and target mask at each q and the masses

```text
q0:81,   q3:90,   q7:0,   q11:81.                       (18)
```

The q0 and q11 rectangles have equal mass but lie in distinct split
translation orbits.  The q7 mask is empty.  Thus an allocation clutch must
transport rectangle factors, endpoint presence, and `E3` truth; q degree
alone discards exactly the obstructing data.

## 4. The formal/rootwise collar algebra is a punctured `V_4` groupoid

Let

```text
G=V_4=F_2^2,
U={00,11,01}=G\{10},
R=00, M1=11, M2=01.                                     (19)
```

Restrict the translation groupoid of the free `G`-torsor to `U`.  It is the
three-object pair groupoid.  If `E_vu` denotes its matrix unit and is graded
by `v-u`, then

```text
A=End(K[U])=M3(K),
dim A_00=3,
dim A_g=2 for every nonzero g in V_4.                    (20)
```

For `g=h!=0`, `A_gA_h` misses one diagonal matrix unit.  For distinct
nonzero `g,h`, `A_gA_h` contains one of the two orientations in
`A_(g+h)` and misses the other.  Therefore each of the nine ordered products
of nonzero degrees has a one-matrix-unit support defect.  Tensoring with the
587-dimensional root space makes every defect rank 587.

A nonzero scalar cocycle, sign twist, or Schur multiplier rescales supported
products and cannot create a missing matrix unit.  In fact the three-object
pair groupoid is contractible, so every scalar two-cocycle on it is a
coboundary.  The obstruction is support, not phase.

There is no transitive three-point `V_4`-set: its orbit sizes are `1,2,4`.
Keeping the three existing objects distinct therefore forces the free
four-point torsor.  Up to equivariant isomorphism fixing `U`, this is the
unique minimal enveloping `V_4`-set/action.  If
`p=e_00+e_11+e_01`, its algebraic form is

```text
B=M4(K)=End(K[V_4]),
A=pBp=pM4(K)p=M3(K).                                     (21)
```

It adds the missing object `X` of root rank 587 and seven block matrix units,
each of rank 587.  Equivalently, in the labelled-root groupoid this is
`7*587=4,109` rank-one decorated arrows; it is not 4,109 independent block
algebra generators.

The same invoice has an exact mapping-cone form.  Write the collar
isomorphisms as

```text
d:R->M1,   a:M1->M2,   S=ad:R->M2.
```

With `X` a copy of `R`,

```text
D0(v)=(dv,v),
D1(m,x)=a(m)-S(x)                                       (22)
```

gives the contractible exact sequence

```text
0 -> R --D0--> M1 direct-sum X --D1--> M2 -> 0.          (23)
```

Conversely, any factorization through a filler `X` has
`dim X>=rank S=587`, proving minimality.

This completion is not physical in the inherited category.  All 587 roots
have empty source carrier and at least one native-factor hole; all `M1,M2`
objects have source carrier `delta0` and all six factors.  The root semantic
population is 573 `(live,dead,live)` triples versus 14
`(dead,live,dead)` triples, so no semantic-flipping permutation of the
existing bank supplies the missing copy.  The formal `X` moves the
empty-to-present boundary to a new edge; it does not remove it.

Graphically, `(4)` is the three-vertex path obtained by deleting one corner
from the `V_4` Cayley square.  Its composite `S` has partial-cube distance
two.  Declaring `S` a new unit edge makes a triangle, which is not a partial
cube.  Within this fixed two-direction `V_4` square, the unique partial-cube
repair is the fourth vertex, not a “teleporting” diagonal.

## 5. Split endpoints versus the nonsplit residue/first-carry quotient

Write a joint residue/first-carry state as the base-thirteen digit pair

```text
n=13a+q.
```

The natural residue lift

```text
L_h(a,q)=(a+floor((q+h)/13), q+h mod13)                  (24)
```

satisfies

```text
L_k L_h=T^floor((h+k)/13)L_(h+k mod13),
L_1^13=T,
ord(L_1)=169.                                            (25)
```

Consequently the sharp abstract state object is the nonsplit extension

```text
0 -> C_13 -> C_169 -> C_13 -> 0,
C_169=(W_2(F_13),+).                                     (26)
```

The inherited endpoint plane `F_13^2` has the same 169 points but the wrong
operation.  For every odd prime `p`, every `p`-element of
`AGL(2,F_p)` has exponent at most `p`: in homogeneous coordinates it is a
unipotent `3 x 3` matrix `I+N`, with `N^3=0`, and
`(I+N)^p=I`.  Hence `AGL(2,F_13)` has no element of order 169.  Exact
enumeration gives 28,561 affine 13-primary candidates: the identity once
and order thirteen in the other 28,560 cases.

The endpoint value in `(10)` is consistent with, but much smaller than,
`(26)`.  If a future action identifies endpoint `Z` with the first-carry
kernel, the presentation

```text
<L,T | T^13=1, [L,T]=1, L^13=T^8>                       (27)
```

is still cyclic of order 169 because `8` is a unit modulo thirteen.  The
current arrow merely has the generator value `Z^8`; it does not supply the
kernel action, the identification `Z=T`, the lifted order-169 element `L`,
or its q/E3 descent.

This also exposes an exact phase interface with THM-2857.  That theorem has
the free Galois orbit

```text
c_r=A-B omega^(3r),   sigma_1(c_r)=c_(r+1).
```

For the physical displacement `66h`, the exact source- and target-phase
exponents reduce modulo thirteen to

```text
source x-sweep:       7,
target endpoint sum:  6.                                (28)
```

Both identities are verified in the two exact cyclotomic specializations;
the source and target characters are inverse, and the source/target overlap
is unchanged.  If `b` is the vertical endpoint coordinate in `(10)`, the
phase-matching rechart is `r=10b`.  The physical move `b -> b+8` gives
`r -> r+2`, and hence

```text
c_(r+2)-A=omega^6(c_r-A).                               (29)
```

Thus, after this chosen rechart, the target phase equals the faithful
THM-2857 Galois character under `sigma_2`, while the source carries its
inverse and the paired phase cancels.  This is an exact character
comparison, not an identification of the physical coordinate `b` with the
Galois torsor coordinate `r`.  Equations `(28)--(29)` are not the missing
semilinear clutch: the proved ancestry action is coefficient-linear and
fixes `c`.  A future dual edge operator must transport the physical fibre by
`b -> b+8`, act on the target coefficient line by `sigma_2`, and act on the
source line by `sigma_(-2)` so that the paired phase cancellation is
respected.

THM-2863 subsequently supplies a unique normalized `C_13`-equivariant
isomorphism between the one-dimensional endpoint character-three line and
the corresponding projection of the THM-2851 carry derivative.  That is a
coefficient-line interface, not the physical endpoint/ancestry action,
simultaneous full current, signed projection, or `E3` transport required
here.

### 5a. The literal `Z^8`-to-Prony attachment fails before Prony

The cheapest proposed composition has now been tested on the same twelve
labelled step-2-to-68 edges.  For each edge and `q in {3,11}`, translating
both endpoints by the centered q displacement

```text
q3:  +171366h,             q11: -114244h
```

produces 48 state rows and 24 edge rows.  None of the 48 translated states
is a vertex of the same labelled cell's complete collar forest, so the
strict same-sheet bank has `0/24` edges; none is a forest vertex even after
forgetting the cell label.  This is the literal obstruction already
foreshadowed by `(7)`, now checked on both ends of `(10)`.

A hostile relaxation keeps the physical base-root edge and merely attaches
the q-shifted current.  All `48/48` states and `24/24` edges retain all six
factor signatures, including `E3`, and both weighted carriers.  Their
q-shifted ancestry is constant across `(10)`, but its `U` component is
empty, so this is not a positive full-ancestry fibre.  The full endpoint
rectangles are decisive: at both endpoint frames and for both q values the
exhaustive translation set is empty.  Imposing the expected `Z^8` gives
Hamming defects

```text
q3: 80,                    q11: 72.                     (29a)
```

These numbers are not claimed as minimum Hamming distances.  Structurally,
the step performs the same local vertical surgery

```text
B -> (B\{3,4}) union {6,7}.                             (29b)
```

At q0 this happens to be the `Z^8` translate of the special factor `B` in
`(12)`.  At q3 and q11 the before/after vertical gap necklaces differ, so
`(29b)` is not any cyclic translation.  For q3 the two length-two
complement blocks merge into one length-four block; for q11 their cyclic
separation changes.

The named THM-2863 origins have vertical coordinate zero, outside this
surgery, and their occupancies remain `(1,0)` at q3 and `(1,1)` at q11.
Consequently each of the four multiplier rows
`Y=38,64,90,116` has an empty typed edge bank; no endpoint coefficient is
declared zero, and the Prony split and character-three comparison are
skipped.  This excludes only the literal same-labelled-cell or relaxed
base-root/full-positive-mask attachment.  It does not obstruct THM-2863's
proved coefficient-line interface, THM-2868's reserved signed/projective
route, or a new external correspondence.

## 6. Why the two invoices do not cancel

The `Z^8` edge supplies one generator-labelled morphism over the common
`E3` object, not kernel isotropy.  The missing `V_4` degree asks for a new
base object.  If the germ first extended to a physical `C_13` isotropy
action, abstract conjugation would spread that action over the three-object
pair groupoid and produce

```text
U times C_13 times U,
M3(K[C_13]).                                              (30)
```

The inherited data do not permit this conjugation across `S:R->M2`: the
source endpoint fibres have cardinalities zero and 81, so no bijective
intertwiner exists on which `S^(-1)Z^8S` could act.  Central carry
extensions decorate
composable arrows; they cannot create an absent intermediate object.

The formal grade coincidence remains a useful locator.  The horn changes
`(outer word,E3 truth)` by `(1,1)` and the even collar changes
`(semantic,carrier presence)` by `(0,1)`.  Only after a proposed sidecar maps
both differently typed grading pairs into one common `V_4` does their sum
select the missing degree

```text
(1,1)+(0,1)=(1,0).                                      (31)
```

This is a mnemonic locator, not a common grading.  Equations `(7)`, `(18)`,
and `(20)--(23)` show why equality of grade names is not composition: it
forgets physical interval, endpoint rectangle, macro truth, and object
support.

## 7. Exact remaining obligation

The result changes the live task to a typed descent problem:

1. anchor the proved ancestry-preserving endpoint `Z^8` reference inside the
   q0/common-`E3` block;
2. construct physical identity, inverse, composition, and all endpoint
   powers, or otherwise prove an endpoint-to-ancestry carry intertwiner;
3. replace the failed literal attachments in Section 5a by a genuinely new
   common-gauge full-current correspondence, then combine THM-2863's
   normalized character-three interface and THM-2861's Hermitian edge with
   a lawful signed projection and semilinear lift of `(28)--(29)`;
4. extend to the nonsplit joint residue/first-carry quotient `(24)--(27)`;
5. transport the resulting object through q3 and q11 to q7 while preserving
   the rectangle factors in `(18)` and the `E3`/complement semantic word;
6. supply the rank-587 fourth cofiber object required by `(21)--(23)`.

The cheapest hostile cell is `(clock,sigma,target)=(1,0,5)`.  It contains
both the positive long-path reference and the minimal empty-to-81-point
hinge, so a proposed descent can be rejected before any all-row computation.
That downstream four-probe test is now closed negatively by Section 5a.
The four-label orbit fragment `(17a)`, absence of `(0,4)`, lack of a
physical `8+9=4` triangle, and empty q3/q11 typed attachment are the current
stopping reasons for the inherited collar.

## 8. Verification boundary

Eight companions use explicit runtime failure gates and no executable Python
`assert` statements.

- The hinge companion checks every q0 through q12 interval membership and
  independently reproduces the twenty-cell physical/factor/carrier hinge.
- The endpoint/carry companion reconstructs all twenty hinge cells, both
  endpoint sides, the q-allocation masks, the exact convolution ranks, all
  physical decorations of `(10)`, and the literal ancestry sets `(11)`.
- The orbit/action companion evaluates the complete 685-path forest,
  separates labelled from physical occurrences, and checks every vertical
  translate of `H`, the directed action laws, and fixed `+66h` powers.
- The rootless-boundary companion explains the 98 common-only paths' zero
  orbit hits by mass, horizontal factor, and three exact non-orbit vertical
  subsets.
- The typed-descent companion tests all 24 q3/q11 attachments, the relaxed
  factor/carrier/ancestry bank, full rectangle translation orbits, and the
  four THM-2863 probe gates.
- The globalization companion independently rebuilds all 587 root triples
  from the THM-2818 physical constructor and checks the graded product
  defects, semantic/factor/carrier types, sign cohomology, partial-cube
  repair, and mapping-cone ranks.
- The dependency-free Witt companion exhausts the natural carry law and all
  affine 13-primary candidates.
- The endpoint companion independently evaluates all full masks from the
  periodic inequalities, reproduces the canonical root/M1/M2 controls, and
  searches every within-path pair.

Normal, optimized, and stored outputs are byte-identical after LF
normalization.  The exact script and output hashes are recorded in the
results index.  The finite universes above are complete only for the named
THM-2825/2847 bank; no affine-only no-go is asserted for arbitrary
permutations of 169 abstract states.
