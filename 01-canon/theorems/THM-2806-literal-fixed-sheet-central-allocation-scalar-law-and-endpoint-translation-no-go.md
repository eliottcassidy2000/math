---
id: THM-2806
title: "Literal fixed-sheet central allocation scalar law and endpoint-translation no-go"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  A marked
  wandering sheet realizes the raw carrier-allocation
  arrays (B,P,Q,H)=w(1,delta_0,delta_0,delta_(0,0)).  Their only fourfold
  primal point is flat, while their central Fourier coefficients are
  w(13^2,13,13,1) with nonzero D3=144w.  The central face is the marginal
  of a 12-by-12 bare-only complement, not a nonzero pre-Fourier common-atom
  face.  Literal allocation is projectively scalar with zero endpoint step;
  the separately nondegenerate moving-sheet parallelogram has no affine
  covariance with it.  The marked representative gauge is exact, but the
  fixed rail-eight bank has common support on only clocks 1,2,3.  A
  separately scoped FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED addendum
  proves that THM-2807's specified tau-twelve address simplex is
  target-only: its common allocation carrier is empty before address
  restriction.  This addendum closes only item 1 of THM-2807's proposed
  test on that fixed simplex.
source: root/lrc-fixed-sheet-allocation-2026-07-28
audit: >
  d3-torsor-descent and jc-chart-entry-2026-07-28 (independent type/scope
  audits, raw-versus-central hostile, finite-abelian Fourier qualification,
  normal/-O/stored/hash replay: ACCEPT); absent-present-constructor-2026-07-28
  (fresh implementation of the literal allocation, marked gauge, affine
  covariance exhaustion, and stored-output replay: ACCEPT)
depends_on:
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
related:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2802-normalized-unit-bispectrum-and-projective-cyclic-orbit-reconstruction
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
script: 04-computation/lrc14_literal_fixed_sheet_allocation_thm2806.py
output: 05-knowledge/results/lrc14_literal_fixed_sheet_allocation_thm2806.out
script_sha256: 311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e
output_sha256: a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f
secondary_script: 04-computation/lrc14_literal_fixed_sheet_allocation_independent_audit_thm2806.py
secondary_output: 05-knowledge/results/lrc14_literal_fixed_sheet_allocation_independent_audit_thm2806.out
secondary_script_sha256: 90386206fa441edaa121b54b501347f14005e997b4c0e8f95947f2fac14050b4
secondary_output_sha256: 36e4d0f2e48bb8bf3aa2673f71e82e5d30e6c84830dc8ea7ff35c6d2b8e61569
addendum_status: >
  FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the fixed
  (clock,sigma,tau)=(1,0,12) rail-eight sheet, the common carrier is empty
  before selecting any of THM-2807's three addresses, while the right
  cofiber contains one whole positive cylinder at each address.  The
  addendum has its own primary and independent implementations and does not
  inherit the theorem's earlier audit badge.
addendum_script: 04-computation/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.py
addendum_output: 05-knowledge/results/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.out
addendum_script_sha256: 7fea18161046f8de35b2e6ef04c88a13485de61045bc363f89c5ebfec8f76480
addendum_output_sha256: cba545e7beff4fe889c76ae681c47c806969ccd33aa79def54527709a6ffafc6
addendum_audit_script: 04-computation/lrc14_tau12_simplex_allocation_support_no_go_independent_audit_thm2806.py
addendum_audit_output: 05-knowledge/results/lrc14_tau12_simplex_allocation_support_no_go_independent_audit_thm2806.out
addendum_audit_script_sha256: fd45e0959f27021dcf67a6f22a1be426db3849d4c0a806ba5135e17e5a0d57fe
addendum_audit_output_sha256: e18b61016df622283b8bfccb4a565d43d6a6a13eecaa21f893c89de46e1b7bfe
hash_basis: LF-normalized bytes
---

# THM-2806 -- the fixed-sheet allocation square is central but raw-flat

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

THM-2772 isolated three objects that must not be conflated: a Boolean
bare/source/target/both allocation square, four endpoint translations, and a
nonzero mixed coefficient on one common atom before Fourier allocation.
THM-2791 then retained one literal Boolean ancestry sheet but stopped before
endpoint origin and allocation.

This theorem performs that lift for the most literal possible allocation:
keep the ancestry sheet fixed and either omit or insert each endpoint
carrier.  The resulting answer is exact and sharply split.

```text
positive:  all four central Fourier coefficients are nonzero at one marked
           endpoint origin, with D3=144w;

negative:  the sole primal twist point supporting all four states is flat;
           the central D3 sums 144 bare-only points;

negative:  literal insertion acts by the scalar 1/13 and zero endpoint
           translation;

separate:  moving the sheet through its orbit gives a determinant-one
           endpoint parallelogram, but is not the allocation toggle.       (0)
```

The mechanism is not peculiar to the LRC constants.  It is the tensor square
of the augmentation pair `(constant function, delta at the identity)` for a
finite twist group.  This identifies the exact object that was missing from
the earlier view and also explains why it still does not pay the physical
root/Cech invoice.

## 1. The wandering-selector allocation lemma

Let `G` be a finite abelian group of order `m`, let `R` be a commutative
ring, and
fix `w in R`.  Suppose a marked physical sheet `xi` is wandering for the
carrier orbit:

```text
1_xi C_e=1_xi,
1_xi C_g=0                         for g!=e.             (1)
```

Carrier absence omits the `G` variable, so comparison with a present state
extends its value constantly over that variable.  Carrier presence retains
`1_xi` and therefore has the delta mask in `(1)`.  On `G x G` the four raw
allocation arrays are consequently

```text
B(a,b)=w,
P(a,b)=delta_e(a)w,
Q(a,b)=delta_e(b)w,
H(a,b)=delta_e(a)delta_e(b)w.                           (2)
```

Their pointwise Boolean Mobius face factors as

```text
Omega=B-P-Q+H
     =w(1-delta_e) tensor (1-delta_e).                  (3)
```

Thus `Omega` is supported on `(G minus {e})^2`, where only `B` is present,
and

```text
sum_(a,b) Omega(a,b)=(m-1)^2w.                         (4)
```

The only point at which all four arrays in `(2)` are nonzero is `(e,e)`.
There their vector is `(w,w,w,w)` and their pointwise mixed face is zero.
This proves, before any computation, that a pure idempotent wandering-sheet
allocation cannot supply a nonzero fourfold-co-supported raw face.

After scalar extension to a splitting field whose characteristic does not
divide `m`, use the ordinary one-dimensional character transform of `G`.
The constant function has value `m` at the trivial character and zero
elsewhere, while `delta_e` has value one at every character.  At the central
character the four coefficients and their Segre--Hadamard coordinates are
therefore

```text
v=w(m^2,m,m,1),
D=w((m+1)^2,m^2-1,m^2-1,(m-1)^2).                     (5)
```

In particular `D3=(m-1)^2w`, agreeing with `(4)`, and

```text
D0 D3=D1 D2.                                           (6)
```

At a pair of nontrivial characters only `Hhat` survives.  Hence `(2)--(5)`
are Fourier-dual versions of the same one-corner phenomenon:

```text
primal complement:   only B survives on (G\{e})^2;
dual mixed sector:   only Hhat survives off both trivial characters.       (7)
```

The central address is exceptional only in that all four transformed
marginals are nonzero there.  It does not turn `(e,e)` into a nonflat raw
common atom.

## 2. The retained LRC sheet and the full clock census

Use THM-2782's rail-eight physical cell

```text
(s,t,clock)=(0,4,1).                                    (8)
```

Align the source one-sided carrier `A` and the pulled-back target carrier
`B`, and put `C=A intersect B`.  Exact interval reconstruction contains the
weighted common piece

```text
I=[142004992589460,142005019034340),
weight(I)=27581135604.                                  (9)
```

Its target copy is `I+431933040`.  The unit selector `1_I` and its target
copy lie in all six native and all six adjacent factors, in the order

```text
(E3, clock, q1, q2, c2, c3).                           (10)
```

Both first-failure labels are therefore
`common_after_all_factors`, not an order-dependent earlier seam.  The unit
semantic values are

```text
source carry 12: (0,103478815440),
target carry 6:  (0,103478815440).                     (11)
```

THM-2791's literal pre-marginal sheet has label sets of sizes

```text
966606 and 28534,
966606*28534=27581135604.                              (12)
```

The exact label

```text
(a,b,e')=(59162,26,56658)                              (13)
```

is active on all of `I` and its target copy.  Its six inverse-branch
conditions give the joint chamber

```text
[138281416853580,159555049051860),                     (14)
```

so `(9)` lies strictly inside it.  This is a retained labelled sheet, not
one scalar unit chosen from an aggregate weight.

The all-configuration control reconstructs every

```text
9*9*7=567                                               (15)
```

lawful rail-eight cell.  The source and target exclusive wings are nonempty
in `193` cells each and every opposite-wing intersection is empty.  The
nonempty **common** cells by clock are

```text
clock 0,1,2,3,4,5,6: (0,81,56,56,0,0,0).              (16)
```

Thus this whole fixed configuration bank has common support on only three
of seven clocks.  Even a lawful central marginal extracted below cannot by
itself provide seven physical clock faces.

## 3. Literal absence and presence on the marked sheet

The physical carrier-twist unit is

```text
T/13=22910530602960.                                   (17)
```

Exact half-open interval intersection with the **full** source and target
one-sided carriers gives

```text
I intersect C_0=I,              I intersect C_h=empty for h!=0,
J intersect C'_0=J,             J intersect C'_h=empty for h!=0,           (18)
```

where `J=I+431933040`.  Hence `(1)` applies with `m=13`.

For endpoint representative `ell`, evaluate the omitted-carrier factors
directly:

```text
p_ell=Phi_L(1_I E_ell),          q_ell=Phi_R(1_J E_ell). (19)
```

No carried coefficient is divided to define `(19)`.  Selector idempotence
first proves that `(19)` equals the literal present raw value at twist zero.
Constant extension of the omitted variable then derives the factor thirteen
under unnormalized central Fourier transform:

```text
P_bare(L)=13p(L),                P_present(L)=p(L),
Q_bare(R)=13q(R),               Q_present(R)=q(R).       (20)
```

Equation `(20)` holds on every one of the `169` endpoint origins in each
certified field.  It is a theorem about the fixed selector.  Replacing
`I,J` by each translated sheet in turn is a different orbit construction.

## 4. Raw flatness and the nonzero central marginal

At a fixed endpoint origin put `w=p(L)q(R)`.  Equations `(2)--(5)` specialize
to

```text
(B,P,Q,H)=w(1,delta_(a=0),delta_(b=0),delta_(a=b=0)),

(v00,v10,v01,v11)=w(169,13,13,1),

(D0,D1,D2,D3)=w(196,168,168,144).                      (21)
```

The primal support census is

| state | support in `F_13^2` |
|---|---:|
| `B` | `169` |
| `P` | `13` |
| `Q` | `13` |
| `H` | `1` |
| all four | `1` |
| nonzero raw `Omega` | `144` |

At the sole fourfold point `(0,0)`, the vector is `(w,w,w,w)` and
`Omega=0`.  On the `12*12` complement, the vector is `(w,0,0,0)` and
`Omega=w`.  Consequently

```text
D3=sum_(a,b) Omega(a,b)=144w.                          (22)
```

This is the exact typing boundary.  The four **central transformed**
coefficients share one harmonic address and are all nonzero below, but
their mixed face is assembled from different primal supports.  It is not
THM-2772's required nonzero mixed face on one fourfold-co-supported atom.

At `L=R=(0,0)`, the first exact-order specialization gives

```text
(P0,P1,Q0,Q1)=
(161934023005863791,175075409527953449,
 300649489018742460, 50230041473974177),

v=(346389504894336138,351883238969953710,
   351883238969953710,216790045382338969),

D=(209922877787817004,129599459511997169,
   129599459511997169,211754122479689528)               (23)
```

modulo `352341050142921841`.  The second gives

```text
(P0,P1,Q0,Q1)=
(155560747474790541,11966211344214657,
 879469786637801433,361914377113479889),

v=(17884554354397974,295638590014756546,
   295638590014756546,96307143767239679),

D=(705468878151150745,877931689546517576,
   877931689546517576,479268797051483842)               (24)
```

modulo `956354278959359281`.  Every displayed factor, product, and Hadamard
coordinate is nonzero, and `(6)` holds.  Both fields contain roots of the
required exact order and evaluate the same cyclotomic-integer constructor;
either nonzero image certifies characteristic-zero nonvanishing, while the
second is an independent specialization control.

## 5. The marked representative gauge

The phrase “fixed sheet” means fixed in the marked carrier torsor, not the
same untransported interval in every representative chart.  THM-2763's
gauge is

```text
(ell,a,b)~(ell+W,a+1,b-1).                              (25)
```

The literal selector transports as

```text
source: (ell,I,delta_0) -> (ell+W,T I,delta_1),
target: (ell,J,delta_0) -> (ell+W,T J,delta_12).         (26)
```

Direct evaluation verifies `(26)` coefficientwise at all `169` endpoint
addresses, on both sides and in both fields.  The central harmonic is the
orbit sum, so the cyclic reindexing in `(26)` leaves it invariant.

Thus the construction descends on the object retaining

```text
(marked ancestry sheet, twist origin, endpoint representative).           (27)
```

It does not descend after forgetting that sidecar and holding `I,J` fixed.
Then `27/169` rows fail on each side in each field, beginning at address
`(0,0)` with the nonzero-to-zero hostiles

```text
field 352341050142921841:
  source 254455016269350867 -> 0,
  target 231164267889491750 -> 0;

field 956354278959359281:
  source 318932490657369324 -> 0,
  target 630230755085920022 -> 0.                       (28)
```

The numeric origin in `(23)--(24)` is therefore the distinguished `e_0`
chart of a lawful marked-sheet pullback, not an unmarked quotient current.

## 6. Projective collapse and the endpoint-translation no-go

On the entire endpoint plane, `(20)` gives

```text
P_present=13^(-1)P_bare,
Q_present=13^(-1)Q_bare.                               (29)
```

Exhaustion of every endpoint shift and character proves that the complete
affine covariance sets are exactly

```text
source: (((0,0),(0,0),13^(-1)),),
target: (((0,0),(0,0),13^(-1)),).                      (30)
```

Hence literal allocation has zero source and target endpoint steps.  All
four flags retain

```text
(L,R,q,Delta)=((0,0),(0,0),(0,0),0),                   (31)
```

and their determinant mixed face is zero.

This also explains the coefficient face conceptually.  Dividing the four
central entries by their twist-orbit multiplicities gives

```text
(v00/169,v10/13,v01/13,v11)=(w,w,w,w).                 (32)
```

Their projective cross-ratio is one.  THM-2802's scalar-invariant
bispectrum likewise cannot distinguish bare from present in `(29)`.  The
nonzero additive `D3` is the orbit-cardinality contrast `(13-1)^2w`, not a
projectively nontrivial endpoint interaction.

The four allocation flags form an abstract `V4`, but the intrinsic
projective endpoint observable ties all four vertices after the known
normalizations.  Any tournament orientation on them would therefore come
from an extra label, not from this physical observable.

## 7. The moving-sheet parallelogram is a different object

Now replace the fixed selector by its thirteen translated copies before
forming the endpoint bank.  This moving-orbit bank has full endpoint support.
At

```text
L=R=(0,0),             s=(0,1),             t=(12,0), (33)
```

its four sampled products and their coefficient `D3` are nonzero in both
fields.  Independently, the endpoint labels are

```text
q=(0,0),(0,1),(1,0),(1,1),
Delta=(0,0,0,1),
Delta00-Delta10-Delta01+Delta11=det(s,t)=1.            (34)
```

But the moving bank is not literal presence on `I`.  For every one of the
`169` shifts, all `169` characters, and a globally consistent scalar, the
companion tests

```text
orbit(x)=c chi_u(x) bare(x+s).                          (35)
```

The covariance set is empty on both endpoint sides in both fields.  A
cheaper invariant gives the same no-go: the inverse-address support sizes
are `81` for the literal bare/present banks and `121` for the moving-sheet
bank, while translation, character multiplication, and nonzero scaling
preserve support size.

Thus `(34)` is a genuine determinant-one **orbit sample**, but neither its
endpoint steps nor its nonzero coefficient face are the fixed-sheet
allocation toggles.  The two positive observations in `(21)` and `(34)` are
adjacent, not identified.

## 8. The filtered survivor and the exact next object

The central vector in `(21)` has the relative integral filtration

```text
(v_13(v00/w),v_13(v10/w),v_13(v01/w),v_13(v11/w))
 =(2,1,1,0).                                            (36)
```

Moreover

```text
D3/v11=(13-1)^2=144=1 mod13.                            (37)
```

Naive reduction modulo thirteen sends `(21)` to `(0,0,0,w)`, exactly the
virtual one-corner hostile in THM-2772.  Therefore `(37)` is not yet a
root-deck map: it obtains the tempting unit by collapsing three vertices.

There is nevertheless more structure than in an arbitrary one-corner
vector.  The marked sheet supplies a canonical integral lift and the full
Rees/valuation profile `(36)`.  A two-dimensional associated-graded or
Bockstein sidecar could retain that filtration before extracting the unit
in `(37)`.  Such a construction would have to provide:

```text
source:    the marked rank-one coefficient line and its (2,1,1,0) square;
map:       a gauge-covariant double Bockstein retaining all four lifts;
target:    one F_13 root-deck coefficient;
preserve:  ancestry, clock, endpoint origin, allocation flags, and marker;
test:      distinguish (36) from an arbitrary post-Fourier one-corner.     (38)
```

THM-2771's intrinsic coefficient Bockstein and THM-2807's positive graded
address two-simplex are close analogues, not supplied maps.  In the present
fixed bank, `(16)` also blocks a direct seven-clock realization.  A future
proof must either build the filtered boundary/clock transport in `(38)` or
use a non-idempotent same-atom carrier whose raw fourfold face is already
nonzero.  Another pass through the moving-sheet shortcut is excluded by
`(35)`.

## 8.1 Finite-exact addendum: the tau-twelve simplex is target-only

**FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED ADDENDUM.  This badge applies
only to this subsection and its separately pinned companions; it does not
inherit the theorem's earlier audit badge.**

THM-2807 ended with an ordered six-step allocation test.  Its first item
asked for a nonempty common allocation atom over the explicit base vertex
`n_0`.  On the specified simplex the answer is exactly negative, and the
failure occurs before endpoint origin, Fourier allocation, or the choice of
address.

Keep precisely the rail-eight sheet

```text
(physical clock,sigma,tau)=(1,0,12)                    (38a)
```

with both relative-present safeties, target root one, delayed carry six,
and the three THM-2807 addresses

```text
n_0=3454614,       n_+=3454627,       n_a=4143978.     (38b)
```

As in THM-2782, before integration write

```text
A=source one-sided carrier,
B=pullback(target one-sided carrier),
M=A intersect B,       L=A\B,       R=B\A.             (38c)
```

The primary constructor first forms the fixed-label carrier face, before
the narrow address restriction.  Its interval-piece census is

```text
             A    B    M    L    R
pieces      240  241   0   240  241.                  (38d)
```

In particular

```text
M=empty,                 L=A,                 R=B.      (38e)
```

The definition of these five objects depends on `(38a)` but not on an
address.  Therefore `(38e)` holds uniformly over every later address
cylinder on this fixed sheet, not merely at a sampled centre.

After the two relative safeties, target push, root-one cut, and each narrow
cylinder in `(38b)`, the source-side faces disappear and the target-side
face is one whole weighted interval:

```text
at each n in {n_0,n_+,n_a}:

A=M=L=empty,
B=R=one piece,
mass(R)=60781651775958960/371293,
coefficient(R)=790161473087466480,
weight(R)=27581135604.                                  (38f)
```

The three pieces translate exactly as whole weighted cylinders.  The
independent implementation begins instead with the full one-sided source
and target carriers and inserts the complete semantic section before taking
differences.  In that factor order it obtains

```text
             A    B    M    L    R
pieces        0  241   0     0  241,                  (38g)
```

and independently reproduces every value in `(38f)`.  Equations `(38d)` and
`(38g)` are not competing raw universes: the latter has already imposed the
full semantic section that kills the `240` source pieces in the former.
Their common conclusion is the load-bearing one:

```text
the fully typed tau-twelve simplex is target-only.       (38h)
```

As a positive/hostile control, the two selected right-cofiber tables are

```text
                 n_0       n_+       n_a
tau=3              c         0         c
tau=12             c         c         c,               (38i)

c=(60781651775958960/371293, 790161473087466480).
```

Thus both implementations reproduce THM-2807's positive triangle and its
deleted diagonal vertex while locating the failure strictly in allocation
support, not in the address or coefficient computation.

### The separate common bank is not the simplex

The common-cell clock census `(16)` belongs to the distinct label bank

```text
sigma in {0,1,2,3,8,9,10,11,12},
tau   in {3,4,5,6,7,8,9,10,11}.                         (38j)
```

In particular `tau=12` is outside `(38j)`.  Direct reconstruction of all
`9*9*7=567` cells gives, for both common-carrier and right-cofiber
nonemptiness,

```text
clock 0,1,2,3,4,5,6:       (0,81,56,56,0,0,0).         (38k)
```

More sharply, their nonemptiness indicators agree cell by cell:

```text
(M nonempty,R nonempty)=(false,false) in 374 cells,
                         (true,true)  in 193 cells.      (38l)
```

Nevertheless

```text
M intersect R=empty                         in 567/567 cells. (38m)
```

This is the exact support lesson.  Identical clock support, even
cell-by-cell, does not produce a common atom.  The common carrier and the
right cofiber are two disjoint physical pieces with the same nonemptiness
shadow.

### Consequence for the THM-2803 decoder

The address digits remain

```text
n_0:(7,6),       n_+:(7,7),       n_a:(7,7),
(n_0,n_+,n_a)=(85,98,98) mod169.                        (38n)
```

THM-2803 therefore still proves the conditional statement that a
common-scalar allocation of the first two values would decode the
determinant fibre, and equality of the last two would test high-digit
descent.  But `(38e)--(38h)` show that its hypotheses are not instantiated:
there is no shared physical allocation current or scalar on this simplex.
The coefficient ratio cannot be taken across a missing common atom.

Consequently this addendum closes **only item 1** of THM-2807's proposed
test, negatively, on the fixed simplex `(38a)--(38b)`.  Items `2--6` are not
performed.  The exact next alternatives are:

1. find a positive address simplex inside the genuine common bank `(38j)`;
2. construct a typed boundary map from that common bank to the tau-twelve
   right cofiber, preserving ancestry, clock, endpoint origin, allocation
   flags, and the Abel boundary.

The Rees profile `(36)` and THM-2771 Bockstein remain candidates for the
second route, not supplied maps.

The addendum replays by

```text
python 04-computation/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.py
python -O 04-computation/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.py

python 04-computation/lrc14_tau12_simplex_allocation_support_no_go_independent_audit_thm2806.py
python -O 04-computation/lrc14_tau12_simplex_allocation_support_no_go_independent_audit_thm2806.py
```

Each normal/optimized pair byte-matches its stored output.  The primary
source/output LF hashes are

```text
7fea18161046f8de35b2e6ef04c88a13485de61045bc363f89c5ebfec8f76480
cba545e7beff4fe889c76ae681c47c806969ccd33aa79def54527709a6ffafc6,
```

and the independent source/output hashes are

```text
fd45e0959f27021dcf67a6f22a1be426db3849d4c0a806ba5135e17e5a0d57fe
e18b61016df622283b8bfccb4a565d43d6a6a13eecaa21f893c89de46e1b7bfe.
```

Both scripts have zero Python `assert` nodes.

## 9. Exact evidence and boundary

The primary companion reconstructs the physical geometry, all `567` cells,
the literal ancestry label sets, the raw allocation supports, both endpoint
fields, marked representative gauge, all affine covariances, moving-orbit
hostile, and determinant control.  The independent companion rebuilds the
full one-sided carriers, evaluates absence without importing the primary
result, checks the support-size hostile and all marked gauge rows, and
rederives `(21)--(30)`.

Run

```text
python 04-computation/lrc14_literal_fixed_sheet_allocation_thm2806.py
python -O 04-computation/lrc14_literal_fixed_sheet_allocation_thm2806.py

python 04-computation/lrc14_literal_fixed_sheet_allocation_independent_audit_thm2806.py
python -O 04-computation/lrc14_literal_fixed_sheet_allocation_independent_audit_thm2806.py
```

Both normal/optimized pairs byte-match their stored outputs and the declared
LF hashes.  The two implementations contain no Python `assert` nodes.

Exactly proved are the wandering-selector lemma, the
literal marked LRC sheet, raw-flat/central-nonzero dichotomy, zero-step scalar
allocation law, clock census, and moving-orbit covariance no-go.  Not proved
are a nonzero pre-Fourier common-atom mixed face, an allocation-to-endpoint
translation, the filtered map `(38)`, a target-role-to-root-deck
intertwiner, a seven-clock physical Cech correction, row exclusion, or
LRC(14).

**QED.**
