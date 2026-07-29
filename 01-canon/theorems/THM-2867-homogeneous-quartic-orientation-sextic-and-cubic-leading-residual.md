---
id: THM-2867
title: "Homogeneous quartic orientation sextic and cubic leading residual"
status: >
  PROOF-COMPLETE CANDIDATE; AWAITING EXACT COMPANION AND INDEPENDENT HOSTILE
  AUDIT.  A general quartic has integral edge and oriented-cycle sextics.
  The edge discriminant is
  2^30 A^12 Q^2 Disc(f)^2, while the correctly blown-up orientation
  discriminant is 2^6 Q^12 Jcal^4 Disc(f)^3.  On A=0 the edge primitive
  collapses as (Z^2-B^2)^3 but has an exact cubic matching blow-up; the
  orientation primitive remains a sextic of the three pair differences of
  the surviving cubic, with square class its cubic discriminant.  All extra
  A, Q, and Jcal factors have an exact index interpretation on the stated
  transverse DVR charts, but are not proved Keller ramification or globally
  maximal outside those charts.
source: root/homogeneous-quartic-orientation-boundary-2026-07-28
depends_on:
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
script: 04-computation/quartic_homogeneous_orientation_boundary_thm2867.py
output: 05-knowledge/results/quartic_homogeneous_orientation_boundary_thm2867.out
---

# THM-2867 -- the nonmonic orientation boundary remembers the leading cubic

**PROOF-COMPLETE CANDIDATE; AWAITING EXACT COMPANION AND INDEPENDENT
HOSTILE AUDIT.**

THM-2598 shows that the usual cubic resolvent of a nonmonic quartic has a
universal `1+2` leading shadow.  THM-2758 constructs the six pair-sum/edge
polynomial for a general monic quartic and THM-2769 warns that its square
discriminant does not remove an even Kummer boundary word.  THM-2864 adds
the inequivalent oriented-cycle sextic in a depressed chart.  The natural
question is whether either sextic survives integrally when the quartic
leading coefficient vanishes, and whether that survival transfers the
repo's grade-three anatomy.

The answer is exact but sharply limited:

```text
edge primitive:       collapses, with a removable matching blow-up;
orientation primitive: has a nonzero integral blow-up governed by the
                       leading cubic's pair differences;
Keller consequence:   none without maximal-order/primitive-generation data.
                                                                    (1)
```

Work over a characteristic-zero field.  Let

```text
f(X)=A X^4+B X^3+C X^2+D X+E,                D4=Disc_X(f). (2)
```

Define the integral depression covariants

```text
P=8AC-3B^2,
Q=B^3-4ABC+8A^2D,
R=-3B^4+16AB^2C-64A^2BD+256A^3E.             (3)
```

Over `A !=0`, put `Y=4AX+B`.  Direct expansion gives

```text
g(Y)=256A^3 f((Y-B)/(4A))
    =Y^4+2P Y^2+8Q Y+R,                       (4)
Disc_Y(g)=2^24 A^6 D4.                        (5)
```

The right sides of all formulas below are integral polynomials, so the
identities extend across `A=0` after they are proved in the localization
`Z[A,B,C,D,E,A^(-1)]`.

## 1. The integral edge sextic

Let `x_i` be the roots of `f`, let `y_i=4Ax_i+B`, and, for an edge
`{i,j}`, put

```text
Z_ij=(y_i+y_j)/2=2A(x_i+x_j)+B.               (6)
```

The squared values `v=Z_ij^2`, one for each opposite-edge matching, are
the roots of

```text
T(v)=v^3+Pv^2+Kv-Q^2,                          (7)
K=(P^2-R)/4
 =3B^4-16AB^2C+16A^2(C^2+BD)-64A^3E.          (8)
```

Thus the edge sextic is

```text
E_edge(Z)=T(Z^2)=Z^6+PZ^4+KZ^2-Q^2.           (9)
```

When `A` is a unit, THM-2758's opposite-edge product in monic coordinates
is `-Q/A^3`.  Thus `Q=0` is exactly its pair-sum collision wall in
homogeneous form; `(9)` is the integral nonmonic extension of that theorem.

Applying THM-2864 to `(4)` with

```text
p=2P,                    q=8Q,                    r=R
```

and rescaling its matching variable by four gives

```text
Disc_v(T)=2^12 A^6 D4,                         (10)
Disc_Z(E_edge)=2^30 A^12 Q^2 D4^2
             =D4^2 (2^15 A^6 Q)^2.            (11)
```

Equation `(11)` is the homogeneous form of edge-parity erasure.  Its large
square cofactor is not automatically a branch divisor of the normalized
edge cover.

## 2. The integral orientation blow-up

For an oriented four-cycle `gamma=(i j k l)` define

```text
Omega_gamma=(x_i-x_k)(x_j-x_l)
             ((x_i-x_k)^2-(x_j-x_l)^2),        (12)
W_gamma=A^3 Omega_gamma.                        (13)
```

The raw THM-2864 orientation variable for `g` is

```text
Omega_Y=(4A)^4 Omega_gamma=256A W_gamma.        (14)
```

At first `(13)` appears to divide the raw variable by `A`.  The following
four exact divisibilities show that its monic sextic is nevertheless
integral:

```text
P^2+3R       =64A^2 I,
PR-9Q^2      =16A^2 L,
R^2-P^2R+12PQ^2
              =128A^2 H,
R^3+3P^2R^2-36PQ^2R+108Q^4
              =256A^3 Jcal.                    (15)
```

Here

```text
I=C^2-3BD+12AE,                                (16)

L=128A^2CE-36A^2D^2-48AB^2E+4ABCD
  +3B^3D-B^2C^2,                               (17)

H=512A^4E^2-256A^3BDE-128A^3C^2E+48A^3CD^2
  +160A^2B^2CE+14A^2B^2D^2-16A^2BC^2D
  -30AB^4E-10AB^3CD+4AB^2C^3
  +3B^5D-B^4C^2.                               (18)
```

For later use, `Jcal` is the homogeneous degree-nine polynomial

```text
Jcal=65536A^6E^3
 -49152A^5BDE^2+49152A^5C^2E^2-18432A^5CD^2E+1728A^5D^4
 -24576A^4B^2CE^2+19200A^4B^2D^2E-6144A^4BC^2DE
 +1152A^4BCD^3
 +4608A^3B^4E^2+768A^3B^3CDE-1888A^3B^3D^3
 +1536A^3B^2C^3E-96A^3B^2C^2D^2
 -576A^2B^5DE-960A^2B^4C^2E+696A^2B^4CD^2
 -96A^2B^3C^3D
 +288AB^6CE-63AB^6D^2-48AB^5C^2D+12AB^4C^4
 -27B^8E+9B^7CD-2B^6C^3.                      (19)
```

Equations `(15)` are often the cleaner definition of `(18)--(19)`.

Substitute `(14)` into the orientation sextic of THM-2864 and divide by
`(256A)^6`.  Equations `(15)` cancel every apparent denominator and give

```text
E_or(W)=W^6+2H W^4-Q^2 I L W^2-Q^4 D4.         (20)
```

The homogeneous weights are

```text
wt(P,Q,R,I,L,H,Jcal,D4,W)=2,3,4,2,4,6,9,6,3.  (21)
```

The exact discriminant is

```text
Disc_W(E_or)=2^6 Q^12 Jcal^4 D4^3
            =D4^3 (8Q^6 Jcal^2)^2.             (22)
```

Thus the orientation sextic retains the quartic discriminant square class.
Its `Q` and `Jcal` walls occur with even exponent and can be
primitive-generator or order-index collisions.  They are not proved branch
components of a maximal orientation cover.

## 3. The homogeneous `D8` product

In the matching cubic algebra `K[v]/(T)`, put

```text
c(v)=12v^2+8Pv+P^2-R.                          (23)
```

Rescaling the radicand identity of THM-2864 gives

```text
c(v)^2 v W_gamma^2=2^12 A^4 Q^2 D4.            (24)
```

On an open set where the displayed factors are units, `(24)` gives

```text
[v] [W_gamma^2]=[D4]                           (25)
```

in the square-class group of the matching algebra.  As before,
`W_gamma^2` denotes the orientation radicand lying in that cubic algebra;
`W_gamma` itself lies in its quadratic lift.  Equation `(25)` is the
homogeneous `D8` character product, not an independence assertion.

The factor `A^4` in `(24)` is a square.  Its vanishing at `A=0` warns that
this particular radicand equation is not itself the correct boundary
normalization.  The next two sections give the normalization explicitly.

## 4. The edge collapse has an exact cubic blow-up

Let

```text
delta_3=C^2D^2-4BD^3-4C^3E-27B^2E^2+18BCDE   (26)
```

be the discriminant of the leading cubic

```text
h(X)=B X^3+C X^2+D X+E.
```

On the leading-coefficient divisor,

```text
D4|_(A=0)=B^2 delta_3,                         (27)
T(v)|_(A=0)=(v-B^2)^3,                        (28)
E_edge(Z)|_(A=0)=(Z^2-B^2)^3.                 (29)
```

The triple matching collision in `(28)` has a lossless first blow-up:

```text
T(B^2+At)=A^3 U(t),                            (30)

U(t)=t^3+8Ct^2+16(C^2+BD-4AE)t
    +64(BCD-B^2E-AD^2),                       (31)

Disc_t(U)=2^12 D4.                             (32)
```

At `A=0`, this blow-up is not merely a cubic with the right discriminant.
It is an affine copy of the leading cubic itself:

```text
U_0(t)=-64B^2 h(-(t+4C)/(4B)),                 (32a)
t_i=-4(C+B alpha_i),                           (32b)
```

where `alpha_i` are the roots of `h`.  Thus the matching blow-up recovers
the degree-three source polynomial up to an explicit affine coordinate
change.

Thus, over a discrete valuation ring with `A` a uniformizer and
`B delta_3` a unit, the three matching roots separate after the single
coordinate blow-up `(30)`.  The six edge roots form two three-root clusters
near `+-B`, and

```text
val_A Disc(E_edge)=12.                         (33)
```

This is coordinate clustering at infinity.  Equation `(32)` proves that
the normalized matching cubic need not ramify there.

## 5. The orientation boundary is the cubic difference sextic

Define the other classical leading-cubic invariant

```text
I_0=C^2-3BD,
tau_3=2C^3-9BCD+27B^2E.                        (34)
```

The familiar identity

```text
4I_0^3-tau_3^2=27B^2delta_3                    (35)
```

is included only to identify the invariant boundary; it is universal cubic
algebra, not a Keller equation.

Specializing `(16)--(20)` gives

```text
I=I_0,                   L=-B^2I_0,
H=-B^4I_0,               Jcal=-B^6tau_3,       (36)

E_or(W)|_(A=0)
 =W^6-2B^4I_0W^4+B^8I_0^2W^2-B^14delta_3.     (37)
```

If `alpha_1,alpha_2,alpha_3` are the roots of `h`, then

```text
E_or(W)|_(A=0)
 =prod_(i<j)(W^2-B^6(alpha_i-alpha_j)^2).       (38)
```

Thus the six oriented quartic values do not disappear with the fourth root.
After the precise scaling `(13)`, they limit to the six signed pair
differences of the surviving cubic, multiplied by `B^3`.

The boundary discriminant is

```text
Disc_W(E_or|_(A=0))
 =2^6 B^66 tau_3^4 delta_3^3.                  (39)
```

Consequently `(37)` is squarefree exactly when

```text
B delta_3 tau_3 !=0.                           (40)
```

This is the rigorous `4 -> 3` connection suggested by the quartic
`V4/S3` anatomy.  The orientation lift retains the leading cubic's
discriminant character; the extra `tau_3=0` wall is again a collision of
the chosen primitive values.

There is an exact group-theoretic refinement.  The six factors in `(38)`
are indexed by the ordered pairs

```text
(i,j),                         i != j,             |set|=6. (40a)
```

The generic `S3` action on `(40a)` is regular: an element fixing one ordered
pair fixes all three cubic labels.  A transposition acts as `2^3`, while a
three-cycle acts as `3^2`.  Hence its ambient six-point sign is exactly the
cubic sign character.  The quartic `S4/C4` orientation carrier therefore
degenerates, after the correct infinity blow-up, to the regular `S3`
ordered-pair carrier:

```text
six oriented four-cycles  --A=0/blow-up-->  six directed cubic differences.
                                                               (40b)
```

For a leading cubic with full `S3` Galois group and with `(40)` satisfied,
one directed difference has trivial stabilizer.  It therefore generates the
entire degree-six cubic splitting field, and `(37)` is a primitive polynomial
for that Galois closure.  If the cubic Galois group drops to `C3`, the six
directions split into the two regular `C3` orbits.  This is another
lost-origin boundary.  Equation `(32a)` is the degree-three cubic source;
equation `(38)` is its degree-six Galois closure.  Passing between them
requires exactly the origin/normalization sidecar emphasized by THM-2598.

This is a literal finite `C2/C3` co-occurrence, not a tournament imposed on
six anonymous vertices.  It explains why the odd discriminant class
survives the degree drop.

For comparison, the raw `g`-orientation coordinate `(14)` has all six roots
of valuation one and

```text
val_A Disc(raw orientation sextic)=30          (41)
```

under the unit hypotheses in `(40)`.  Dividing the root coordinate by
`256A` gives `(13)` and removes this entire valuation.  Therefore the
power `A^30` is not intrinsic ramification.

## 6. Two sharp hostiles

The degree-drop boundary can be completely unramified after the blow-up.
Take

```text
(A,B,C,D,E)=(0,1,0,-1,1).                      (42)
```

Then

```text
delta_3=-23,                    tau_3=27,
E_or(W)=W^6-6W^4+9W^2+23,                      (43)
```

which is squarefree by `(39)`, even though the raw edge and orientation
coordinates have collapsed.

Conversely `Jcal=0` can collide the orientation primitive away from
infinity and away from quartic branching.  For

```text
f=X^4+X^2+4X-3                              (44)
```

one has

```text
A=1,                    D4=-22000,            Jcal=0,
E_or(W)=(W^2-1280)^2(W^2+14080).              (45)
```

The quartic is separable, while the orientation sextic has a repeated
primitive value.  Hence `Jcal=0` cannot be called a source branch divisor.

## 7. The local maximal-order and index ledger

A simple transposition of four sheets has cycle type `2^2,1^2` on the six
edges and `2^3` on the six oriented cycles.  The expected tame branch
exponents of the normalized covers are therefore respectively two and
three.  Equations `(11)` and `(22)` have exactly those odd-information
parts:

```text
edge:        D4^2,
orientation: D4^3.                           (46)
```

The remaining factors are squares:

```text
(2^15 A^6Q)^2,              (8Q^6Jcal^2)^2.  (47)
```

On the universal `S4` coefficient cover, or on an `S4` specialization with
the indicated degree-six algebras, these squares have an exact local meaning.
Let `R_v` be an excellent henselian discrete valuation ring with normalized
valuation and with `2` invertible.  Assume the coefficients lie in `R_v`,
the generic quartic splitting algebra has the stated separable `S4`
realization, and the edge and orientation algebras have generic rank six.
Write

```text
O_prim subset O_tilde
```

for a displayed monogenic order and its finite integral closure, and define
its index length to be

```text
length_(R_v)(O_tilde/O_prim).                              (47i)
```

One may equivalently make the inertia calculations after strict
henselization.  The normalized algebras may split locally; they need not be
fields.

### 7.1 The true quartic branch is primitive-maximal

Assume `A,Q,Jcal` are units and

```text
val_v(D4)=1.                                                (47a)
```

Assume in addition that `(47a)` is a simple tame branch with inertia a
transposition.  For a tame permutation cover on a six-set `X`, the exact
normalization formula is

```text
val_v(disc O_tilde)=|X|-|<inertia>\X|.                     (47j)
```

The edge and orientation actions have respectively four and three inertia
orbits, so their normalized degree-six discriminant exponents are

```text
6-4=2,                         6-3=3.                       (47b)
```

Equations `(11)` and `(22)` give the same exponents for the two polynomial
orders.  The order discriminant is the maximal-order discriminant times the
square of the index.  Hence both indices have length zero:

```text
R_v[Z_edge] and R_v[W_orientation] are maximal at (47a).    (47c)
```

Thus the odd `D4^3` in `(22)` is genuine orientation ramification, not a
primitive artifact.

### 7.2 The `Q` and `Jcal` walls are exact index collisions

If `A,D4,Jcal` are units and `val_v(Q)=1`, the normalized associated covers
inside the finite-etale `S4` Galois closure are etale.  Their maximal-order
discriminants are units, while `(11)` and `(22)` have valuations two and
twelve.  Therefore

```text
length(index edge order)=1,
length(index orientation order)=6.                          (47d)
```

If `A,Q,D4` are units and `val_v(Jcal)=1`, the normalized orientation
algebra in the same finite-etale closure is again etale, while `(22)` has
valuation four.  Hence

```text
length(index orientation order)=2.                          (47e)
```

The hostile `(44)--(45)` is the visible special fibre of `(47e)`.

### 7.3 The leading-coefficient costs are removable

At `A=0`, under the unit hypotheses `2Bdelta_3tau_3 !=0`, equation `(32)`
gives an etale matching normalization.  Adjoining `Z` by

```text
Z^2=B^2+At
```

is also etale because `2Z` specializes to `+-2B`.  Thus the normalized edge
cover has unit discriminant, while `(11)` has valuation twelve:

```text
length(index edge order at A=0)=6.                          (47f)
```

The natural orientation polynomial `(20)` has unit discriminant by `(39)`,
so on this unit chart its order is etale and already maximal.  The raw
coordinate `(14)` is `256A` times the natural one.  Scaling the power basis
`1,W,...,W^5` has determinant `(256A)^(1+2+3+4+5)`, and therefore

```text
length(index raw orientation order at A=0)=15.              (47g)
```

This independently recovers the discriminant valuation thirty in `(41)`.

The exact generic ledger is therefore

| divisor | unit side conditions | normalized exponent | polynomial exponent | index length |
|---|---|---:|---:|---:|
| `D4=0` | `A Q Jcal`, simple tame branch | edge `2`; orientation `3` | `2`; `3` | `0`; `0` |
| `Q=0` | `A D4 Jcal` | `0` | edge `2`; orientation `12` | `1`; `6` |
| `Jcal=0` | `A Q D4` | `0` | orientation `4` | `2` |
| `A=0` | `2 B delta_3 tau_3` | `0` | edge `12`; raw orientation `30` | `6`; `15` |
| `A=0`, natural `(13)` | same | `0` | orientation `0` | `0` |

Outside these transverse unit charts, intersections and higher valuations
need their own local calculation.

## 8. The Keller transfer and stopping rule

The index ledger is genuine universal cover information.  It still does not
put the cover on a Keller source.  In particular:

- `(30)--(33)` show that the `A`-valuation can be removed by normalization;
- `(44)--(45)` show that the `Jcal` wall can be only a bad primitive;
- `Q=0` can likewise collapse the chosen edge sign without proving
  ramification of the normalized associated field;
- the sextic fields live in a Galois closure and are not automatically
  polynomial functions on one Keller source sheet.

The next valid Keller question is therefore:

```text
Does a graph-quartic Keller chart identify the intrinsic D4 divisor and
the maximal oriented-cycle order with a source/Jelonek owner divisor? (48)
```

The primitive is already maximal at a simple intrinsic quartic branch by
`(47c)`.  What remains missing is the affine realization: a proof that this
branch/order belongs to the graph chart with the required source ownership,
different, and Jelonek data.  Without that sidecar, the theorem closes no
`A4/S4` branch.

## 9. Exact-companion and status ledger

The exact companion must verify every divisibility in `(15)`, both sextic
formulas, the discriminants `(10)--(11)` and `(22)`, the `D8` product,
the edge blow-up `(30)--(32)`, the cubic factorization `(38)--(39)`, and
both hostiles.  It must use exact arithmetic and explicit exceptions in
normal and optimized mode.

```text
PROVED IN THE CANDIDATE: integral depression and edge sextic;
                         integral orientation blow-up and sextic;
                         exact homogeneous discriminants;
                         homogeneous D8 radicand product;
                         edge matching blow-up at A=0;
                         cubic pair-difference orientation residual;
                         exact A-adic and collision boundaries;
                         transverse maximal-order/index ledger.

NOT PROVED:              global maximal-order generation outside the
                         transverse charts (47a)--(47g);
                         interpretation of square cofactors as ramification;
                         graph-chart or Jelonek owner compatibility;
                         exclusion of A4/S4 Keller monodromy;
                         JC(2), DC(2), G1, or LRC(14).             (49)
```

Promotion requires exact normal/optimized replay with fixed hashes and an
independent audit of every scaling, divisibility, valuation, and boundary
scope.
