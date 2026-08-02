---
id: THM-3059
title: "Quartic two-jet even-Jelonek C3 escape counterexample and parity-inertia law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Reversing a primitive
  quartic at a leading-coefficient divisor makes the infinity branches an
  integral quartic order.  In the tame strict-henselian setting, the pole-
  clearing exponent of the monic discriminant is odd exactly when local
  inertia is an odd permutation; order index and leading-coefficient taxes
  are even and cannot change this parity.  The dominant field-degree-four
  two-jet map (x,xz^2+y,xyz^2+z) has generic Galois group S4, exact Jelonek
  set {u=0}, C3 inertia there, and cleared discriminant exponent 8.  It
  belongs to an exact two-parameter S4 family whose Newton wall a=2b
  separates ternary C3 escape from binary transposition escape, and it has a
  simpler planar S4 shadow (x,x^m y^4+y) with cleared exponent 4m.  This
  refutes HYP-9027's general-dominant odd-exponent clause.  The nonconstant
  Jacobians put the separate odd discriminant factors on ordinary critical
  branches, so the Keller-restricted odd-inertia question and JC remain open.
source: codex-jc-resolvent-bridge-2026-08-01
audit: >
  root/2026-08-01 independently rederived quartic and resolvent
  irreducibility, generic S4, the discriminant minima, reciprocal Newton
  wall, tame inertia cases, parity law, exact Jelonek sets, planar shadow,
  and branchwise Jacobian residues.  Ordinary, optimized, and stored exact
  transcripts byte-match.  ACCEPT.
depends_on:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
related:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/quartic_twojet_even_jelonek_c3_escape_thm3059.py
output: 05-knowledge/results/quartic_twojet_even_jelonek_c3_escape_thm3059.out
script_sha256: 16a1d4813950bc623a48f99e0d140f5aff69b465ff872378cd14815084abb694
output_sha256: af579b31fd427579cfcf6797e9d10a7c794646d8452e6a65d60ba7809fd50b26
hash_basis: LF-normalized bytes
---

# THM-3059 -- odd discriminant clearing is odd inertia, not nonproperness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and result

[THM-3057](THM-3057-tame-quartic-inertia-clutch-index-resonance.md)
shows that a tame integral quartic order has

```text
v(Disc)=4-number_of_inertia_orbits+2*order_index.          (1)
```

Its sharp hostile also shows why integer-normalized matching data must not be
used without the base valuation.  The missing operation for a Jelonek
leading-coefficient divisor is **reciprocal reversal**: roots escaping in the
original coordinate become integral roots near zero.  Reversal turns (1)
into an exact parity test for the pole-clearing exponent in
[HYP-9027](../../05-knowledge/hypotheses/HYP-9027-twojet-disc-jelonek-odd-exponent-law.md).

That test exposes a counterexample to the hypothesis's current
general-dominant formulation.  Define

```text
F:C^3 -> C^3,
F(x,y,z)=(x, xz^2+y, xyz^2+z).                            (2)
```

This is a two-jet map in `z`, is dominant of field degree four, and has
generic Galois group `S4`.  Its Jelonek set is exactly the plane `{u=0}`.
For the primitive quartic fiber polynomial `N`, however,

```text
leading_coefficient(N)=u^2,
v_u(Disc N)=4,
6 v_u(u^2)-v_u(Disc N)=12-4=8.                           (3)
```

Thus its unique Jelonek component has **even**, not odd, cleared exponent.
The local mechanism is one finite sheet plus a three-cycle of escaping
sheets.  The example is not Keller: its separate odd discriminant factor is
the image of the ordinary critical locus.  Consequently (2) refutes the
general-dominant clause, not a Keller-restricted parity law.

## 2. General reciprocal parity--inertia lemma

Let `R` be a strictly henselian DVR with fraction field `K`, normalized
valuation `v`, and residue characteristic different from `2` and `3`.  Let

```text
N(T)=a_4T^4+a_3T^3+a_2T^2+a_1T+a_0 in R[T]               (4)
```

be primitive and separable over `K`, with `a_0` a unit and `a_4` nonzero.
Assume the splitting algebra is tame.  Put

```text
rho=v(a_4),
N^vee(X)=X^4 N(1/X),
g(X)=a_0^(-1)N^vee(X).                                   (5)
```

Then `g` is monic and integral.  Let `i` be the length of its monogenic
order inside the normalization of its total quotient algebra, and let
`sigma` generate tame inertia on its four roots.  Define

```text
d_sigma=4-number_of_orbits(sigma).                        (6)
```

If `q=N/a_4` is the monic polynomial in the original coordinate and

```text
E=6rho-v(Disc N)=-v(Disc q),                              (7)
```

then

```text
v(Disc N)=d_sigma+2i,
E=6rho-d_sigma-2i,
E mod 2=d_sigma mod 2.                                   (8)
```

Equivalently,

```text
E is odd  iff  sigma is an odd permutation.               (9)
```

The equality in (7) is an identity even when `E` is not positive; the
"pole-clearing" terminology is used only when `Disc(q)` has a pole.

To prove (8), first note the exact reversal identity

```text
Disc(N^vee)=Disc(N).                                      (10)
```

It follows either from the homogeneous binary-quartic discriminant or by
writing the roots of `N^vee` as the reciprocals of those of `N`; every root
and leading-coefficient factor cancels.  Dividing by the unit `a_0` does not
change the valuation.  Over a strict-henselian base, each tame inertia orbit
of length `ell` contributes `ell-1` to the maximal-order discriminant.
Summing gives `d_sigma`.  The order-discriminant formula contributes the
square of the index, hence `2i`, proving the first equation of (8).

For a quartic, scaling `N` by `a_4^(-1)` scales its discriminant by
`a_4^(-6)`.  This proves (7) and the second equation of (8).  Both `6rho`
and `2i` are even.  Finally the sign of a permutation with `r` cycles on
four letters is `(-1)^(4-r)`, which proves (9).

This gives a useful exact replacement for the false slogan
"Jelonek means odd":

```text
leading drop + reciprocal integral order
    -> cleared-exponent parity = infinity-inertia sign.   (11)
```

Nonproperness by itself does not choose the sign.

## 3. Fiber quartic and generic `S4`

Write target coordinates as `(u,v,w)`.  From (2),

```text
x=u,
y=v-uz^2,
w=uvz^2-u^2z^4+z.                                      (12)
```

Hence `z` obeys the primitive quartic

```text
N(T)=u^2T^4-uvT^2-T+w.                                  (13)
```

Over `C(u,v)`, equation (13) is the generic fiber of the polynomial map

```text
T |-> -u^2T^4+uvT^2+T,                                  (14)
```

which has degree four.  Therefore

```text
[C(u,v,T):C(u,v,N(T)-w)]=4,                              (15)
```

so (13) is irreducible over `C(u,v,w)` and `F` has field degree four.
This also proves dominance without appealing merely to a nonzero Jacobian.

The monic depressed quartic has

```text
p=-v/u,                 q=-1/u^2,              r=w/u^2. (16)
```

Its standard matching resolvent becomes, after clearing denominators,

```text
R(W)=u^4W^3-2u^3vW^2+u^2(v^2-4w)W-1.                   (17)
```

Put `X=uW`.  Solving (17) for `w` gives

```text
w=[uX(X-v)^2-1]/(4uX).                                  (18)
```

The numerator and denominator in (18) are coprime and the rational map has
degree three.  Thus the cubic resolvent is irreducible over `C(u,v,w)`.

Direct calculation gives

```text
Disc_T(N)=u^4 H,                                         (19)

H=16u^2v^4w-128u^2v^2w^2+256u^2w^3
  +4uv^3-144uvw-27.                                     (20)
```

As a rational function of `w`, `H` has an odd pole of order three at
`w=infinity`; it is therefore not a square.  The quartic is irreducible, its
resolvent is irreducible, and its discriminant is nonsquare.  The standard
quartic Galois criterion now gives

```text
Gal(N/C(u,v,w))=S4.                                      (21)
```

There is also a one-line arithmetic control.  At `(u,v,w)=(1,0,1)`, the
quartic is `T^4-T+1`, irreducible modulo two; its cubic resolvent is
`W^3-4W-1`, irreducible by the rational-root test; and its discriminant is
the nonsquare prime `229`.  Thus an explicit `S4` specialization independently
checks that the generic argument has not selected the wrong transitive row.

This makes (2) a stronger hostile than a cyclic or reducible cover: even a
generic `S4` two-jet may have an even Jelonek component.

## 4. The Jelonek set is exactly one plane

Over `u!=0`, divide (13) by `u^2`.  It becomes a monic equation for `z`
over

```text
C[u,u^(-1),v,w],                                         (22)
```

and (12) then recovers `y`.  Thus the restriction of `F` over `{u!=0}` is
finite, so its nonproper-value set is contained in `{u=0}`.

Conversely fix any `(v_0,w_0) in C^2`, put `u=s^3`, and introduce
`z=s^(-2)Z`.  Equation (13), multiplied by `s^2`, becomes

```text
G_s(Z)=Z^4-sv_0Z^2-Z+w_0s^2=0.                          (23)
```

At `s=0`,

```text
G_0(Z)=Z(Z^3-1).                                        (24)
```

Each cube root `zeta` of one is simple because
`partial_Z G_0(zeta)=3`.  The analytic implicit-function theorem therefore
gives a convergent branch `Z(s)` with `Z(0)=zeta`.  Set

```text
x=s^3,
z=s^(-2)Z(s),
y=v_0-s^(-1)Z(s)^2.                                     (25)
```

Then (23) and (12) give exactly

```text
F(x,y,z)=(s^3,v_0,w_0).                                  (26)
```

Both `y` and `z` escape while (26) tends to `(0,v_0,w_0)`.  Every point of
the plane is a nonproper value, and hence

```text
A_F={u=0}.                                                (27)
```

The root `Z=0` in (24) is the one finite branch; the other three are the
escaping branches.

## 5. Exact local anatomy and the even exponent

Reverse (13):

```text
N^vee(X)=wX^4-X^3-uvX^2+u^2.                             (28)
```

At the generic point of `{u=0}`, both `v` and `w` are units.  The lower
Newton polygon has vertices

```text
(0,2), (3,0), (4,0),                                     (29)
```

while `(2,1)` lies strictly above the first edge.  Its slopes are `-2/3`
and `0`.  After `u=s^3`, equation (23) splits the first face into three
simple branches.  The tame inertia action `s -> omega s` cycles them and
fixes the unit branch.  Thus

```text
sigma has type (123),                 d_sigma=2.          (30)
```

Equations (19)--(20) give `v_u(Disc N)=4`, since `H=-27 mod u`.  Applying
(8) to (28) yields

```text
4=2+2i,                             i=1.                  (31)
```

In the splitting valuation normalized by `w_L(u)=3`, the three small-root
edges have valuation two and the three edges to the fixed unit root have
valuation zero.  Hence THM-3057's matching vector is

```text
lambda=(2,2,2)=(1+i,1+i,1+i).                            (32)
```

This is precisely the index-one `C3` cone, not a hidden transposition.

Finally, the leading coefficient in (13) has order two.  Therefore

```text
Disc(N/u^2)=H/u^8,
E=6*2-4=8.                                               (33)
```

The equality `(33)` is the promised even Jelonek exponent.

## 6. An infinite binary--ternary escape family

The hostile is not isolated.  For integers `a>=1`, `b>=0`, put

```text
F_(a,b)=(x, x^a z^2+y, x^b y z^2+z).                    (33a)
```

Exactly as in (12)--(15), this is a dominant two-jet map of field degree
four with primitive fiber polynomial

```text
N_(a,b)=u^(a+b)T^4-u^b vT^2-T+w.                         (33b)
```

It is finite over `u!=0`.  At `u=0`, the specialization is `-T+w` when
`b>0`, of degree one, and `-vT^2-T+w` when `b=0`, of degree at most two.
For every fixed `(v,w)`, projective root specialization therefore puts at
least two roots at infinity; Newton--Puiseux continuation supplies Laurent
branches above the exact target curve `(u,v,w)=(s^e,v,w)`.  Hence every point
of the plane is nonproper and the Jelonek set is again exactly `{u=0}`.

The same degree-four polynomial-map argument proves quartic irreducibility.
For the monic quartic, the standard resolvent equation solves as

```text
w=u^(a+b)(W-vu^(-a))^2/4-1/(4W u^(a+b)).                 (33b1)
```

After putting the right side over a common denominator, its numerator has
degree three in `W`; the nonzero constant term prevents a common `W` factor.
Thus this is a reduced degree-three rational map and the resolvent is
irreducible.  Together with the odd pole of the discriminant at `w=infinity`,
this proves that the generic Galois group is `S4` for every `(a,b)`.

The full discriminant is

```text
Disc N_(a,b)
 =256u^(3a+3b)w^3-128u^(2a+4b)v^2w^2
  -144u^(2a+3b)vw-27u^(2a+2b)
  +16u^(a+5b)v^4w+4u^(a+4b)v^3.                         (33c)
```

Only the last terms in the two competing groups can be initial.  Hence,
at the generic point of the Jelonek plane,

```text
delta_(a,b)=v_u(Disc N_(a,b))
           =min(2a+2b,a+4b),                            (33d)

E_(a,b)=6(a+b)-delta_(a,b)
 = 4(a+b),             a<=2b;
 = 5a+2b,              a>=2b.                           (33e)
```

The formulas agree when `a=2b`.  The wall has a direct Newton meaning.  The
reciprocal polynomial has points

```text
(0,a+b), (2,b), (3,0), (4,0).                            (33f)
```

- If `a<2b`, the first face has length three and slope `-(a+b)/3`.
  Its inertia is `C3` when `3` does not divide `a+b`, and is trivial when
  the slope is integral.  Formula (33e) is always even.
- If `a>2b`, the first face has length two and slope `-a/2`, followed by
  the integral singleton slope `-b`.  Its inertia is a transposition exactly
  when `a` is odd.  Formula (33e) has exactly the same parity.
- If `a=2b`, all three lost roots lie on an integral-slope face; its generic
  residual cubic is separable over the strict-henselian residue field, so
  inertia is trivial and (33e) is even.

Thus `a=2b` is an exact phase boundary between the repo's recurring ternary
and binary tower grammars.  It does not identify those grammars abstractly:
the common object is the reciprocal Newton polygon, and the preserved
predicate is precisely inertia sign.  The original example is the first
ternary cell `(a,b)=(1,1)`.

There is a sharper unit-level stopping boundary.  For the whole family,

```text
J_(a,b)=1+2x^b y z-2x^(a+b)z^3,
(partial_T N_(a,b))|_F=-J_(a,b).                         (33g)
```

In the ternary regime `a<2b`, put `u=s^3` and
`z=s^(-(a+b))Z`.  After multiplying the fiber equation by `s^(a+b)`,
the escaping face is

```text
Z^4-s^(2b-a)vZ^2-Z+ws^(a+b)=0.                          (33h)
```

At an escaping root `Z^3 -> 1`, while

```text
partial_T N_(a,b)
 =4Z^3-2s^(2b-a)vZ-1 -> 3.                              (33i)
```

The finite root has `T->w` and derivative `-1`.  Hence the Jacobian residues
on the four branches are

```text
finite branch: 1;                 C3 orbit: -3,-3,-3.    (33j)
```

All four Jacobian valuations are zero.  Thus a valuation-only or matching-
clutch test cannot see why this family is non-Keller; the obstruction is the
incompatible **unit residues** in (33j).  In a general primitive-element
chart the chain rule inserts the cofactor

```text
det(partial(u,v,z)/partial(x,y,z))                       (33k)
```

between `partial_T N` and the physical Jacobian.  Any attempt to exclude a
Keller `C3` component must therefore retain that cofactor and its residue on
each normalization branch.  This is a more precise sidecar than another
discriminant valuation.

### 6.1 The planar quartic shadow

The same object occurs before any three-variable two-jet presentation.  For
`m>=1`, define the planar polynomial map

```text
G_m:A^2 -> A^2,             G_m(x,y)=(x,x^m y^4+y).       (33l)
```

Its primitive fiber quartic is

```text
P_m(T)=u^mT^4+T-v,                                       (33m)
```

and direct calculation gives

```text
Disc(P_m)=-u^(2m)(27+256u^m v^3),
Disc(P_m/u^m)=-(27+256u^m v^3)/u^(4m).                  (33n)
```

Thus the primitive discriminant order is `2m` and the cleared exponent is
`4m`, always even.  The map is finite over `u!=0`; at `u=0` its fiber degree
drops from four to one, and the same Puiseux argument proves that its Jelonek
set is exactly the line `{u=0}`.

Irreducibility follows from the degree-four polynomial map
`T -> u^mT^4+T`.  The cubic resolvent is irreducible because solving it gives

```text
v=(1-u^(2m)W^3)/(4u^mW),                                 (33n1)
```

a reduced rational function of degree three.  The second factor in (33n)
has odd degree three in `v`, hence is nonsquare.  Therefore the generic
Galois group is again `S4`.  At `m=1`, the arithmetic control

```text
T^4+T-1,        W^3+4W-1,        Disc=-283              (33o)
```

has irreducible quartic modulo two, irreducible cubic by the rational-root
test, and nonsquare discriminant.

The reciprocal Newton points are `(0,m),(3,0),(4,0)`.  Hence the escaping
inertia is `C3` exactly when `3` does not divide `m`, and is trivial when the
slope is integral.  Finally,

```text
J_(G_m)=1+4x^m y^3=(partial_T P_m)|_(u,v)=G_m(x,y).      (33p)
```

The finite branch again has Jacobian residue `1`, while the three escaping
branches have residue `-3`.  This planar family proves that the
parity/inertia and unit-residue phenomena belong to the quartic inverse
object itself.  They are not artifacts of suspension or of the vector
two-jet frame.

## 7. Where the odd discriminant went

The Jacobian of (2) is

```text
J_F=1+2xyz-2x^2z^3,                                      (34)
```

and direct substitution in (13) gives

```text
(partial_T N)|_(u,v,w)=F(x,y,z)=-J_F.                    (35)
```

Thus, away from `u=0`, the remaining factor `{H=0}` in (19) is the ordinary
critical discriminant of the finite map.  It is separate from the Jelonek
plane because

```text
H mod u=-27.                                             (36)
```

The complete mechanism is therefore

```text
u=0:       nonproper C3 infinity inertia, primitive order 4 / cleared pole 8;
H=0:       finite critical branch carrying odd S4 sign ramification.      (37)
```

The Jacobian in (34) is nonconstant.  A Keller map cannot use this critical
branch to carry its sign character.  This explains exactly why the example
does not settle the Keller-restricted question, while also showing that
dominance, field degree four, two-jet form, and generic `S4` are insufficient.

## 8. Correction and surviving frontier

The map (2) refutes clause 2 of HYP-9027 with its present quantifier over
general dominant two-jet maps.  Clause 1 survives this hostile because
`{u=0}` divides the leading coefficient `u^2`.  The strongest uniform
replacement is the proved equivalence (9).

For a hypothetical degree-four Keller map, the repaired question is:

```text
must every Jelonek component have odd infinity inertia?                 (38)
```

THM-2633 guarantees a fixed sheet but permits both a transposition and a
three-cycle.  THM-2465 and purity force an odd component somewhere in an
`S4` lane, but not yet at every component.  Excluding the local `C3` escape
pattern (29), or proving that it must meet a finite critical divisor, would
establish the uniform Keller version.  Neither statement is proved here.

The exact connection ledger is

```text
source:       primitive nonmonic quartic at a leading divisor;
map:          reciprocal reversal -> integral local order;
preserved:    discriminant, inertia sign, order-index parity;
destroyed:    affine ownership and the distinction critical/nonproper;
sidecar:      original-source Jacobian and boundary regularity;
hostile:      F=(x,xz^2+y,xyz^2+z).                       (39)
```

```text
PROVED HERE:       reciprocal parity--inertia formula;
                   explicit dominant degree-four two-jet map;
                   irreducible quartic and resolvent, generic S4;
                   exact Jelonek plane and C3 escape;
                   even cleared exponent 8;
                   separation from the ordinary critical branch.

REFUTED HERE:      HYP-9027 clause 2 for general dominant two-jet maps.

NOT PROVED:        existence or exclusion of a degree-four two-jet Keller map;
                   uniform odd inertia for Keller Jelonek components;
                   exclusion of A4, S4, G1, JC(2), or DC(2).              (40)
```

## 9. Exact companion

Run

```text
python3 04-computation/quartic_twojet_even_jelonek_c3_escape_thm3059.py
python3 -O 04-computation/quartic_twojet_even_jelonek_c3_escape_thm3059.py
```

Both modes must LF-byte-match the stored transcript.  The companion verifies
the eliminant, Jacobian derivative, primitive and monic discriminants,
quartic-resolvent discriminant equality, rational-map degree certificates,
reciprocal polynomial and Newton data, exact Puiseux chart, index invoice,
all five tame cycle-type parity rows, and the discriminant-order formula on
the exact grid `1<=a<=5`, `0<=b<=5`.  It also checks the planar formulas for
`1<=m<=6` and the arithmetic `S4` specialization (33o).  Every truth-bearing
check uses an explicit runtime exception rather than a Python assertion.
