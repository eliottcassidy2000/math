---
id: THM-3212
title: "Centered heptic source Morse obstruction and off-center clutch"
status: >
  PROVED + VERIFIED-EXACT.  For every one of the four unmarked S7 covers in
  THM-3123 and every polynomial E_0, the canonical centered source polynomial
  P=V^2 z^4+A_src z+E_0 has five forced, pairwise distinct, reduced Morse
  critical points.  Hence it cannot be one coordinate of a polynomial Keller
  pair.  The obstruction is owned by the squarefree degree-five divisor
  g=ST: both V and A_src vanish on g, while A_src' is a unit there.  In the
  off-center family P=(Vz^2+Bz+C_0)^2+A_src z+E_0, a unit B replaces each
  whole g-fibre by one candidate point, and those five candidates disappear
  exactly when B E_0'-A_src' C_0 is also a unit modulo g.  The control
  B=1,C_0=0,E_0=x passes this local gate only; it is not a Keller mate, a
  Faber-chart entry, or a proof of JC(2).
audit: >
  A self-contained exact companion reconstructs both THM-3123 accessory
  algebras, checks D+A_0=SE^2, the balanced response identities, all
  squarefree/disjointness gates for S,E,T, and the unit resultant of A_src'
  against g.  It then checks the universal centered Hessian calculation and
  the off-center one-jet clutch symbolically.  Ordinary, optimized, and stored
  transcripts are byte-identical; the script has no assert node.
source: root/frontier-synthesis-cont-2026-08-02
depends_on:
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
related:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2822-sextic-response-centered-lift-mod-three-faber-obstruction
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
  - THM-3172-shear-invariant-differential-owner-filtration-and-transverse-recurrence
script: 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
output: 05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out
script_sha256: 03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448
output_sha256: 729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85
hash_basis: LF-normalized bytes
---

# THM-3212 -- centered heptic source Morse obstruction and off-center clutch

**PROVED + VERIFIED-EXACT.**

## 1. Scope and inherited response data

Work over an algebraically closed field of characteristic zero.  THM-3123
classifies two heptic pole passports,

```text
(4,1,1,1),                         (3,2,1,1),              (1)
```

and obtains respectively one and three unmarked `S_7` covers.  For either
passport write

```text
D=x^a(x-1)^b(x^2-ux+v),
T=x(x-1)(x^2-ux+v),

E={a(x-1)(x^2-ux+v)+b x(x^2-ux+v)
     +x(x-1)(2x-u)}/7.                                  (2)
```

The accessory algebra, the linear factor `S`, and the nonzero constant `A_0`
are as follows:

```text
(4,1,1,1):
 q=100u^3+244u^2+237u+44,
 v=(8u^2+9u+8)/7,
 S=x+5(u+1)/7,
 A_0=80v^2(u+1)/343;

(3,2,1,1):
 q=75u^3-89u^2-31u+61,
 v=(24u^2-16u-16)/21,
 S=x+(5u-4)/7,
 A_0=9v^2(5u-4)/343.                                  (3)
```

All identities below are taken in the corresponding reduced algebra
`k[u]/(q)`.  THM-3123 proves

```text
D+A_0=SE^2.                                             (4)
```

Put

```text
C=-7A_0,                 M=SET,
F=SE^2/D,                V=4SDT^2/C^2,
G=CE/(2DT),              A_src=VG=2SET/C.               (5)
```

Then

```text
F=VG^2,                  2VG'+V'G=2,
T_F=A_src^2/V=F.                                          (6)
```

Thus `(5)` is the canonical polynomial coefficient suggested by the balanced
response, not an arbitrary choice of source lift.  What THM-3123 did **not**
decide was whether this abstract response enters a two-dimensional polynomial
Keller chart.

## 2. The degree-five owner divisor

Define

```text
g=ST.                                                    (7)
```

For both accessory algebras, `g` is squarefree of degree five and

```text
V/g=4DT/C^2,                     A_src/g=2E/C.            (8)
```

Moreover

```text
gcd(g,A_src')=1.                                         (9)
```

Here is the complete gate behind `(9)`.  The discriminants of `g,T,E` and the
three pairwise resultants of `S,T,E` are units in `k[u]/(q)`.  Hence the five
roots of `ST` are simple, and none is a root of `E`.  At a root `alpha` of
`g`, differentiating `A_src=2gE/C` gives

```text
A_src'(alpha)=2g'(alpha)E(alpha)/C !=0.                  (10)
```

This proves `(9)` and also identifies the missing coordinate: the obstruction
is not merely common support of `V` and `A_src`; it is that their common
degree-five divisor has a transverse `A_src` one-jet.

## 3. Five forced Morse critical points

For an arbitrary polynomial `E_0(x)` consider the centered source family

```text
P_(E_0)(x,z)=V(x)^2 z^4+A_src(x)z+E_0(x).                (11)
```

If `alpha` is any of the five roots of `g`, then `V(alpha)=A_src(alpha)=0`.
Consequently

```text
P_z(alpha,z)=0,
P_x(alpha,z)=A_src'(alpha)z+E_0'(alpha).                 (12)
```

Equation `(10)` makes the second expression affine with nonzero slope.  Thus
there is exactly one critical point over each root of `g`, namely

```text
(alpha, -E_0'(alpha)/A_src'(alpha)).                     (13)
```

The five points are distinct because their `x`-coordinates are distinct.  At
each point,

```text
P_zz=0,                    P_xz=A_src'(alpha),
det Hess(P)=-A_src'(alpha)^2 !=0.                         (14)
```

They are therefore reduced Morse critical points for **every** choice of
`E_0`.  If a polynomial `Q` satisfied `Jac(P_(E_0),Q)=kappa in k*`, the
gradient of `P_(E_0)` could vanish nowhere.  The points `(13)` contradict
that necessary condition.  Hence no member of `(11)` can be a coordinate of
a polynomial Keller pair.

In the Faber notation of THM-3133, this centered lane has

```text
d_F=0,                          s_F=-E_0.                 (15)
```

Allowing arbitrary `E_0` therefore does not repair the obstruction.

## 4. The exact off-center clutch

Now retain the same response coefficients `V,A_src` but allow

```text
P=(Vz^2+Bz+C_0)^2+A_src z+E_0,                           (16)
```

where `B,C_0,E_0 in k[x]`.  At a root `alpha` of `g`, write
`H=Bz+C_0`.  Since `V=A_src=0` there,

```text
P_z(alpha,z)=2B(alpha)(B(alpha)z+C_0(alpha)).             (17)
```

Suppose `B` is a unit modulo `g`.  The unique zero of `(17)` is

```text
z_alpha=-C_0(alpha)/B(alpha).                            (18)
```

At this point `H=0`, so the square in `(16)` contributes nothing to `P_x`
and

```text
P_x(alpha,z_alpha)
 ={B E_0'-A_src'C_0}(alpha)/B(alpha).                    (19)
```

Therefore, inside the unit-`B` lane, all five forced `g`-fibre critical
points disappear **if and only if**

```text
Delta_B=B E_0'-A_src'C_0
```

is a unit modulo `g`.  This is a source, target, map, and sidecar statement:
the centered source `(11)` is the source; `(16)` is the off-center target;
adjoining `B,C_0` is the map; `T_F=A_src^2/V=F` is preserved; the centered
coordinates `d_F,s_F` are changed; and the cheapest decisive sidecar is the
pair of residue classes `(B,Delta_B)` in `k[x]/(g)`.

The smallest local control is

```text
B=1,                    C_0=0,                    E_0=x. (20)
```

Then `Delta_B=1`; at every root of `g`, `P_z=2z` and its only zero `z=0`
has `P_x=1`.  Its Faber coordinates are

```text
d_F=-1/(4V),                     s_F=G/2-x.              (21)
```

Thus `(20)` genuinely removes the five-point local obstruction but moves the
candidate off the centered Faber lane.  It may have critical points away from
`g`, need not satisfy the remaining Faber fluxes, and supplies neither a
polynomial mate nor global chart entry.  This is a precise next test for the
owner filtration of THM-3172, not a proof of `JC(2)` or `DC(2)`.

## Exact reproduction

Run

```text
python3 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
python3 -O 04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py
```

Both modes must reproduce the declared output byte for byte.  The companion
uses exact rational, polynomial, resultant, gcd, and accessory-algebra
arithmetic only.  It imports no repository executable and uses explicit
`require` gates, so optimized execution retains every check.  QED.
