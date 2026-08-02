---
id: THM-3123
title: "Heptic e3 remaining accessory classification and S7 monodromy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the normalized
  maximal-pole response layer N=7,e=3,s=1,h=4,
  the remaining pole multisets (4,1,1,1) and (3,2,1,1) have radical,
  already-saturated accessory algebras of length six.  They give respectively
  one and three unmarked covers, all with monodromy S7.  Together with
  THM-2840 this would complete the abstract heptic accessory layer, but it
  does not supply Keller-chart entry, the missing Faber fluxes, JC(2), or
  DC(2).
audit: >
  An independent algebra path derived R'/E from 7D/T-E-2SE', recovered both
  equation pairs and cubic eliminants, rebuilt the two triangular quotient
  rings, and checked radicality, every saturation unit, the response
  identities, and the Jacobian.  A separate 5040-permutation census recovered
  the 105 involutions, 35 planar rows, 7/21/7 passport histogram, 1/3/1
  unmarked orbit counts, and monodromy orders 5040/5040/5040/5040/14.
  Normal, optimized, and stored exact transcripts agree.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2816-maximal-pole-clean-nielsen-ribbon-tree-prufer-vandermonde-count
  - THM-2840-heptic-e3-2221-chebyshev-accessory-classification
related:
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
script: 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py
output: 05-knowledge/results/jc_heptic_e3_remaining_accessory_classification_thm3123.out
script_sha256: 12e6d3ae5e0ec253d79eb127987150d59819685db633f4c147fb7058f01accc2
output_sha256: 42077492ad769c1c09eeb0af6fb83a071898c043ff1a8cc563146b80239c483a
hash_basis: LF-normalized bytes
---

# THM-3123 -- the remaining heptic accessory strata

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Scope and inheritance

Work over an algebraically closed field of characteristic zero.  In the
balanced response notation of THM-2796 fix

```text
N=7,                 e=3,                 s=1,                 h=4.
```

The zero inertia is `(2,2,2,1)` and the third inertia is the full cycle
`(7)`.  THM-2840 classifies the pole multiset `(2,2,2,1)` as the single
unmarked Chebyshev cover.  The two remaining unordered pole multisets are

```text
(4,1,1,1),                 (3,2,1,1).                 (1)
```

The closest mechanism is THM-2840's equal-critical-value remainder ideal.
The necessary sidecar is the marked Nielsen class: THM-2796 already shows
that from `e=2` onward an integer passport need not determine a map, and the
present `(3,2,1,1)` passport has three unmarked classes.

This theorem concerns only the normalized one-variable response layer.  It
does not show that a planar Keller pair enters the polynomial exact-prefix
chart or satisfies the remaining Faber flux equations.

## 2. The universal two-equation mechanism

For one of the canonical ordered pole passports

```text
p=(a,b,1,1)=(4,1,1,1) or (3,2,1,1),
```

normalize the pole points to `(0,1,lambda,mu)` and put

```text
D=x^a(x-1)^b(x-lambda)(x-mu),
T=x(x-1)(x-lambda)(x-mu),                              (2)

K=
 a(x-1)(x-lambda)(x-mu)
 +b x(x-lambda)(x-mu)
 +  x(x-1)(x-mu)
 +  x(x-1)(x-lambda),

E=K/7.
```

Then `E` is a monic cubic and

```text
D'=7(D/T)E.                                             (3)
```

Divide by the prospective double-zero square:

```text
D=S E^2+R,             deg(S)=1, S monic, deg(R)<6.    (4)
```

Differentiating `(4)` and using `(3)` shows `E|R'`.  Since `deg(R')<=4`,
the quotient `R'/E` is linear.  Consequently

```text
R is constant
  iff both coefficients of R'/E vanish
  iff D+A=S E^2 for the nonzero constant A=-R.          (5)
```

Thus the whole heptic accessory condition is exactly two equations.  This
is an equivalence, not a critical-value necessary condition followed by an
unproved reconstruction.

Both simple-pole variables enter symmetrically.  Set

```text
u=lambda+mu,                  v=lambda mu.              (6)
```

## 3. Pole multiset `(4,1,1,1)`

For `p=(4,1,1,1)`, the two coefficients in `(5)`, after removal of nonzero
rational scalars, are

```text
8u^2+9u-7v+8=0,
25u^2+25uv+25u+11v=0.                                  (7)
```

Equivalently,

```text
v=(8u^2+9u+8)/7,
q_4111(u)=100u^3+244u^2+237u+44=0.                     (8)
```

The cubic is separable:

```text
Disc(q_4111)=-480200000 !=0.                            (9)
```

The quotient and accessory constant are

```text
S=x+5(u+1)/7,
A=80v^2(u+1)/343.                                     (10)
```

Eliminating `mu` instead gives the triangular ordered-chart basis

```text
51700mu
 +157500lambda^5+446800lambda^4+727623lambda^3
 +734179lambda^2+550948lambda+236700=0,

h_4111(lambda)=
 2500lambda^6+6100lambda^5+10121lambda^4+12558lambda^3
 +10121lambda^2+6100lambda+2500=0.                    (11)
```

Hence the ordered accessory algebra has length six.

## 4. Pole multiset `(3,2,1,1)`

For `p=(3,2,1,1)`, the two equations are

```text
24u^2-16u-21v-16=0,
40u^2+50uv-32u-61v=0.                                  (12)
```

Equivalently,

```text
v=(24u^2-16u-16)/21,
q_3211(u)=75u^3-89u^2-31u+61=0.                       (13)
```

Again the cubic is separable:

```text
Disc(q_3211)=-149361408 !=0.                            (14)
```

Here

```text
S=x+(5u-4)/7,
A=9v^2(5u-4)/343.                                     (15)
```

The ordered triangular basis is

```text
616875lambda^5-1459525lambda^4-1163843lambda^3
 +1256649lambda^2+1397860lambda-261792mu+844544=0,

h_3211(lambda)=
 1875lambda^6-2225lambda^5-1967lambda^4-443lambda^3
 +2072lambda^2+1472lambda+512=0.                      (16)
```

This ordered accessory algebra also has length six.

## 5. Radicality and exact saturation

For either cubic, a root `u` determines `v` by `(8)` or `(13)`, and
`lambda,mu` are the roots of

```text
z^2-uz+v.                                               (17)
```

The following gate polynomials have nonzero resultants with the relevant
cubic.  Harmless nonzero rational scalars have been cleared.

| gate | `(4,1,1,1)` resultant | `(3,2,1,1)` resultant |
|---|---:|---:|
| `v` | `3430000` | `27783` |
| `(1-lambda)(1-mu)=1-u+v` | `68600000` | `36006768` |
| `(lambda-mu)^2=u^2-4v` | `68600000` | `781396875` |
| `Disc_x(E)` | `837402343750000` | `593373251953125` |
| extra factor of `A` | `49` | `-2205` |
| `Res_x(S,E)` | `3430000` | `-20003760` |

For reference, after imposing `v=v(u)` and reducing modulo the cubic,

```text
(4,1,1,1):
 Disc_x(E) = -32(7767u^2+12916u+6192)/21875,
 Res_x(S,E)=5(8u^2+9u+8)/28;

(3,2,1,1):
 Disc_x(E) = -16(8817u^2-5332u-4408)/65625,
 Res_x(S,E)=(10u^2+5u-23)/21.                         (18)
```

Thus the poles are distinct and avoid `0,1`, the quadratic `(17)` is
separable, `E` is squarefree, `A` is nonzero, and the simple zero is disjoint
from all double zeros.  The separable cubic followed by the separable
quadratic gives six reduced points.  Direct reduction in each triangular
ordered algebra also makes the pole-collision product, `Disc_x(E)`, `A`,
`Res_x(S,E)`, and the two-equation Jacobian units.  Therefore

```text
I : Gamma^infinity = I,                               (19)
```

where `Gamma` is the full product of those geometric gates.  There are no
discarded collision, confluent-critical, zero-value, or nonreduced branches.

## 6. Exact response reconstruction

At every accessory point put

```text
B=D+A,                  C=-7A,
F=B/D,
G=C E/(2DT),
V=4 S D T^2/C^2.                                        (20)
```

Equations `(3)` and `(5)` give

```text
F'=-A D'/D^2=C E/(DT)=2G,
F=VG^2.                                                 (21)
```

Differentiating the second identity and using the first yields

```text
2VG'+V'G=2.                                             (22)
```

The exact finite fibres are

```text
zero fibre:       (2,2,2,1),
pole fibre:       (4,1,1,1) or (3,2,1,1),              (23)
```

and

```text
F-1=A/D
```

has a zero of order seven at infinity, giving third inertia `(7)`.  The
disjoint simple factor `S` makes `V` genuinely nonsquare in the sense of
THM-2796.  Thus every accessory point constructs a genuine balanced response,
and every normalized response in the two stated passports entered the same
two-equation ideal by `(5)`.

## 7. Marked classes and monodromy

THM-2816 gives exactly

```text
(e-1)! binom(N-e-1,e-1)=2! binom(3,2)=6               (24)
```

ordered Nielsen classes for every ordered heptic maximal-pole passport.
The six algebraic points in each of `(11)` and `(16)` therefore saturate the
topological count.

An independent permutation census fixes

```text
c=(0 1 2 3 4 5 6)
```

as the third inertia and enumerates all `105` three-edge matchings.  Exactly
`35` are planar; their raw pole-passport histogram is the THM-2840 control

```text
(4,1,1,1): 7,        (3,2,1,1): 21,        (2,2,2,1): 7.
```

After quotienting by the rotations generated by `c`, the first two rows have
respectively one and three unmarked classes.  Labelling equal pole cycles
restores six fixed-order charts in each row.

All four new unmarked maps have full symmetric monodromy.  Explicit zero
involutions may be chosen as follows:

```text
(4,1,1,1):
 tau=(0 1)(2 3)(5 6),                 (tau c^3)^5=(1 4);

(3,2,1,1), three unmarked representatives:
 tau_1=(1 2)(3 6)(4 5),               (tau_1 c)^3=(3 5),
 tau_2=(0 1)(3 6)(4 5),               (tau_2 c)^3=(3 5),
 tau_3=(0 5)(1 2)(3 4),               (tau_3 c)^3=(5 6). (25)
```

Conjugating any displayed transposition by the seven-cycle gives a connected
transposition graph on `Z/7`, so each generated group is `S_7`.  Exact group
closure independently returns order `5040` in all four cases.  By contrast,
THM-2840's unique `(2,2,2,1)` Chebyshev map has order-`14` dihedral monodromy.

## 8. Complete abstract heptic consequence

THM-2840 plus the present result complete the normalized abstract
`N=7,e=3,h=4` accessory layer:

```text
pole multiset       fixed-order charts       unmarked maps      monodromy
(4,1,1,1)                    6                     1               S_7
(3,2,1,1)                    6                     3               S_7
(2,2,2,1)                    6                     1               D_14
```

Thus the layer has five unmarked covers: four full-symmetric covers and one
Chebyshev/dihedral cover.  The operation moving pole multiplicity away from
the symmetric `(2,2,2,1)` wall changes monodromy from `D_14` to `S_7`; it
does not create accessory dimension.

## 9. Quotient loss and Keller boundary

The symmetric quotient

```text
(lambda,mu) -> (u=lambda+mu,v=lambda mu)
```

preserves the factorization, response identity, and unmarked simple-pole
pair, but forgets the order of `lambda,mu`.  A chosen root of `(17)`, or the
marked Nielsen chart, restores that label.  For `(3,2,1,1)` this exact
two-to-one loss turns six marked points into three unmarked maps.  For
`(4,1,1,1)`, relabelling all three simple poles further identifies the three
quadratic pairs into one unmarked map.

The map from accessory data to `(F,V,G)` in `(20)` preserves the balanced
square-potential ODE, passport, and nonsplit squareclass.  It forgets the
two-dimensional source prefix, `A_src,B_src,d,s`, and the remaining Faber
flux equations.  To enter the inherited Keller chart one must still form

```text
M=SET
```

and solve the typed conditions

```text
phi_Q=0,             psi_Q in C,             A_src K_Q=lambda M.
```

THM-2827 already excludes these small pole parts from its high-degree
polynomial exact-prefix survivor lanes.  The present classification is
therefore exact accessory mathematics, not a new Keller survivor or degree
descent.

The `S_7` group in Section 7 is the dessin monodromy of these degree-seven
rational maps.  It has no proved transport to the quartic `S_4/C_3`
branchwise-cofactor problem; the covers and preserved predicates differ.

## 10. Exact companion and nonclaims

Run

```bash
python3 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py
python3 -O 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py
cmp 05-knowledge/results/jc_heptic_e3_remaining_accessory_classification_thm3123.out \
    <(python3 04-computation/jc_heptic_e3_remaining_accessory_classification_thm3123.py)
```

The exact companion:

1. rederives `(3)--(16)` by polynomial division and symmetric elimination;
2. checks the exact resultants and every actual quotient-ring saturation
   gate, including the Jacobian;
3. reconstructs `(20)--(22)` coefficientwise;
4. reproduces THM-2840's sextic basis as a positive control;
5. enumerates all `105` matchings, marked/unmarked Nielsen classes, and exact
   monodromy orders; and
6. contains no Python `assert` node, with normal, optimized, and stored
   transcripts byte-identical.

```text
script_sha256 = 12e6d3ae5e0ec253d79eb127987150d59819685db633f4c147fb7058f01accc2
output_sha256 = 42077492ad769c1c09eeb0af6fb83a071898c043ff1a8cc563146b80239c483a
hash_basis   = LF-normalized bytes
```

This theorem does not classify higher `N` or `e`, nonmaximal-pole layers,
nonpolynomial prefixes, or arbitrary Keller charts.  It proves neither
`JC(2)` nor `DC(2)`.
