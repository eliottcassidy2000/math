---
source: codex-2026-07-25-jc2-central-wall-s0
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; SUPERSEDED AS
  A CLOSURE STATEMENT BY THM-2345. In the genuine
  nonsplit polynomial exact-square-prefix degree-eighteen branch of
  THM-2262/2297, every point with D/B^2=25/126, B!=0, and S=0 is empty.
  An exhaustive one-parameter spectral normalization has genus four
  generically, genus one at one ordinary-triple value, and genus three at
  one nodal value and seven algebraic nodal values. This closes one central
  repeated-branch face. The later audited
  `jc2-degree18-central-twall-closure-opus-20260725` packet combines it
  with R=0 and T=0 to close the full central ratio, not JC(2).
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2313-degree-eighteen-bd-linear-ratio-closure
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-central-wall-r0-closure-opus-20260725
  - jc2-degree18-central-twall-closure-opus-20260725
---

# The degree-eighteen central `S=0` face is empty

> **CANONICAL SUPERSESSION.** THM-2345 closes the entire common-root
> wall `126D=25B^2`.  This file is retained for its exhaustive
> one-parameter genus and singularity atlas.

## 1. The invariant and exhaustive weighted normalization

On the central wall

```text
D=25B^2/126,
```

the named weight-ten invariant is

```text
S
 =2888B^5+108864B^2C^2+571536BCW+750141W^2

 =1701(21W+8BC)^2+2888B^5.                         (1)
```

Assume `B!=0`. Over the algebraically closed constant field choose a
weighted scale `a` with

```text
B=-7a^2/54.
```

Changing `a` to `-a` absorbs the two square-root signs in `S=0`.
Consequently the entire face, not merely a sample line, has the exact
parameterization

```text
B=-7a^2/54,

D=175a^4/52488,

C=xi a^3/26244,

W=(2xi+399)a^5/1062882.                            (2)
```

Weighted scaling gives an isomorphism of spectral curves, so it is legal
to compute at `a=1`. This does not set a nonconstant Keller coordinate to
one; `a` is a constant weighted-projective scale.

The point `C=W=0` is not on this face:

```text
S=2888B^5!=0.                                      (3)
```

Thus this closure is disjoint from reserved THM-2313's `B--D` line.

## 2. The normalized trigonal curve

Make the affine spectral change

```text
U=-(1701/4)[u+5(4-y^2)/243].
```

Direct substitution into THM-2262's translated cubic gives, up to the
nonzero scalar `64/189`,

```text
F_xi(U,y)=U^3+3pU-7y q_xi,                          (4)

p=245y^2(y^2-1),

q_xi
 =539y^5-1470y^3+7xi y^2+7(133-xi).
```

Put

```text
tau=xi/7.
```

The branch discriminant factors as

```text
-(1/27)Disc_U(F_xi)
 =117649 y^2(y-1)^2 R_tau(y),                       (5)
```

where

```text
R_tau
 =621y^8+1242y^7-297y^6+(22tau-1836)y^5
  +(44tau-975)y^4+(304-16tau)y^3
  +(tau^2-76tau+1083)y^2
  +2(tau-19)^2y+(tau-19)^2.                        (6)
```

The boundary values are

```text
R_tau(0)=(tau-19)^2,

R_tau(1)=(2tau-35)^2.                               (7)
```

The exact residual discriminant is

```text
Disc_y(R_tau)
 =-49682096496000000000000
   (tau-19)^6 P_7(tau),                             (8)
```

with

```text
P_7
 =1242tau^7+65205tau^6+3949300tau^5
  -287555625tau^4+28949231250tau^3
  -1626234553125tau^2+31352481000000tau
  -183687836484375.                                 (9)
```

The polynomial `P_7` is squarefree and irreducible. A compact exact
certificate is its monic reduction modulo `37`:

```text
z^7+34z^6+5z^5+20z^4+6z^3+6z^2+20z+10.            (10)
```

The Rabin test gives

```text
gcd(P_7,z^37-z)=1,

z^(37^7)=z mod P_7.                                 (11)
```

Also

```text
P_7(19)!=0,                  P_7(35/2)!=0,          (12)
```

so the two rational exceptional parameters below are disjoint from the
seven algebraic ones.

## 3. Generic genus four

For generic `xi`, the fibres `y=0` and `y=1` are smooth totally ramified
fibres of the cubic projection and contribute two ramification units each.
The eight roots of `R_tau` are simple and contribute one each. Thus the
total ramification is

```text
2+2+8=12.
```

At infinity put `V=U/y^2`. The limiting cubic is

```text
V^3+735V-3773,

Disc=-1972620783!=0.                                (13)
```

There are exactly three smooth unramified points over infinity.
Riemann--Hurwitz therefore gives

```text
2g-2=3(-2)+12=6,

g=4.                                                (14)
```

## 4. The ordinary-triple value `xi=133`

Here `tau=19` and

```text
R_tau=y^4 R_4,

R_4=621y^4+1242y^3-297y^2-1418y-139.               (15)
```

The polynomial `R_4` is squarefree, with

```text
R_4(0)=-139,                  R_4(1)=9.
```

At the origin the tangent cubic is

```text
a^3-735a-6517,

Disc=441536697!=0.                                  (16)
```

Thus the origin is an ordinary triple point. Its three normalization
branches are unramified for the `y` projection. The totally ramified fibre
at `y=1` contributes two and the four simple roots of `R_4` contribute
four, so

```text
ramification=2+4=6,

g=1.                                                (17)
```

## 5. The rational nodal value `xi=245/2`

Here `tau=35/2` and

```text
R_tau=(1/4)(y-1)R_7,                               (18)

R_7
 =2484y^7+7452y^6+6264y^5+460y^4
  -360y^3-264y^2-27y-9.
```

It is squarefree, with

```text
R_7(0)=-9,                    R_7(1)=16000.
```

At `(U,y)=(0,1)` the tangent cone is

```text
(735/2)Y(4U-35Y),                                  (19)
```

so the point is an ordinary node. Its vertical branch has

```text
Y=-U^2/1470+O(U^3)
```

and contributes exactly one ramification unit; the other node branch is
unramified. The central triple fibre contributes two and the seven simple
roots contribute seven:

```text
ramification=2+1+7=10,

g=3.                                                (20)
```

## 6. The seven algebraic nodal values

At every root of `P_7`, the valuation of (8) is one. Hence `R_tau` has
exactly one double root and six simple roots. The double root is away from
`0,1`.

The common critical-section identity places the double branch root on an
affine singularity of (4). It is nontriple, so `F_UU!=0`; the local
critical-section normal form makes it an ordinary node. Both node branches
are unramified. The two smooth triple fibres contribute four and the six
simple residual roots contribute six:

```text
ramification=4+6=10,

g=3.                                                (21)
```

## 7. Connectedness and the Keller contradiction

The cubic (4) is irreducible for every `xi`. Any polynomial root over
`C[y]` must have the form

```text
U=a y(y-1).
```

Coefficient comparison first forces `xi=133` and then

```text
a(a^2+245)=0;
```

each possibility contradicts the remaining leading-coefficient equation.
Gauss's lemma therefore gives irreducibility over `C(y)`, so every
normalization used above is connected.

Every case has positive genus: respectively `4`, `1`, or `3`. A rational
Keller trajectory from `P^1` to that normalization is constant. If the
constant `y` is nonzero, THM-2262 equation (11) makes `Z=T^2` constant,
then `T` and the nonzero deck-odd `q` are constant, contradicting the
genuine deck. If `y=0` identically, THM-2262's retained third flux, Keller
one-form, and whole-Faber polynomial sidecar close the exceptional wall.

Therefore

```text
D/B^2=25/126,               B!=0,               S=0

has no survivor.                                      (22)
```

## 8. Combined central-wall effect

The independently audited companion result
`jc2-degree18-central-wall-r0-closure-opus-20260725.md` closes

```text
R=20BC+21W=0,                  S!=0.
```

Because `S=0` is already empty by (22), the **whole `R=0` central wall is
empty**. In the current factorization, every remaining central-wall
survivor must lie on

```text
R!=0,                     S!=0,                     T=0,        (23)
```

where `T` denotes the primitive weight-thirty repeated-branch factor.
Equation (23), the other discriminant components, split deck, even-leading
descent, and other Newton edges remain open. Neither this theorem nor its
combination proves JC(2) or DC(2).

## 9. Information ledger

```text
source:
  the central weight-ten factor S and THM-2262's trigonal spectrum;

map:
  parameterize the weighted S=0 face exactly, normalize every singular
  trigonal curve, and retain the first-flux/deck sidecars;

preserved:
  the full S=0 orbit, branch multiplicities, normalization component,
  infinity, first flux, genuine deck, and exceptional y=0 closure;

destroyed by raw discriminant:
  singularity type and normalization genus;

restoring sidecar:
  ordinary-triple/node local models plus Riemann--Hurwitz;

next decisive target:
  reconstruct and normalize the primitive weight-thirty T=0 factor while
  retaining Z=T_flux^2 and the Keller one-form.                       (24)
```

## 10. Exact reproduction

The standard-library companion and stored transcript are

```text
04-computation/jc2_degree18_central_s0_closure_probe.py

05-knowledge/results/jc2_degree18_central_s0_closure_probe.out.
```

Their LF byte hashes are

```text
script:
  cbab90cc024eec11e4685557f367eac57cc26cb6a87504d310d0496eb951d982;

output:
  a7687d39a9c986dd674b890d0c2bf1a3fc4d794d4b121ac4cc01cec9b2ba8965.
```

Normal and optimized runs are byte-identical to the stored output. The
companion reconstructs the `15 x 15` Sylvester determinant over
`Q[tau]`, verifies (1), (6)--(12), both rational specializations,
squarefreeness of `R_4,R_7`, the exact mod-37 Rabin certificate, infinity,
and the complete genus ledger. The local singularity and Keller/deck
arguments remain the mathematical proof above rather than computer
assumptions.
