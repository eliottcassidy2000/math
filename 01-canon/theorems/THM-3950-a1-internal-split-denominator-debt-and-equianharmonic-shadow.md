---
id: THM-3950
title: "A1 internal splits pay denominator debt and an equianharmonic shadow"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the normalization A1 of any nondegenerate one-factor internal-split
  branch, write the two cusp rows as (r^2,2r^3) and (s^2,2s^3). After
  removing gcd(r,s), a nonconstant ratio forces two exact denominator
  factors into the gcd and gives a universal factor parametrization. The two
  undisplayed roots form a connected positive-genus pullback of one fixed
  j=0 elliptic double cover. Equivalently, the assigned-factor ratio is a
  degree-three rational map whose S3 Galois-closure curve is that elliptic
  cover and whose branch values are 0,1,-omega,infinity, exactly the scalar
  collision packet of THM-3947. A degree-one-ratio realization has a reduced
  A1 graph plus an irreducible genus-one residual, a normal quadratic
  surface, and two independent Cardano classes. It is not a Keller map.
source: jc-cohn3709 / normalized-cusp pullback factorization, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-zero-debt-lift, with an independent
  square-class reconstruction by a1_pullback_squareclass, 2026-08-24). The
  audit rederived the UFD cusp extraction and every exceptional linear-form
  case; the residual factorization, nonsquare and genus argument; the
  degree-three S3 closure and isogeny; the explicit irreducible residual and
  squarefree normal surface; and both local A2 generators and their global
  Kummer independence. Normal and optimized runs byte-match the frozen
  51-gate output, all hashes agree, documentation checks pass, and no repair
  was required.
depends_on: []
related:
  - THM-3942-affine-linear-double-torus-factor-split-one-place-obstruction
  - THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse
  - THM-3946-affine-internal-factor-split-two-end-conductor-collision-dichotomy
  - THM-3947-scalar-weighted-repeated-square-split-trichotomy
  - THM-3949-coprime-one-variable-internal-factor-splits-are-reducible-or-multi-ended
  - THM-3941-all-degree-centered-cubic-pole-carrier-routing
script: 04-computation/jc2_a1_internal_split_equianharmonic_shadow_thm3950.py
output: 05-knowledge/results/jc2_a1_internal_split_equianharmonic_shadow_thm3950.out
script_sha256: dd40f6a2fa3e978990811aaab0ff9d690688a01940e7230150ad979874644732
output_sha256: e3936926cba7ddb0dc3e5fa9748d5de9074f467a599eb23ff4d62c19449d3966
semantic_sha256: 7b027947c393486b42e7b8fb30c0cfdbf9edb2af03c26ebd2121816e5f542402
hash_basis: raw LF bytes
---

# THM-3950 -- every nonconstant A1 split carries the same j=0 shadow

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Fix

```text
omega^2+omega+1=0,       delta=omega-omega^2,       delta^2=-3.       (1)
```

The theorem first classifies the pullback of any nondegenerate internally
split double-torus row to a branch normalized by `A1`. It then realizes the
degree-one-ratio packet as an explicit plane discriminant. That example
simultaneously attains a reduced one-place component, a normal quadratic
surface, and two independent Cardano characters. Its unavoidable companion
is elliptic, so this is a sharp structural survivor rather than a planar
Jacobian counterexample.

## 1. Universal A1 pullback classification

Let `Gamma` be an irreducible affine curve with normalization
`Gamma^nu=A1_t`. Suppose polynomial functions on `Gamma` obey

```text
q0^2=4p0^3,                         q1^2=4p1^3,
p1-p0=A B,
q1-q0=2 lambda A(p1-omega p0),
q1+q0=2 lambda^(-1)B(p1-omega^2 p0),                  (2)
```

where `lambda in k*`, and assume the pullbacks of `p0`, `p1`, and `p1-p0`
are nonzero in `k[t]`. Put

```text
a=lambda A|Gamma,                 b=lambda^(-1)B|Gamma.               (3)
```

The UFD `k[t]` turns the two cusp equations into

```text
p0=r^2,       q0=2r^3,       p1=s^2,       q1=2s^3                  (4)
```

for nonzero `r,s in k[t]`; signs are absorbed into `r,s`. Write

```text
r=gR,              s=gS,              gcd(R,S)=1.                    (5)
```

If `S/R` is constant, then `R,S` are constants and `(4)` is a scalar
repeated-square packet on the branch:

```text
p1-p0=(S^2-R^2)g^2.                                      (6)
```

This is only the normalized pullback shadow of THM-3947's collision grammar;
no global repeated-square factorization is inferred from `(6)`.

Suppose instead that `S/R` is nonconstant. Define

```text
U=S+omega^2 R,                         V=S-omega R.       (7)
```

Then, for some nonzero `c in k[t]`,

```text
g=cUV,
a=c(S-R)V^2,
b=c(S+R)U^2.                                           (8)
```

Conversely, `(7)-(8)` produce every nonconstant-ratio solution of `(2)` on
`A1`: with `r=cRUV` and `s=cSUV`, all four identities `(2)` hold.

For necessity, substitute `(4)-(5)` into the last two rows of `(2)` and
factor the sixth-power difference. Cancellation gives

```text
a U = g(S-R)V,                  b V = g(S+R)U.            (9)
```

Distinct linear forms in the coprime pair `(R,S)` are coprime. Hence

```text
gcd(U,(S-R)V)=1,                gcd(V,(S+R)U)=1.          (10)
```

Thus `U|g`, `V|g`, and `gcd(U,V)=1`, proving `(8)`. This is the exact
**denominator debt**: a nonconstant root ratio is possible only because both
missing linear denominators are paid by the common square root `g`.

## 2. The two remaining roots are an equianharmonic shadow

Attach to `(8)` the cubic in an auxiliary cusp root `h`

```text
D=(1-omega^2)b-(1-omega)a,
E=(b-a)ab,
C_t(h)=E+D h^2-2h^3.                                  (11)
```

Its known root is `h=r`, and

```text
C_t(h)=(h-r)Q_t(h),
Q_t(h)=-2h^2+(D-2r)h+r(D-2r).                          (12)
```

The residual discriminant is

```text
disc_h(Q_t)=(D-2r)(D+6r)
 =c^2/3 (R-S)(R+S)(R+delta S)(3R+delta S)^3,            (13)
```

because

```text
D=3r+delta cS^3,
D-2r=-c(R-S)(R+S)(R+delta S),
D+6r=-(c/3)(3R+delta S)^3.                              (14)
```

Modulo squares, `(13)` is

```text
F(R,S)=(R-S)(R+S)(R+delta S)(3R+delta S).               (15)
```

This is the pullback under `x=R/S` of

```text
mathcal E: y^2=f(x),
f(x)=(x^2-1)(x+delta)(3x+delta).                        (16)
```

The four branch points are distinct. In the standard binary-quartic
convention, `I=0` and `J=1728!=0`, so `mathcal E` is a smooth genus-one
curve with `j=0`.

Let `C_res` be the smooth projective normalization of the residual double
cover. The map

```text
(t,z) |-> (x=R(t)/S(t), y=z/S(t)^2)                    (17)
```

extends to a nonconstant morphism `C_res -> mathcal E`. If `(15)` were a
square in `k(t)`, a component `P1_t` would map nonconstantly to `mathcal E`,
contradicting Riemann--Hurwitz. Hence `Q_t` is irreducible over `k(t)`, the
residual cover is connected, and another Riemann--Hurwitz application gives

```text
genus(C_res)>=1.                                        (18)
```

Equality means `(17)` is unramified, hence an isogeny. The minimal-genus
ratios are therefore the Lattes/isogeny pullbacks of this fixed four-color
cover; generic ratios have larger residual genus.

## 3. The degree-three map behind all four colors

The assigned-factor ratio is

```text
a/b=phi(R/S),
phi(x)=((1-x)(1-omega x)^2)/((1+x)(1+omega^2 x)^2).     (19)
```

The numerator of `phi(x)-z` has discriminant

```text
-6z(z-1)[(delta+3)z+(delta-3)],                         (20)
```

so this degree-three map has branch values

```text
{0,1,-omega,infinity}.                                  (21)
```

They are exactly the compactified scalar collision values in THM-3947. The
ramification signature is `C2^4`, with no index-three point. Thus `phi` is
not Mobius-conjugate to the polynomial cubic projections routed by THM-3941:
it is the complementary rational, split-infinity cubic carrier. Its `C3`
appears in the `S3` Galois closure and three-isogeny, not in the inertia of
`phi` itself. The
off-diagonal equation `phi(y)=phi(x)`, after deleting `y=x`, is

```text
x^2+y^2+4xy-delta(x^2y+xy^2)+delta(x+y)=0,              (22)
```

whose discriminant in `y` is `-f(x)`. Thus `(16)` is the `S3`
Galois-closure curve of `(19)`, not an accidental genus-one byproduct. The
two standard `j=0` models are linked by the exact three-isogeny square

```text
z(z-1)(z+omega)
 =(3+delta)/6
  [x(x-omega^2)/((x+1)^2(x+omega)^3)]^2 f(x),
z=phi(x).                                               (23)
```

The scalar wall, denominator payment, residual elliptic sheet, and recurrent
`C3` character are four views of this one degree-three cover.

## 4. A reduced one-place, normal, rank-two survivor

Take

```text
R=1,                    S=Y,                    c=1,
g=(Y+omega^2)(Y-omega)=Y^2-delta Y-1,
A=(Y-1)(Y-omega)^2,
B=(Y+1)(Y+omega^2)^2.                                  (24)
```

In `k[P,Y]`, put

```text
p0=P,                         p1=P+AB,
q1-q0=2A(p1-omega p0),
q1+q0=2B(p1-omega^2 p0),
H=q0^2-4p0^3=q1^2-4p1^3.                               (25)
```

The graph

```text
Gamma: P=g(Y)^2                                         (26)
```

is a reduced component of `H=0`, isomorphic to `A1_Y`, with one
normalization place at infinity. Along it,

```text
(p0,q0)=(g^2,2g^3),
(p1,q1)=((Yg)^2,2(Yg)^3).                              (27)
```

With

```text
D=(1-omega^2)B-(1-omega)A,
N=4P^2-(D^2-4g^2)P+g^2(D-2g)^2,                        (28)
```

one has

```text
H=(g^2-P)N,
disc_P(N)=(D-2g)^3(D+6g).                              (29)
```

The square class in `(29)` is

```text
(1-Y)(1+Y)(1+delta Y)(3+delta Y),                       (30)
```

with four simple roots. Hence `N` is irreducible and its projective
normalization has genus one. The graph factor has multiplicity one because

```text
H_P(g^2,Y)=4delta g^3Y^3!=0.                            (31)
```

Thus the degree-four graph and degree-ten residual are distinct irreducible
factors of the degree-fourteen full discriminant.

The quadratic surface

```text
S=Spec k[P,Y,W]/(W^2-H)                                (32)
```

is normal. Equations `(29)-(31)` make `H` squarefree, so `(32)` is a domain;
its singular locus lies over the finite singular set of the reduced plane
curve. A hypersurface is `S2`, and this codimension-two singular locus gives
`R1`.

Define height-one primes

```text
D0=(P,q0+W),                    D1=(p1,q1+W).            (33)
```

The norm identities give

```text
div(q0+W)=3D0,                    div(q1+W)=3D1.         (34)
```

There are two exclusive local witnesses. First,

```text
B-A=-delta Y-1                                         (35)
```

has the simple root `Y=delta/3`; there `P=q0=W=0` while `p1!=0`. The
completed surface is `uv+4P^3=0`, an `A2` point, and `D0` is its local class
generator while `D1` is absent. Second,

```text
A+omega B=Y^2[(1+delta)Y+3-delta]/2.                   (36)
```

At the simple nonzero root of the bracket, put `P=-AB`; then
`p1=q1=W=0` while `p0!=0`. This is another `A2` point, with `D1` the local
generator and `D0` absent. Restriction to the two local class groups proves

```text
< [D0],[D1] > ~= (Z/3)^2  inside Cl(S).                 (37)
```

Since `S` is normal, `Cl(S)=Pic(Sreg)`. The radicands in `(34)` therefore
define independent etale `mu_3` classes on `Sreg`. Unlike THM-3944, neither
class is merely supported on a nonnormal conductor.

## 5. Meaning, scope, and reproduction

The example crosses three gates that had previously failed one at a time:

```text
reduced A1 branch  +  normal quadratic surface
                   +  two independently extendable Cardano classes.       (38)
```

It pays for the one-place sheet with a positive-genus `j=0` companion
controlled by the same degree-three map as the scalar collision wall. This
is the correct successor to merely trying more coefficients in THM-3946 or
THM-3949.

No affine-plane normalization of either natural cubic, no polynomial Keller
pair, and no counterexample to `JC(2)` is constructed here. Those require a
separate same-field cubic-atlas and boundary audit. Distributions across
multiple `Li`, gcd overlap between distinct cube factors, and non-`A1`
primary branches remain open.

Run

```bash
python3 04-computation/jc2_a1_internal_split_equianharmonic_shadow_thm3950.py
python3 -O 04-computation/jc2_a1_internal_split_equianharmonic_shadow_thm3950.py
```

for the exact 51-gate companion. It freezes `(19)-(23)`, the degree
`14=4+10` factorization, the reduced-factor derivative, and both exclusive
`A2` witness resultants.
