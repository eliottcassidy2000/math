---
id: THM-3940
title: "I7 rank-two linear-cross-term resolvents still have a unique cubic character"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For q=t^3+t(r+beta t)X-X^2 with beta nonzero, the quadratic resolvent is
  normal, all vertical fibres are integral, and its generic elliptic curve
  has boundary sections Q,-2Q. It has no rational three-torsion for any
  parameters. The non-Cartier Cardano prime and the two-boundary gate force
  Cl(B)[3]=Z/3 and scalar units, so the Cardano cover is the unique smooth-
  locus C3 character up to inversion. On a nonempty parameter open the
  elliptic surface has I7+5I1 and Mordell--Weil rank two, proving that higher
  free rank alone does not create the missing twist. The natural cubic is
  normal and monogenic for every parameter.
source: root / post-THM-3939 rational-three-torsion escape test, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASSES (jc_degree6_one_place and root, 2026-08-24).
  The audits reconstructed the vertical square charts and normality; elliptic
  inverse, boundary sections, `I7+5I1` rank-two open; uniform no-three-torsion
  ladder; Cardano non-Cartier class, two-boundary quotient and scalar units;
  Kummer/descent; and cubic irreducibility, normality and monogenicity. Normal
  and optimized 45-gate runs byte-match the frozen LF transcript; all hashes
  and documentation checks pass.
depends_on:
  - THM-3939-two-boundary-elliptic-resolvent-three-character-rank-one-gate
related:
  - THM-3935-linear-conic-resolvent-class-group-unique-cubic-character
  - THM-3937-linear-conic-fold-three-family-uniform-resolvent-rigidity
script: 04-computation/jc2_i7_rank_two_resolvent_unique_character_thm3940.py
output: 05-knowledge/results/jc2_i7_rank_two_resolvent_unique_character_thm3940.out
script_sha256: 8738427b2d35791c5b4b27a802676f275a07c7bb1638936b2acc6cfa70f05998
output_sha256: 7d8a92fab0549a3690134b0ead2300bdad54975237b156152c1d1c6e6997e7f1
semantic_sha256: b679bee12062d06ac88b036dc222c9bf8f2414a45d8cdf653f7b67ccb1ff3747
hash_basis: raw LF bytes
---

# THM-3940 -- rank two without rational three-torsion is still one line

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Fix
arbitrary `r,beta in k` with `beta!=0`, and put

```text
R=r+beta t,
q=t^3+tRX-X^2,
B=k[t,X,W]/(W^2-(q^2-4X^3)),             S=Spec(B).       (1)
```

Then `B` is normal, every closed fibre over `A1_t` is integral, and

```text
B^*=k*,                              Cl(B)[3]=Z/3.        (2)
```

The generator is the non-Cartier Cardano prime

```text
D+=(X,q+W),                                                (3)
```

and on the regular locus `U=Sreg`,

```text
H^1_et(U,mu_3)=Z/3.                                      (4)
```

The unique nontrivial cyclic cover up to inversion is

```text
z^3=(q+W)/2,             v^3=(q-W)/2,             zv=X.  (5)
```

It recovers the cubic

```text
u^3-3Xu-(t^3+t(r+beta t)X-X^2)=0,                        (6)
```

whose affine ring is normal and monogenic for every `r,beta` in scope.
Consequently no second normal nonmonogenic `S3` cubic completion exists with
this exact quadratic resolvent.

The geometric point of the family is sharp. Its generic elliptic completion
has an `I7` fibre at infinity for every `beta!=0`; on a nonempty Zariski-open
parameter set its remaining fibres are `5I1`, so its Mordell--Weil rank is
two. Nevertheless it has no rational three-torsion section, and the single
boundary relation can create only the Cardano three-class. Thus increasing
free Mordell--Weil rank from THM-3937's one to two does not advance the cubic
counterexample mechanism.

This theorem does not assert that the plane branch has one-place rational
normalization. It closes the proposed higher-free-rank character escape before
that additional branch gate, not nonlinear quadratic coefficients, rational
three-torsion families, higher-genus resolvents, or JC(2).

## 1. Every vertical fibre is integral and the surface is normal

Fix `t=lambda`, and abbreviate

```text
A=lambda^3,                     B_0=lambda(r+beta lambda). (7)
```

The fibre quartic is

```text
(-X^2+B_0X+A)^2-4X^3.                                  (8)
```

If it were the square `(X^2+pX+c)^2`, its `X^3` and constant rows would
give

```text
p=-B_0-2,                         c=+A or -A.             (9)
```

In the plus chart the remaining rows give `A+B_0+1=0` and
`A(B_0+1)=0`; in the minus chart they give `B_0=-1,A=0`.
Either way `A=0,B_0=-1`. But `A=lambda^3=0` implies `lambda=0` and
`B_0=0`, a contradiction. Hence every fibre is a reduced irreducible curve.
In particular, the hypersurface in `(1)` is a domain.

Its partial derivatives, up to nonzero scalar factors, are

```text
W,
q(2X-tR)+6X^2,
q(3t^2+(r+2beta t)X).                                  (10)
```

The `X=0` singular row is only the origin. For `X!=0`, choose

```text
X=a^2,                  q=2 epsilon a^3,
t=u a,                  epsilon in {+1,-1}.              (11)
```

The last two equations of `(10)` become

```text
3u^2+r+2beta u a=0,
epsilon[a(2-beta u^2)-ru]+3=0.                           (12)
```

For `u!=0`, eliminating `a` gives

```text
3beta epsilon u^4-(beta epsilon r+6epsilon)u^2
 +6beta u-2epsilon r=0.                                  (13)
```

Its leading coefficient is nonzero. The row `u=0` is finite as well, so the
singular locus of the two-dimensional hypersurface is zero-dimensional.
It is `R1+S2`, hence normal.

Since every vertical prime is the sole reduced divisor of `t-lambda`, Weil
localization gives

```text
Cl(B)=Pic(C),                                             (14)
```

where `C` is the generic affine fibre.

## 2. The two boundary sections and the I7 deformation

On the generic curve over `K=k(t)`, put

```text
z=(q+W)/(2X^2),
x=t z,
y=Xz(z+1)+(1-tRz)/2.                                    (15)
```

This is birational to

```text
E: y^2=x^3+(t+R^2/4)x^2-(R/2)x+1/4.                     (16)
```

An exact inverse is

```text
z=x/t,
X=t^2(2y-1+Rx)/(2x(x+t)),
W=2X^2z-(t^3+tRX-X^2).                                  (17)
```

The two quartic infinity branches map to

```text
Q=(0,-1/2),
-2Q=(-t,-(1+tR)/2).                                     (18)
```

Indeed the tangent at `Q` has slope `R/2`, so

```text
2Q=(-t,(1+tR)/2).                                       (19)
```

Thus

```text
C=E minus {Q,-2Q},                  (-2Q)-Q=-3Q.          (20)
```

For `(16)`,

```text
b2=4t+R^2,              b4=-R,        b6=1,       b8=t,
c4=(4t+R^2)^2+24R,
c6=-(4t+R^2)^3-36R(4t+R^2)-216,
Delta=-t(4t+R^2)^2-R^3-36Rt-27.                         (21)
```

The finite discriminant has degree five with leading coefficient
`-beta^4`. At infinity, the standard integral scaling has

```text
v_s(c4)=0,                 v_s(c6)=0,       v_s(Delta)=7, (22)
```

so the fibre is `I7`. The generic finite ledger really is `5I1`: the exact
specialization `r=0,beta=1` has

```text
Delta=-t^5-8t^4-17t^3-36t^2-27,
disc_t(Delta)=126759838761,
Res_t(Delta,c4)=12131289.                                (23)
```

Both integers are nonzero, so squarefreeness and multiplicativity hold on a
nonempty parameter open. The global Weierstrass coefficients have the degree
bounds of a rational elliptic surface. On that open, Shioda--Tate gives

```text
rank E(K)=10-2-rank(A6)=2.                               (24)
```

## 3. The three-division polynomial has no rational root

The nonzero three-torsion points of `(16)` have `x`-coordinates among the
roots of

```text
psi_3(x)=3x^4+(4t+R^2)x^3-3Rx^2+3x+t.                   (25)
```

Any root in `K` is integral over `k[t]`, since the leading coefficient is a
unit, and therefore belongs to `k[t]`. Let its degree be `d`. For `d>2`, the
term `3x^4` has unique top degree; for `d=1`, the term `R^2x^3` has unique
top degree; and a constant root is immediately impossible. Hence `d=2`.

Writing its leading coefficient as `a`, the degree-eight row gives
`a=-beta^2/3`. Put

```text
x=-(beta^2/3)t^2+c t+d.                                 (26)
```

The next three coefficient rows are successively

```text
t^7: c=-(2beta r+4)/3,
t^6: d=-r^2/3,
t^5: -beta^5/3=0.                                       (27)
```

The last row contradicts `beta!=0`. Thus

```text
E(K)[3]=0                                                 (28)
```

for every parameter in the theorem, not only on the five-`I1` open.

## 4. The genuine Cardano class forces the unique character

Over `X=0`, the surface has the two primes

```text
D+=(X,q+W),                         D-=(X,q-W),           (29)
```

and

```text
div(X)=D++D-,
div(q+W)=3D+,
div(q-W)=3D-.                                             (30)
```

The ideal of `D+` is reflexive: its quotient is `k[t]`, so the depth lemma
applies. At the origin its two generators have independent linear initials
`X,W`, hence it is not locally principal. Therefore `[D+]` is a nonzero
class of exact order three in `Cl(B)`.

Let `M=E(K)` and let `D=(-2Q)-Q=-3Q`. Equations `(14),(20)` identify

```text
Cl(B)=M/<D>.                                              (31)
```

The difference `D` must have infinite order. If it had finite order, `(28)`
would make that order prime to three; quotienting `M` by such a finite cyclic
subgroup cannot create the nonzero three-torsion class just found in `(30)`.
THM-3939 now applies. It bounds the three-torsion of `(31)` by one line, while
`[D+]` supplies that line. Hence

```text
Cl(B)[3]=Z/3.                                             (32)
```

The same infinite-order boundary difference gives `C^*=K^*`; the integral
vertical ledger then gives `B^*=k^*`. Normal Hartogs and Kummer theory yield
`(4)`. Notice also that `Pic(S)[3]=0`: the only nonzero Weil three-class is
not Cartier. Thus the Cardano cover is quasi-etale over `S`, not etale at the
origin.

## 5. The forced cubic field is normally monogenic

The two nonzero characters in `(4)` are inverse presentations of `(5)` and
define the same cyclic function-field cover. As in THM-3935, the quadratic
deck involution acts by inversion and all `S3` descent lifts are conjugate.
The normal cubic field is therefore unique.

Its natural affine order is the monic ring from `(6)`. The cubic has no
rational root: a polynomial root has total degree at most one; its `X^3`
coefficient forces the `X` coefficient to vanish, after which its `t^2X`
coefficient is `-beta`, a contradiction. The derivative equations, after
putting `X=u^2`, are

```text
2u^2-3u-r t-beta t^2=0,
3t^2+(r+2beta t)u^2=0.                                  (33)
```

Their resultant in `t` has degree six in `u` and leading coefficient
`-8beta^3`, so the singular locus is finite. The cubic hypersurface is a
normal domain and is visibly monogenic over `k[t,X]`.

Therefore the unique normal cubic completion with this exact quadratic
resolvent is the natural monogenic one. The attempted escape has raised the
generic free Mordell--Weil rank, but THM-3939 identifies why that information
is destroyed when only one boundary relation controls three-torsion.

## Reproduction

```bash
python3 04-computation/jc2_i7_rank_two_resolvent_unique_character_thm3940.py
python3 -O 04-computation/jc2_i7_rank_two_resolvent_unique_character_thm3940.py
```

Both runs must byte-match
`05-knowledge/results/jc2_i7_rank_two_resolvent_unique_character_thm3940.out`.
