# Sextic normal-strip composite frontier scout

**Status: FINITE-EXACT + SCOUT.  Not a sextic-strip theorem.  `JC(2)` remains
open.**

This note records the exact first obstruction after
[`THM-3871`](../../01-canon/theorems/THM-3871-quintic-normal-strip-keller-pairs-are-automorphisms.md).
The computation is reproducible, symbolic over `QQ`, optimization-safe, and
uses no sampled coefficients or asserted roots.

## 1. Inheritance and exact top split

The closest proved mechanism is THM-3871: normalize the top target direction,
factor the highest Wronskian in the UFD `k[s]`, depress the lower target, and
retain every conserved penultimate bucket.  Its canonical hostile is the
degree-five arms-only cancellation family, which is killed only by conserved
buckets.  The repaired near miss for degree six is sharper: in `(6,4)`, even
the arms **and** both conserved leading equations retain a nonreduced cusp.
The relevant method card is “divide exceptional multiplicity before judging a
wall.”

Write

```text
A=sum_(i=0)^6 a_i(s) z^i,       C=sum_(j=0)^6 c_j(s) z^j.
```

There are exactly twelve possible Jacobian buckets:

```text
E_m=sum_(i+j=m+1) (i a_i c_j'-j a_i' c_j),    0<=m<=11.    (1)
```

The top row is

```text
E_11=6(a_6 c_6'-a_6'c_6).                              (2)
```

Thus a constant target `SL_2(k)` change makes `c_6=0`, `a_6=w!=0`.
If `j=deg_z C` and `q=c_j`, the highest surviving row is

```text
6wq'-jw'q=0.                                             (3)
```

For `g=gcd(6,j)`, prime valuations give

```text
w=R h^(6/g),          q=Q h^(j/g),       R,Q in k*.       (4)
```

Consequently the exact row table is

| row | reduced top packet | consequence |
|---|---:|---|
| `(6,0)` | -- | empty by the highest nonconstant bucket |
| `(6,1)` | `(6,1)` | shear by a scalar multiple of `C^6` |
| `(6,2)` | `(3,1)` | shear by a scalar multiple of `C^3` |
| `(6,3)` | `(2,1)` | shear by a scalar multiple of `C^2` |
| `(6,4)` | `(3,2)` | new cusp/Kummer row |
| `(6,5)` | `(6,5)` | new coprime Kummer row |

The first three shears lower `deg_z A` to at most five, so THM-3871 makes
those pairs automorphisms.  Only `(6,4)` and `(6,5)` remain genuinely new.

## 2. The `(6,4)` packet and its nonreduced cusp

Equation (4) gives `w=Rh^3`, `q=Qh^2`.  After adjoining `eta` with
`eta^2=h`, put `y=eta z`; then the rational bracket is `lambda/eta`.
Depression gives

```text
A_*=R x^6+D x^5+U x^4+M x^3+L x^2+P x+a,
C_*=Q x^4+B x^2+N x+b,                    D'=0.           (5)
```

The six high rows integrate uniquely, the next two rows are derivatives of
two explicit conserved polynomials, and the constant row is

```text
P b'-a'N=lambda/eta.                                     (6)
```

All coefficients and both conserved polynomials are verified exactly in the
companion.  No division by `B`, `N`, `b`, a conserved value, or a discriminant
is used.

On the principal balanced pole face, normalize

```text
X=B/(Qg^2),          Y=N/(Qg^3),          Z=b/(Qg^4).
```

Regularity of both arms and conservation give

```text
1+X+Y+Z=0,                                             (7)
Y(4Z-X^2)=0,                                           (8)
(X^2-4Z)^2-8XY^2=0,                                    (9)
6Y^2-(X+2)^3=0.                                       (10)
```

Set `u=X+2`, `v=Y`, and eliminate `Z=1-u-v`.  A reduced Groebner basis of
the **scheme ideal** is equivalent to

```text
6v^2-u^3,            u^2(2u+3v),            u^4.        (11)
```

Hence the reduced support is the unique point

```text
(X,Y,Z)=(-2,0,1),                                      (12)
```

but the local algebra has basis

```text
1,u,u^2,u^3,v,uv                                      (13)
```

and length six.  This is a genuine nonreduced survivor, not a missed scalar
resultant.

The normalized leading constant bracket is

```text
F_0=Y(X^3-4XZ-6Y^2).                                   (14)
```

It vanishes in the local algebra.  More explicitly, after (12),

```text
F_0=(-v+2u/3-4/3)(6v^2-u^3)
    -(2/3)u^2(2u+3v)+(2/3)u^4.                          (15)
```

Thus the leading face cannot pay a nonzero Keller constant.  It does **not**
yet follow that every strict transform is empty: one must divide the proved
exceptional multiplicity and compute the first surviving normal jet.

There is already one exact strict-transform step.  In the primitive
simple-pole chart put `t=g^-1`, allow a regular `C` arm to begin at `t^4`, and
write

```text
u=u_1t+u_2t^2+...,             v=v_1t+v_2t^2+....         (16)
```

Successive arm/conservation coefficients are

```text
t^1:  A_arm=3D/8,       K_2=-5D/4,
t^2:  A_arm=3v_1^2/8,
t^3:  K_2=3F/2,         A_arm=-u_1^3/16  after F=0,
t^4:  K_2=-3v_2^2,
t^5:  A_arm=I,
t^6:  K_2=-3v_3^2.                                      (17)
```

Hence this chart forces

```text
D=F=I=0,                 u=O(t^2),        v=O(t^4).       (18)
```

This is the first strict-transform dividend: the three odd sextic constants
die, and the remaining ratios enter at the cusp weights.  It is still not a
closure.  The next chart contains `u_2,v_4`, the regular target-arm coefficient
at `t^4`, the even constants `E,G,A0`, and both conserved constants.  Moreover,
for a pole of order greater than one there can be intermediate uniformizer
jets not visible in the primitive `t=g^-1` chart.

The support point is the exact cusp composition

```text
T=s^2z^2+2z,                  A=T^3, C=T^2, J(A,C)=0.    (19)
```

It is the canonical hostile against deleting the constant bucket.  Arms
alone are much weaker: `(X,Y,Z)=(4,6,-11)` cancels both arms but leaves the
two conserved residuals `-360` and `2448`.

### Quadratic-extension loss ledger

The map to (5) preserves the bracket, transverse degrees, arms, and all
coefficient buckets.  It forgets descent to `k(s)` and original polynomial
integrality.  The required sidecar is the involution

```text
(eta,y)->(-eta,-y)
```

together with parity and valuation of every original coefficient.  Any
future `(6,4)` closure must audit this sidecar after the strict transform.

## 3. The `(6,5)` packet and the decisive extra bucket

Here (4) is `w=Rh^6`, `q=Qh^5`, so `y=hz` needs no algebraic extension.
After depressing the quintic,

```text
A_*=R x^6+D x^5+U x^4+M x^3+L x^2+P x+a,
C_*=Q x^5+B x^3+N x^2+Vx+b,                 D'=0.        (20)
```

The six high rows integrate uniquely.  The next two rows are exact
derivatives of conserved polynomials.  Unlike degree five, the remaining
`x^1` row is a genuinely nonexact one-form; only after it comes the constant
one-form.

On the principal balanced pole face use

```text
X=B/(Qg^2), Y=N/(Qg^3), Z=V/(Qg^4), W=b/(Qg^5).
```

The target arm, two conserved rows, the penultimate row, source arm, and
constant leading bracket are respectively

```text
C_0=1+X+Y+Z+W,                                           (21)
F_3=-2X^2Y+5XW+5YZ,                                      (22)
F_2=3X^4-20X^2Z-20XY^2+50YW+25Z^2,                       (23)
F_1=25WX^2+225WZ+8X^3Y-65XYZ-30Y^3,                      (24)
A_0=-4X^3+15X^2+30XY+30XZ+15Y^2-25,                      (25)
F_0=125W^2+25WXY+4X^3Z-30XZ^2-15Y^2Z.                    (26)
```

The ideal `(C_0,F_3,F_2,A_0)` is a reduced zero-dimensional scheme of length
fourteen.  Its lexicographic basis is triangular, with a squarefree degree-14
terminal eliminant in `X`.  Thus deleting the `x^1` row creates fourteen
honest algebraic apparent branches, not a positive-dimensional artifact.

Exact Groebner reduction gives

```text
(C_0,F_3,F_2,A_0,F_1)=(1),
(C_0,F_3,F_2,A_0,F_0)=(1).                               (27)
```

Therefore `F_1` kills all fourteen apparent branches, while `F_0` is nonzero
at every one of them.  This closes the principal balanced `(6,5)` face and
shows precisely why the extra nonexact bucket is indispensable.  The
primitive terminal eliminant has coefficient-list SHA-256
`0e3c6807f6d25f289fd3e1116a59d4aa01749e838efd6e6a901a051df7eab4ef`.

An arms-only rational hostile is

```text
(X,Y,Z,W)=(1,0,7/15,-37/15),                             (28)
```

with residuals

```text
(F_3,F_2,F_1,F_0)=(-37/3,-8/9,-962/3,6803/9).             (29)
```

## 4. What is proved, and what is not

**FINITE-EXACT:**

- all twelve universal buckets and the gcd/top-factor table;
- target-shear closure of `(6,1)`, `(6,2)`, `(6,3)` via THM-3871;
- complete depressed rational packets for `(6,4)` and `(6,5)`;
- the length-six nonreduced principal balanced `(6,4)` cusp scheme and the
  leading-bracket ideal certificate (15);
- the reduced fourteen-point principal balanced `(6,5)` apparent scheme and
  its exact annihilation by `F_1`.

**OPEN:**

- strict-transform jets above (18), higher-order pole charts, and quadratic
  descent in `(6,4)`;
- the `(6,4)` channels where the depression shift is regular; only a subset
  of `B,N,b` has a pole; their pole orders are not in the principal `2:3:4`
  ratio; a leading coefficient vanishes; or the Kummer square root is
  ramified, unramified, or already rational;
- the `(6,4)` unit-scale/infinity channel, all-constant channel, and every
  cancellation stratum of the two conserved values;
- in `(6,5)`, the regular-depression channel; every proper subset of
  `B,N,V,b` having poles; every nonprincipal pole-order ray distinct from
  `2:3:4:5`; all ties among leaders; and every zero polynomial among
  `B,N,V,b`;
- the `(6,5)` unit-scale/infinity and all-constant channels, plus finite-zero
  strata of any coefficient combination created during elimination;
- a full pole/unit/finite-zero/identically-zero carrier split analogous to
  the quintic `mathcal W` split.  No canonical sextic carrier of that kind has
  yet been isolated; the `W` in (21)--(26) is merely normalized `b`, not such
  a carrier;
- return from depressed rational coefficients to simultaneous regularity of
  both original arms and then to original polynomial integrality in every
  preceding channel;
- the full degree-six normal-strip theorem;
- arbitrary polynomial Keller pairs and `JC(2)`.

The next generated tasks, in order, are:

1. continue (18) with variables `u_2,v_4,c_4,E,G,A0` and the two conserved
   constants, divide the `u^4` exceptional multiplicity, and compute the first
   nonzero constant bracket;
2. retain the `eta`-parity sidecar through that strict transform;
3. exhaust the nonbalanced valuation fans and all zero/unit/constant-scale
   channels in `(6,5)` before attempting theorem status;
4. independently replay the two Groebner computations with a second CAS or
   coefficient-elimination path.

## 5. Reproduction

```bash
python3 -B 04-computation/jc2_sextic_normal_strip_composite_frontier_scout_20260823.py
python3 -B -O 04-computation/jc2_sextic_normal_strip_composite_frontier_scout_20260823.py
```

Both executions must byte-match
`05-knowledge/results/jc2_sextic_normal_strip_composite_frontier_scout_20260823.out`.

```text
semantic_sha256 = 103b0b839db044fcec56c18e87516dda4aea2da93c4380e7dd5777c72be2ef4f
script_sha256   = 0f23f081efe9359699cfbbe34f1169832ad1d86b042a66710856a1c9131b6e02
output_sha256   = 0ee27143276ee5f0883dc956545d12704a7e6d00c5899e94a559dbaa24a0099a
checks          = 76
```
