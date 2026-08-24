---
id: THM-3921
title: "Quintic genus-collapse decic degeneration, normal order, and persistent class group"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The quintic
  genus-zero locus of THM-3917 is the exact degeneration on which the
  rational one-place decic's origin fibre changes from eight smooth branches
  to four smooth plus two cuspidal branches, still with delta 28. The six
  collapsed-ramification addresses are exactly these six origin branches.
  The other two A2 and two A5 singularities persist and complete the genus
  ledger. Despite the geometric degeneration, every algebraic invoice of
  THM-3915 persists: an explicit normal finite-free rational cubic order is
  globally nonmonogenic with power-basis index A^5, and the quadratic
  resolvent still has scalar units and class group Z^3 direct-sum Z/3. The
  candidate is excluded by THM-3917's six-branch boundary obstruction.
  The reserved THM-3924 candidate supplies a prospective independent
  primitive-class obstruction, but is not used here until its own audit;
  JC(2) remains open.
source: jc_degree6_one_place / post-THM-3917 decic degeneration lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (audit_thm3921/root, 2026-08-23), after
  repairs 893fab9554 and 642ee9738e. Normal, optimized, and frozen runs
  byte-match in 69 active gates and both raw hashes agree. The decic
  normalization, singularity/genus ledger, six-address correspondence,
  integral basis and A^5 index, normality, nonmonogenicity, node-blowup
  boundary span, class group, units, and actual Kummer generator were all
  independently reconstructed. The audit corrected “generate the unit
  ideal” to the exact gcd-one/proper-ideal statement and scoped Kummer
  unramifiedness to codimension-one valuations centred on the affine
  quadratic resolvent. No mathematical gap remains.
depends_on:
  - THM-3917-quintic-parameter-rational-collapsed-cubic
related:
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
  - THM-3914-decic-boundary-three-class-degree-one-isotropic-divisor
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
  - THM-3924-decic-cubic-index-five-ramification-class-obstruction
script: 04-computation/jc2_quintic_decic_degeneration_order_thm3921.py
output: 05-knowledge/results/jc2_quintic_decic_degeneration_order_thm3921.out
script_sha256: 5f0db36702cfdd09296000f39dc5210bd6d8a2e6d385d0c1e07a4873198cf245
output_sha256: 63e0fd607c0db03aa73f579c9b6d53cdb149e37e48997569847a459bf53cf076
semantic_sha256: d273b4c2ebba004daf98098437a91b8c7c358aa728bab8cad90d1353411adf00
hash_basis: raw LF bytes
---

# THM-3921 -- genus collapses while every algebraic invoice persists

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Retain the notation
of THM-3917:

```text
K(b)=2304b^5+10176b^4+4064b^3+996b^2+84b+5=0,

p=(u^2-1)(u^2+b)^2,
h=the degree-nine polynomial part of p^(3/2),
r=p^3-h^2,                    F=z^3-3pz+2h.                (1)
```

The rational branch curve attached to the chart

```text
A=F/4,                         C=uF/4                     (2)
```

is an irreducible one-place decic. Its complete singularity packet is

```text
origin: four smooth branches + two A2-cuspidal branches,
        all six tangent directions distinct, delta=28;

away from origin: two A2 cusps + two A5 contacts.           (3)
```

The six points in which THM-3917's collapsed curve meets its rational
ramification curve are exactly the six normalization addresses above the
origin in `(3)`: the four simple addresses give the four smooth branches,
and the two double addresses give the two cuspidal branches.

This geometric trade leaves the algebraic completion intact. The degree-
three field has an explicit finite-free normal order `B` whose fraction
field is `k(u,z)`, whose power order has index ideal `(A^5)`, and which is
globally nonmonogenic. If

```text
Q_2=Spec k[A,C,W]/(W^2-Delta)                              (4)
```

is the quadratic resolvent, then

```text
Q_2^*=k^*,                       Cl(Q_2)=Z^3 direct-sum Z/3. (5)
```

Thus the quintic condition clears the positive-genus obstruction without
losing normality, nonmonogenicity, or the global three-class. The six-branch
boundary obstruction proved in THM-3917 already shows that this candidate
is not a Keller map. The reserved THM-3924 candidate prospectively gives a
second, primitive-class obstruction, but no claim from that pending audit is
used in the present theorem. The conjecture `JC(2)` remains **OPEN**.

## 1. The decic normalization and the hidden THM-3915 family

Put

```text
G_b(t)=t^8+(12b+3)t^6+(48b^2+12b)t^4
       +(64b^3-48b^2-36b-8)t^2
       -192b^3-96b^2-36b-6.                               (6)
```

At `b=1/4`, this is exactly the polynomial `G` of THM-3915:

```text
G_(1/4)=t^8+6t^6+6t^4-19t^2-24.                           (7)
```

Normalize the conic `v^2=u^2-1` by

```text
u=(t+t^-1)/2,                 v=-(t-t^-1)/2,
z=(u^2+b)v.                                               (8)
```

On the ramification curve `z^2=p`, equation `(2)` becomes

```text
A=tG_b/512,                     C=(t^2+1)G_b/1024.         (9)
```

After the harmless target scaling `(a,c)=1024(A,C)`, the homogeneous map is

```text
[T:S] |-> [2TS G_h:(T^2+S^2)G_h:S^10].                   (10)
```

It has no base point. Its slope is

```text
c/a=(t+t^-1)/2.                                           (11)
```

Two generic parameters with the same slope are equal or reciprocal. The
reciprocal difference computed below is not identically zero, so `(10)` is
generically one-to-one. Its image therefore has degree ten and `(10)` is
the normalization of an irreducible rational decic.

There is one preimage of the infinity line, namely `S=0`, mapping to
`[0:1:0]`. In the local parameter `s=S/T`,

```text
a/c=2s/(1+s^2)=2s+O(s^3),
Z/c=s^10/((1+s^2)G_h(1,s))=s^10+O(s^12).                 (12)
```

Thus the infinity point is smooth and its line contact is exactly ten.

## 2. Why K=0 changes exactly two origin branches

Before imposing `K=0`,

```text
disc_t(G_b)=
 -17915904(4b+1)(8b^2+2b+1)
  (64b^3+48b^2+24b+5)^2 K(b)^2.                           (13)
```

Write `G_b(t)=G_x(t^2)`. In `Q[b]/(K)[x]`, define

```text
xi=(3456b^4+15840b^3+8352b^2+998b-115)/150,

q_1=(2/75)(1728b^4+7920b^3+4176b^2+949b+55),

q_0=-(1/75)(5760b^4+26592b^3+13664b^2+3142b+365).         (14)
```

Then the exact factorization is

```text
G_x=(x-xi)^2(x^2+q_1x+q_0).                               (15)
```

The remaining collision tests are

```text
(x^2+q_1x+q_0)|_(x=xi)
 =-(8064b^4+36192b^3+16192b^2+2894b+205)/30,

disc(x^2+q_1x+q_0)
 =(16128b^4+72384b^3+33104b^2+7588b+815)/45.              (16)
```

Their numerators, the numerator of `q_0`, and `xi` are all nonzero modulo
the irreducible quintic `K`. Hence `G_b` has precisely

```text
two double roots +/-sqrt(xi) and four simple roots.        (17)
```

They are nonzero, avoid `+/-1`, and no two are reciprocal. The last fact is
certified by

```text
Res_t(G_b,t^8G_b(1/t))=
 2^18 b^6(48b^2+16b+3)^12
       (64b^3+48b^2+24b+5)^6 !=0.                         (18)
```

Every root of `G_b` maps to the target origin. A simple root gives a smooth
branch with tangent direction `(2t,t^2+1)`. At a double root `t_0`, write
`gamma=(2tG_b,(t^2+1)G_b)`. Then

```text
gamma''=G_b''(2t_0,t_0^2+1),
gamma'''=G_b'''(2t_0,t_0^2+1)+3G_b''(2,2t_0),

det(gamma'',gamma''')=6(G_b'')^2(t_0^2-1) !=0.            (19)
```

Thus each double root gives an `A2`-cuspidal branch. Equation `(18)` makes
all six tangent directions distinct. Their multiplicities are

```text
2,2,1,1,1,1.                                              (20)
```

Pairwise intersections therefore contribute `26`, while the two cusps
contribute intrinsic delta `1+1`. This proves the origin packet and
`delta=28` in `(3)`.

## 3. The two external cusps, two A5 contacts, and completeness

For `a=2tG_b,c=(t^2+1)G_b`, simultaneous derivative vanishing away from
`t=+/-1` forces `G_b=G_b'=0`, because the coefficient determinant in
`(G_b,G_b')` is `2(1-t^2)`. These are exactly the two cuspidal origin
branches already counted.

At `t=+/-1`, one has

```text
(G_b+G_b')|_(1)=0,               (G_b-G_b')|_(-1)=0.      (21)
```

The second--third derivative determinant at either address is

```text
12288(b+1)^3(64b^3+48b^2+24b+5) !=0,                     (22)
```

so these give two further `A2` cusps away from the origin.

All multiple-preimage singularities are reciprocal. Put

```text
J_b=t^4+(4b+2)t^2+1.
```

The exact collision identity is

```text
t^9(a(t)-a(1/t))=
 2(t-1)^3(t+1)^3 J_b(t)^3.                                (23)
```

Moreover

```text
disc(J_b)=4096b^2(b+1)^2,
Res(G_b,J_b)=16b^2(48b^2+16b+3)^4,                        (24)
J_b/t^2=4(u^2+b).
```

Thus the four simple roots of `J_b` avoid all critical and origin
parameters and form two reciprocal pairs. Since `du/dt !=0` there,
equation `(23)` says that each pair consists of two immersed branches with
intersection multiplicity three. These are the two `A5` contacts in `(3)`.

The arithmetic genus of a decic is `36`, and

```text
36=28+(1+1)+(3+3).                                        (25)
```

Every singular point of a parametrized curve comes from a critical
parameter or a multiple preimage. The preceding classifications and `(25)`
therefore prove that `(3)` is complete.

## 4. The six-address identity

The conceptual core is the product formula

```text
r((t+t^-1)/2)=
 -G_b(t) t^8G_b(1/t)/(2^16 t^8).                          (26)
```

At a root of `G_b`, the reciprocal factor is nonzero by `(18)`. Hence the
multiplicity of `r` at the corresponding `u` is exactly the multiplicity of
`G_b` at `t`. The double residual root `x_0` of THM-3917 and the double
origin root `xi` satisfy the exact relation

```text
x_0=(xi+2+xi^-1)/4.                                       (27)
```

On the ramification normalization `(8)`, `F=0` is equivalent to `A=0`, so
its finite support is `G_b=0`. Equations `(17),(18),(26)` prove the promised
one-for-one dictionary:

```text
four simple roots of r  <-> four simple G_b roots
                         <-> four smooth origin branches,

two double roots of r   <-> two double G_b roots
                         <-> two cuspidal origin branches. (28)
```

This explains geometrically, rather than only numerically, why the genus
collapse creates exactly the six-branch boundary failure of THM-3917.

## 5. The normal nonmonogenic cubic order

Return to the natural coordinates `(A,C)` in `(2)` and homogenize in the
slope:

```text
P=A^6p(C/A),                    H=A^9h(C/A),
Q=2A^10-H,
Delta=A^8r(C/A)+4H-4A^10,

f(Z)=Z^3-3PZ-2Q.                                           (29)
```

These are polynomials, and

```text
disc_Z(f)=108A^10 Delta.                                  (30)
```

The parametrization `(9)` lies on `Delta=0`; since its image is an
irreducible degree-ten curve and `Delta` also has degree ten, `Delta` is its
irreducible defining equation up to a nonzero scalar.

The rational chart gives `k(A,C)=k(u,F)` and `A=F/4`. Since the rational
function `F(u,z)` has degree three in `z`,

```text
[k(u,z):k(A,C)]=3,                                        (31)
```

and `f` is irreducible.

There is a universal global overorder. Define

```text
c=b-1/2,                     d=(4b+1)/8,
L=C^3+cA^2C,                 N=CL-dA^4,

e_1=(Z^2+LZ-2L^2)/A^4,

e_2=(2CL^2-4dA^4L-NZ-CZ^2)/A^5.                          (32)
```

The exact cubic reduction in the companion verifies that every coordinate
of

```text
e_1^2,                       e_1e_2,                       e_2^2
```

in the basis `(1,e_1,e_2)` lies in `k[A,C]`. The only displayed coefficient
denominators divide the nonzero constant `4b+1`. Also `Z,Z^2` have polynomial
coordinates in this basis, and

```text
det((1,e_1,e_2)/(1,Z,Z^2))=-(4b+1)/(8A^5).                (33)
```

Thus

```text
B=k[A,C]1+k[A,C]e_1+k[A,C]e_2                             (34)
```

is a genuine finite-free cubic overorder with power-order index ideal
`(A^5)`. Its discriminant is

```text
disc(B/k[A,C])=(27/16)(4b+1)^2 Delta.                     (35)
```

For an independent local check, at the generic point of `A=0` put

```text
m=L-dA^4/C,                         Z=m+A^5x.              (36)
```

After division of `f` by `A^10`, the residual quadratic is

```text
384C^4x^2-512C-96b^3-48b^2-18b-3=0.                      (37)
```

It is separable and irreducible over `k(C)`: after division by the square
`C^4`, its numerator has a simple zero as a rational function of `C`.
Hence `A=0` is unramified and `(33)` removes exactly the five units of local
index debt.

Away from `Delta`, `(35)` is a unit over every height-one DVR. At the
generic point of the reduced irreducible divisor `Delta`, its exponent is
one, so no further overorder can lower it by a positive even index length.
Therefore `(34)` is `R1`. It is finite free over the regular surface
`k[A,C]`, hence `S2`; Serre's criterion makes it normal. Consequently `(34)`
is the full integral closure and its fraction field is the rational field
`k(u,z)`.

The binary index form is also exact. If

```text
I(U,V)=det(1,Ue_1+Ve_2,(Ue_1+Ve_2)^2),                    (38)
```

its four coefficients are

```text
[U^3]  =-A(-256A^2+48Cb^2+24Cb+3C)/(8(4b+1)),

[U^2V] =-3(256A^2C+32A^2b^3+8A^2b^2
             -16C^2b^2-8C^2b-C^2)/(8(4b+1)),

[UV^2] =3AC(512C+96b^3+32b^2+6b+1)/(16(4b+1)),

[V^3]  =-(64A^2b^3+48A^2b^2+12A^2b+A^2+2048C^3
            +384C^2b^3+192C^2b^2+72C^2b+12C^2)
             /(64(4b+1)).                                 (39)
```

They have gcd one (equivalently, no common irreducible factor), but their
coefficient ideal is a proper ideal contained in `(A,C)`: all four vanish at
`(A,C)=(0,0)`. Hence `I(U,V)` can never be a unit for polynomial `U,V`; no
element generates `B`. This proves global nonmonogenicity.

## 6. The conic bundle and persistent three-class

The quadratic surface `Q_2` is normal. It is a hypersurface and hence `S2`;
the irreducible branch curve has only the finite singular packet `(3)`, so
the double surface is regular in codimension one. Serre's criterion applies.

On a pencil line `C=sA`, equation `(29)` is

```text
Delta(A,sA)=-4A^8(A^2-Ah(s)-r(s)/4).                      (40)
```

The residual quadratic discriminant is the exact cube

```text
h^2+r=p^3=(s^2-1)^3(s^2+b)^6.                             (41)
```

Thus the smooth conic-bundle model has the same four degenerate-fibre
exponents as THM-3915:

```text
3,3,6,6.                                                   (42)
```

The only geometric change occurs in the inverse image `D` of the blowup-
origin exceptional line:

```text
D: V^2=r(s).                                               (43)
```

Its arithmetic genus is still three and `D^2=-2`, but `(26)` shows that it
now has two nodes and normalization genus one. The ambient conic surface is
smooth at either node, because there `h !=0` and the derivative of

```text
V^2-r(s)-4Ah(s)+4A^2
```

with respect to `A` is `-4h !=0`.

There is also a useful exact genus gap. In both the generic THM-3915 packet
and the present degeneration, let `R_odd` be the number of distinct finite
roots of `r` having odd multiplicity. Infinity is unramified in both the
quadratic curve `D` and the normalized cubic curve `F=0`; the odd roots give
the simple normalized ramification. Riemann--Hurwitz therefore gives

```text
g(D^nu)=(R_odd-2)/2,             g(F^nu)=(R_odd-4)/2,
g(D^nu)=g(F^nu)+1.                                      (43a)
```

Thus THM-3915 has genera `(3,2)`, while `K=0` has genera `(1,0)`. The two
double roots remove two branch pairs simultaneously; the arithmetic genus
of `D` stays three because the loss is absorbed by its two nodes.

Blowing up the two boundary nodes replaces the old origin block `[-2]` by

```text
R_node=[-10  2  2]
       [  2 -1  0]
       [  2  0 -1].                                      (44)
```

This is integrally congruent to `diag(-2,-1,-1)` by the determinant-one
change `D_strict=D_total-2E_1-2E_2`. Killing the strict transform and the
two new exceptional curves is therefore equivalent to killing the old
class `D_total`.

The rest of the removed lattice is unchanged:

```text
R=R_node orthogonal_sum K_10 orthogonal_sum A2(-)^2
             orthogonal_sum A5(-)^2,

K_10=[-4  5].                                              (45)
     [ 5 -4]
```

It has determinant `5832` and nonunit Smith factors

```text
3,3,6,6,18.                                                (46)
```

In the unimodular conic-bundle Picard basis, killing the two infinity
sections, `D`, the middle ADE curves, and the two node exceptional curves
again leaves four far fibre endpoints with the sole relation

```text
3L_1+3L_2+6L_3+6L_4=0.                                   (47)
```

Consequently

```text
Cl(Q_2)=Z^4/<(3,3,6,6)>=Z^3 direct-sum Z/3.               (48)
```

Nondegeneracy of `(45)` makes the boundary-divisor kernel zero, so every
unit is scalar. This proves `(5)`.

The surviving three-class is the actual cyclic-layer class, not merely an
available abstract slot. Indeed, `(30)` and irreducibility of `Delta` show
that the cubic discriminant is nonsquare, while `(31)` makes the cubic
connected. Its Galois closure therefore has group `S_3`, and over `Q_2` the
alternating layer is a connected `C_3` extension unramified at every
codimension-one point of the affine surface `Q_2` (equivalently, every
divisorial valuation centred on `Q_2`). Kummer theory writes it as adjoining
a cube root of some `g` with `div(g)=3E`. If `[E]` were zero, then `g` would
be a cube times a unit; scalar units over the algebraically closed field are
cubes, contrary to connectedness. Hence `[E]` is the nonzero element of
`Cl(Q_2)[3]`, and so generates the unique `Z/3` summand in `(48)`.

The comparison is now exact. The `K=0` deformation makes the collapsed
divisor rational while retaining a normal rational nonmonogenic cubic order
and the needed global three-class. It nevertheless fails: the same
double-root mechanism turns two smooth origin branches into cuspidal ones
and sends all six ramification addresses through one point. THM-3917 proves
that this six-branched ramification curve cannot lie in the boundary of an
affine-plane Keller atlas. The reserved THM-3924 candidate identifies an
additional class-divisibility mechanism, explicitly pending its own audit.
No planar Jacobian counterexample is obtained.
