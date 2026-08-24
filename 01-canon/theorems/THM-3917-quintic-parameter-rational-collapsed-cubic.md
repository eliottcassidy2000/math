---
id: THM-3917
title: "Quintic-parameter rational collapsed cubic and six-branch boundary obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the
  irreducible quintic parameter locus displayed below, the polynomial-part
  construction produces an irreducible cubic whose collapsed divisor has
  genus zero: its degree-three projection has four finite ramifications,
  two unramified nodes, and three unramified infinity places. Thus it
  genuinely escapes THM-3916's positive-genus obstruction. Nevertheless,
  the rational ramification divisor has six distinct normalization points
  above the contracted divisor, and all six map to one point of the finite
  normalization. Its image is therefore at least six-branched there. Since
  every irreducible boundary curve of an affine-plane open in a normal
  surface is unibranch, the associated cubic normalization admits no
  same-field affine-plane Keller atlas. This closes the candidate, not
  JC(2).
source: jc_degree6_one_place / post-THM-3915 genus-zero design lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root and jc_zero_debt_lift, 2026-08-23).
  The audits independently recovered the residual quartic, its unique
  double root, the six-point support and p--h avoidance, both cubic-root
  descents, the four-ramification Riemann--Hurwitz ledger, and the split
  infinity expansions. They also checked the normalization factorization
  B subset k[u,z], generic birationality of the ramification curve, the
  contraction of all six addresses to one point, and the relative common-
  resolution proof of boundary unibranchness, including the parallel-edge
  case in the dual multigraph. Normal and optimized runs byte-match the
  frozen output and all 55 active gates pass; no repair was required.
related:
  - THM-3916-positive-genus-collapsed-valuation-keller-obstruction
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
script: 04-computation/jc2_quintic_parameter_rational_collapsed_cubic_thm3917.py
output: 05-knowledge/results/jc2_quintic_parameter_rational_collapsed_cubic_thm3917.out
script_sha256: 66e54dfe3cf1365a18f0df784e1dad9987c0ea74077da219edb5f817a16a583f
output_sha256: e8c8b9802d48d0d1a4ba70afba78bfda1e0081a9d1d0bd79a0ba4da9c29567a7
semantic_sha256: 593e5af9e8099c8a5f3397a2247ba166235dd3710e69296eac3aa63628629926
hash_basis: raw LF bytes
---

# THM-3917 -- genus zero is attainable, but six branches cannot be plane boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Choose `b in k`
satisfying

```text
K(b)=2304b^5+10176b^4+4064b^3+996b^2+84b+5=0.             (1)
```

Put

```text
p=(u^2-1)(u^2+b)^2,                                       (2)

h=u^9+(3b-3/2)u^7+(3b^2-9b/2+3/8)u^5
    +(b^3-9b^2/2+9b/8+1/16)u^3
    +(-3b^3/2+9b^2/8+3b/16+3/128)u,                      (3)

r=p^3-h^2,                 F(u,z)=z^3-3pz+2h.             (4)
```

Then `h` is the polynomial part at infinity of the branch of `p^(3/2)`
asymptotic to `u^9`. The cubic `F` is irreducible. If `C_F` denotes the
normalization of `F=0`, then

```text
C_F is P1 over k,             deg(u:C_F -> P1)=3.          (5)
```

Its finite discriminant packet consists of four simple roots and two double
roots. The simple roots give four tame transposition ramifications. Each
double root is an ordinary node of the plane model and gives two unramified
normalization branches. There are three further unramified places over
`u=infinity`. In particular

```text
total ramification=4,          2g(C_F)-2=-6+4=-2.           (6)
```

This is an exact genus-zero escape from the obstruction of THM-3916. It is
not a planar Jacobian counterexample. In the natural rational chart the
ramification curve is rational, but its image in the finite normalization
has at least six branches at one point. The boundary-unibranch lemma below
therefore rules out every same-field affine-plane Keller atlas.

## 1. The exact repeated-residual packet

The reduction of `K` modulo `11` is irreducible of degree five. Hence `K`
is irreducible over `Q`. Its discriminant is

```text
disc(K)=1170238183833600000000000 !=0.                     (7)
```

In particular all five roots are simple. The quintic is coprime to each of

```text
b, 4b+1, 8b^2+2b+1,
48b^2+16b+3, 64b^3+48b^2+24b+5.                           (8)
```

These are precisely the leading, common-zero, and higher-collision factors
needed below.

Let

```text
S(x)=16384 r(sqrt(x)).                                     (9)
```

This is a quartic. Before imposing `(1)`, its discriminant is

```text
disc_x(S)=
 -27*2^14 (48b^2+16b+3)^6
          (64b^3+48b^2+24b+5)^3 K(b).                    (10)
```

There is a more precise factorization. In `Q[b]/(K)[x]`, put

```text
x_0=(95616b^4+421728b^3+162368b^2+16726b+965)/15120,

q_2=384(4b+1)(8b^2+2b+1),

q_1=(32/5)(4224b^4+1472b^3+192b^2-156b-55),

q_0=-(16/45)(188352b^4+111696b^3+34556b^2+5272b+215).
                                                                    (11)
```

Then

```text
S(x)=(x-x_0)^2(q_2x^2+q_1x+q_0).                          (12)
```

The two remaining collision tests are

```text
(q_2x^2+q_1x+q_0)|_(x=x_0)
 =-(2/3)(35712b^4+41952b^3+13408b^2+2438b+85),

disc(q_2x^2+q_1x+q_0)
 =-(1024/27)(114974016b^4+47144272b^3
             +11949168b^2+978444b+58837).                 (13)
```

Both displayed numerators, as well as the numerator of `q_0`, are coprime
to `K`. Thus `(12)` has exactly one double root and two simple roots, all
nonzero. Moreover

```text
S(1)=-(64b^3+48b^2+24b+5)^2,
S(-b)=b(48b^2+16b+3)^2,                                  (14)

Res_u(p,h)=-b^2(48b^2+16b+3)^4
                 (64b^3+48b^2+24b+5)^2/2^42 !=0.         (15)
```

Consequently no root of `r` is a root of `p`. Over `k`, `(12)` becomes

```text
div(r)=2[sqrt(x_0)]+2[-sqrt(x_0)]
       +[sqrt(x_1)]+[-sqrt(x_1)]
       +[sqrt(x_2)]+[-sqrt(x_2)],                         (16)
```

with all six support points distinct.

At a root `u_0` of `r`, choose a local square root `s^2=p` with `h=s^3`.
Writing `z=s+v` gives

```text
F=v^2(v+3s)+2(h-s^3),
r=(s^3-h)(s^3+h).                                         (17)
```

All omitted factors are units by `(15)`. Thus the double sheet is locally
`V^2=unit*r`. Odd order of `r` gives one simple ramification; order two
gives an ordinary node with two unramified branches. This proves the finite
part of `(6)`, including normal rather than merely power-basis ramification.

## 2. Irreducibility and the reducible boundary

Because `F` is monic in `z`, a root in `k(u)` would lie in `k[u]`. Degree
comparison forces such a root to have degree exactly three. Write it as

```text
c_3u^3+c_2u^2+c_1u+c_0.                                  (18)
```

The coefficient of `u^9` in `F(u,(18))` is

```text
(c_3-1)^2(c_3+2).                                         (19)
```

For `c_3=-2`, the coefficients of `u^8,u^7,u^6` successively force

```text
c_2=0,                 c_1=1-2b,                 c_0=0.    (20)
```

The `u^5` coefficient is then `-9(4b+1)/4`, so `b=-1/4`.
But at that value the remaining `u` coefficient is `27/64`, a
contradiction.

For `c_3=1`, the `u^7` and `u^5` coefficients are

```text
3c_2^2,                 (3/4)(2b-2c_1-1)^2.               (21)
```

They force `c_2=0,c_1=b-1/2`. The `u^3` coefficient is then `3c_0^2`,
and after `c_0=0` the `u` coefficient is

```text
3(4b+1)^2/64.                                              (22)
```

Thus the only possible reducible parameter is again `b=-1/4`. Since

```text
K(-1/4)=81/4,                                              (23)
```

it is disjoint from `(1)`. This proves irreducibility of `F` and also
identifies the exact reducibility boundary that a parameter search must
exclude.

## 3. The three infinity places

Put

```text
t=1/u,                 y=z/u^3,
P(t)=t^6p(1/t)=(1-t^2)(1+bt^2)^2,
H(t)=t^9h(1/t),
s(t)=(1+bt^2)sqrt(1-t^2),
D=(4b+1)(8b^2+2b+1).                                    (24)
```

The definition of the polynomial part is the exact expansion

```text
H(t)-s(t)^3=-(3D/256)t^10+O(t^12).                        (25)
```

The leading equation at `t=0` is `(y-1)^2(y+2)`. At its simple root,
substitution `y=-2s+t^10q` gives the residual linear equation

```text
9q-3D/128=0.
```

At the double root, substitution `y=s+t^5q` gives

```text
3(q^2-D/128)=0.                                           (26)
```

The factor `D` is nonzero by `(8)`, so Hensel lifting gives the three
distinct expansions

```text
y_0=-2s+(D/384)t^10+O(t^12),
y_+= s+sqrt(D/128)t^5+O(t^7),
y_-= s-sqrt(D/128)t^5+O(t^7).                             (27)
```

All use `t` itself as uniformizer, so all three infinity places are
unramified. The first expansion is defined over `Q(b)`. Hence the genus-zero
normalization has a `Q(b)`-rational point and is already `P1` over `Q(b)`,
not merely after algebraic closure. Riemann--Hurwitz now gives `(5),(6)`.

## 4. The rational ramification curve and its six-point collision

The other divisor in the natural Jacobian is

```text
D_1: F_z/3=z^2-p=0.                                       (28)
```

It is irreducible and rational. Indeed, after

```text
z=(u^2+b)v,
```

its normalization is the conic

```text
v^2=u^2-1.                                                 (29)
```

At `p !=0`, the common point of `F=0` and `D_1` is uniquely `z=h/p`, and
exact elimination gives

```text
Res_z(F,F_z/3)=-4r.                                       (30)
```

Thus `(16)` is the complete six-point support of `F intersect D_1`. At the
four simple roots the intersection is transverse. At either double root,
`F` is nodal and the smooth curve `D_1` crosses both normalization branches,
with total local intersection multiplicity two.

Now define the polynomial chart

```text
A=F/4,                       C=uF/4.                       (31)
```

Its Jacobian is

```text
Jac_(u,z)(A,C)=-F F_z/16.                                 (32)
```

The associated degree-three global cubic is also exact. Put

```text
P=A^6 p(C/A),             H=A^9 h(C/A),
Q=2A^10-H,
Delta=A^8 r(C/A)+4H-4A^10,                                (33)

f(Z)=Z^3-3PZ-2Q.                                          (34)
```

Every expression in `(33)` is a polynomial in `A,C`, and

```text
f(A,uA,A^3z)=A^9(F-4A),
disc_Z(f)=108A^10 Delta.                                  (35)
```

Equations `(31),(35)` identify

```text
k(A,C) subset k(u,z),                  [k(u,z):k(A,C)]=3. (36)
```

Let `B` be the integral closure of `k[A,C]` in `k(u,z)` and let

```text
X=Spec(B).                                                 (37)
```

There is an everywhere-defined birational morphism

```text
phi:A2_(u,z) -> X.                                        (38)
```

Indeed, every element of `B` is integral over `k[A,C]`, hence integral over
`k[u,z]`; since `k[u,z]` is integrally closed, `B subset k[u,z]`.
Conversely, `(34)` is monic, so `Z=A^3z` belongs to `B`.

Let `D subset X` be the image of `D_1`. On `D_1`, equation `(31)` reads
`A=(h-pz)/2`. Away from `A=0` and `p=0`, the formulas

```text
u=C/A,                         z=(h-2A)/p
```

recover its generic point, so `D_1 -> D` is generically birational.
Equation `(32)` shows that `D` is an intrinsic ramification divisor of
`X -> A2_(A,C)`. On the other hand, the irreducible curve `F` maps into the
finite fibre of `(37)` over `(A,C)=(0,0)`. Its image is therefore one point,
say `p`. All six distinct smooth normalization points of `D_1` listed in
`(16),(30)` consequently map to `p`. Since they define six distinct
valuations of `k(D)`,

```text
D has at least six normalization branches at p.            (39)
```

This is the global invoice hidden by the successful genus collapse.

## 5. Boundary curves of an affine plane are unibranch

We use the following general lemma.

> **Boundary-unibranch lemma.** Let `X_0` be a normal surface containing a
> dense open `U isomorphic to A2`. Every irreducible curve
> `D_0 subset X_0 minus U` is unibranch at every point of `X_0`.

To prove it, take a normal proper completion of `X_0` and a common relative
log resolution

```text
Y -> Xbar,                         Y --> P2.               (40)
```

Both maps are the identity on `U`. The reduced SNC boundary `Y minus U` has
a tree as its **dual multigraph**. One way to see this is to resolve the
birational map to the standard pair `(P2,line at infinity)`: blowing up a
smooth boundary point adds a leaf and blowing up a boundary crossing
subdivides an edge, so the one-vertex tree remains a tree.

Suppose `D_0` had two branches at a point `q`. On a sufficiently high log
resolution, the strict transform of `D_0` meets the fibre over `q` in at
least two branch attachments. The fibre is connected because `Xbar` is
normal. A path inside that connected exceptional fibre, together with the
two edges to the strict-transform vertex, makes a cycle in the boundary
dual multigraph. If both attachments meet the same exceptional component,
the two parallel edges already make a multigraph two-cycle. Both cases
contradict the tree property. This proves the lemma.

Finally suppose that the same field `k(u,z)` admitted coordinates `x,y`
such that

```text
k(u,z)=k(x,y),            A,C in k[x,y],
Jac_(x,y)(A,C) in k*.                                      (41)
```

The induced map `U=A2_(x,y) -> A2_(A,C)` is etale and quasi-finite. It
factors through `(37)` by normality, and the resulting birational
quasi-finite map `U -> X` is an open immersion by Zariski Main. The
intrinsic ramification divisor `D` is disjoint from this etale open, so it
is a boundary curve. But `(39)` contradicts the boundary-unibranch lemma.
Therefore no coordinates satisfying `(41)` exist.

The construction has achieved something mathematically useful: repeated
residual roots really can lower the collapsed cubic from genus two to genus
zero without making it reducible. What fails is the next global address
invoice--six branches of one ramification component are forced through one
normalization point. This closes this specific cubic design and gives a new
obstruction for successor searches. It does **not** prove the planar
Jacobian conjecture, which remains **OPEN**.
