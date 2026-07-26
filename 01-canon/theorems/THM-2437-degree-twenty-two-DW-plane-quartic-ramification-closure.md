---
id: THM-2437
title: "Degree-twenty-two D-W plane quartic ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane B=C=E=0 is empty. After the choice D=delta^2,
  the invariant ratio lambda=W^2/D^3 and coordinates p=delta/y^2,
  v=u/y^2 turn the first two fluxes into an absolutely irreducible quartic.
  Uniform irreducibility follows from its squarefree odd-degree constant
  term and its degree-two cubic coefficient. Its discriminant is
  wall^4 L_5 K_9(v,lambda). A complete exceptional-parameter audit leaves
  at least twelve simple finite branch values away from the first-flux
  wall in every nonzero fibre, hence normalization genus at least three.
  No rational trajectory survives. This closes the D-W plane, not other
  mixed strata, degree twenty two, JC(2), or DC(2).
source: codex-2026-07-26-degree-twenty-two-dw-plane
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
related:
  - THM-2320-degree-eighteen-dw-ratio-bank-closure
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
script: 04-computation/jc2_degree22_dw_plane_quartic_ramification_thm2437.py
output: 05-knowledge/results/jc2_degree22_dw_plane_quartic_ramification_thm2437.out
script_sha256: 6fdda40f03dde5cb38cb74a5932f482d2f959159c28e334dde9b225d663d39c9
output_sha256: dfd20f71680acda40d0c5af27125796e3ce4340e78033f0b62275a3f838acba8
hash_basis: working-tree bytes (LF)
---

# THM-2437 -- the degree-twenty-two D-W plane is empty

**PROVED + VERIFIED-EXACT.**

THM-2429 closes the `C,W` plane by a complete hyperelliptic fibre
classifier. The next weighted plane is not hyperelliptic: its quotient is
quartic in the scale coordinate. The odd degree of its constant term gives a
different and stronger irreducibility mechanism, while exact branch counting
survives every parameter collision.

The conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=C=E=0
    => contradiction.                                           (1)
```

Thus this is the second complete support-at-most-two coefficient plane in
the open degree-twenty-two chart.

## 1. Weighted quotient and the choice-free ratio

Use the target-translated coordinates and weights of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `D=0` and `W=0` are already the full `W`-axis of THM-2423 and
the `D`-axis of THM-2425. Hence suppose `D,W!=0`. Since the constant field
`C` is algebraically closed, choose

```text
delta in C*,                  D=delta^2,

mu=W/delta^3,                 lambda=mu^2=W^2/D^3.    (3)
```

Changing `delta` to `-delta` changes `(p,mu)` to `(-p,-mu)` below, so the
curve is unchanged. The weighted projective ratio `lambda` is choice-free.

First take `y!=0` in `C(x)` and put

```text
v=u/y^2,             zeta=Z/y^3,             p=delta/y^2.
                                                               (4)
```

Then

```text
D/y^4=p^2,                 W/y^6=mu p^3.              (5)
```

Dividing the first two fluxes `N_1=N_2=0` of THM-2411 by `y^5,y^6`
gives

```text
f_1
 =11979(7-121v)zeta
  +4(922383v^2-25410v+63+511104p^2)
 =0,                                                        (6)

f_2
 =15944049zeta^2-162339408zeta v+2236080zeta
  -1190488992v^3+147581280v^2-1219680v+672
  +(-1978994688v+16355328)p^2
  -1319329792mu p^3
 =0.                                                        (7)
```

The open chart says

```text
121v-7!=0.                                               (8)
```

Thus (6) reconstructs `zeta` uniquely in the function field.

## 2. The exact quartic

Define

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583,                           (9)

A_2(v)
 =-1810903826688v^3+119729178624v^2
   +4329050880v-109716992.                              (10)
```

Exact elimination gives

```text
Res_zeta(f_1,f_2)=2295943056 R_mu(v,p),                 (11)

R_mu
 =29025255424p^4
  -82458112mu(121v-7)^2p^3
  +A_2(v)p^2
  -7L_5(v).                                             (12)
```

There is no `p` term. The leading `p^4` coefficient is a nonzero constant,
and

```text
deg A_2=3,             deg L_5=5,

gcd(L_5,L_5')=1,       L_5(7/121)=-44800.              (13)
```

No root choice or squaring entered (11).

## 3. Uniform absolute irreducibility

For every `mu!=0`, the polynomial `R_mu` is absolutely irreducible in
`C[v,p]`.

Divide by the constant leading coefficient and write

```text
P=p^4+b(v)p^3+c(v)p^2+d(v),                            (14)
```

where `d` is a nonzero constant multiple of the squarefree quintic `L_5`,
`deg c=3`, and

```text
b=constant*mu(121v-7)^2,              deg b=2.         (15)
```

Because (14) is monic in `p`, any factorization over `C(v)` descends by
Gauss to monic factors in `C[v,p]`.

### 3.1 No linear-times-cubic factor

Suppose

```text
P=(p+r)(p^3+sp^2+tp+q).                                (16)
```

The `p` and constant coefficients give

```text
q+rt=0,                     rq=d,

d=-r^2t.                                               (17)
```

Squarefreeness of `d` forces `r in C*`. Then `t=-d/r^2` has degree five.
But the `p^3` coefficient gives `s=b-r`, of degree at most two, so the
`p^2` coefficient `t+rs` has degree five, contradicting `deg c=3`.

### 3.2 No quadratic-times-quadratic factor

Suppose

```text
P=(p^2+ap+r)(p^2+cp+s).                                (18)
```

Then

```text
rs=d,                         as+cr=0.                 (19)
```

The squarefree product makes `r,s` coprime. Hence

```text
a=rh,                         c=-sh
```

for some `h in C[v]`, and the cubic coefficient is

```text
b=h(r-s).                                               (20)
```

Since `deg r+deg s=5`, the two degrees cannot be equal and

```text
deg(r-s)=max(deg r,deg s)>=3.
```

Equation (15) makes `h!=0`, so (20) has degree at least three, contradicting
`deg b=2`. The two possible factor-degree partitions are exhausted; (12) is
absolutely irreducible.

The boundary `mu=0` is exactly the already closed `D`-axis. It is a genuine
structural degeneration, not a gap in this argument.

## 4. Discriminant and complete collision divisor

The quartic discriminant is

```text
Disc_p R_mu
 =-2^36*7*11^18 (121v-7)^4 L_5(v) K_9(v,lambda),      (21)
```

which defines the exact primitive normalization of `K_9` used here. It has
degree nine generically and degree two in `lambda`. Its complete root and
degree-drop controls are

```text
LC_v(K_9)
 =3302590884214385499714 lambda(231lambda+128),         (22)

Disc_v(K_9)
 =2^155 3^30 5^8 7^11 11^148
    lambda^4 S_6(lambda)^3 S_7(lambda),                 (23)

Res_v(K_9,L_5)
 =nonzero_constant*Q_5(lambda),                         (24)

K_9(7/121,lambda)
 =10276044800(385lambda+512).                           (25)
```

Here

```text
S_6
 =71778115591875lambda^6
  -10643296267296000lambda^5
  -45431893495296000lambda^4
  -70165361991352320lambda^3
  -40747428281843712lambda^2
  -589081608192000lambda
  +4306718326521856,                                    (26)

S_7
 =37565749714841455078125lambda^7
  +27921686965239375000000lambda^6
  -13643584211703600000000lambda^5
  -11118982374162432000000lambda^4
  -2775896567285022720000lambda^3
  +3185830828706176696320lambda^2
  +998410901657688735744lambda
  -455539102308525670400,                               (27)

Q_5
 =4932186875lambda^5
  -1123257520000lambda^4
  -5423878240000lambda^3
  -7375509504000lambda^2
  -1274053017600lambda
  +1520839950336.                                       (28)
```

The three displayed parameter polynomials are irreducible over `Q`. The six
polynomials

```text
lambda,  231lambda+128,  385lambda+512,  S_6, S_7, Q_5
```

are pairwise coprime. Thus the exceptional cases below do not overlap.
Multiplying (23) by (22) shows why the complete homogeneous collision
divisor contains

```text
lambda^5(231lambda+128)S_6^3S_7:
```

the ordinary degree-nine discriminant alone does not record a root escaping
to infinity at `231lambda+128=0`.

## 5. Every exceptional fibre retains twelve simple branches

At `lambda=0`,

```text
K_9(v,0)=2^14 H_D(v)^2,                                (29)
```

where `H_D` is the squarefree quartic of the closed `D`-axis in THM-2425.
For `lambda!=0`, equations (22)--(28) give an exhaustive classifier.

| parameter | exact `K_9` behaviour | simple finite branch values counted away from `121v-7=0` |
|---|---|---:|
| generic | degree nine, squarefree, disjoint from `L_5` and the wall | `5+9=14` |
| `231lambda+128=0` | squarefree degree-eight polynomial, disjoint from `L_5` and the wall | `5+8=13` |
| `S_6(lambda)S_7(lambda)=0` | one double root and seven simple roots, disjoint from `L_5` and the wall | `5+7=12` |
| `Q_5(lambda)=0` | squarefree `K_9` and `L_5` share exactly one simple root | `4+8=12` |
| `385lambda+512=0` | one simple `K_9` root is the wall; the other eight are simple and disjoint from `L_5` | `5+8=13` |

For the degree drop, exact substitution gives a nonzero scalar times the
squarefree octic

```text
3721928118949345041v^8-963805077634265658v^7
-251896538372380902v^6+27089409377228814v^5
-719579612757852v^4+5089117224114v^3
+35886232998v^2-398179782v-1740725.                    (30)
```

For a root of `S_6` or `S_7`, the exact subresultant sequence of
`K_9,partial_v K_9` has degrees

```text
9,8,7,6,5,4,3,2,1,0,
```

and only the constant subresultant vanishes modulo the relevant parameter
polynomial. Hence the gcd has degree exactly one: there is one double root,
not a triple root or two double roots. The exponent three on `S_6` in (23)
records tangency of the parameter line to the discriminant, not extra
multiplicity in the special fibre.

For a root of `Q_5`, the `K_9,L_5` subresultants have degrees

```text
9,5,4,3,2,1,0,
```

and again only the constant subresultant vanishes. Thus there is exactly one
common root. At the wall collision,

```text
partial_v K_9(7/121,-512/385)=954932291174400!=0.       (31)
```

The table deliberately ignores every double or shared discriminant value
and every point on the first-flux wall. Only simple zeros of the quartic
discriminant away from that wall are counted.

## 6. Uniform genus and constant-map closure

Let `C_lambda` be the smooth projective normalization of (12). Section 3
makes it a connected degree-four cover of the `v`-line. At a simple zero
of the discriminant, the constant leading coefficient and characteristic
zero imply one simple ramification point. The table gives

```text
sum_P(e_P-1)>=12.                                      (32)
```

No infinity or discarded collision contribution is needed. Riemann--Hurwitz
gives

```text
2g(C_lambda)-2
 =4*(-2)+sum_P(e_P-1)
 >=4,

g(C_lambda)>=3.                                        (33)
```

A nonconstant pair `(v,p) in C(x)^2` satisfying (12) would therefore give a
nonconstant morphism from `P^1` to a positive-genus curve, impossible.
Hence `v` is constant. Since the leading coefficient in (12) is nonzero,
the fixed-`v` quartic is never an identity in `p`; thus `p` is algebraic
over the algebraically closed constant field and is constant.

Now `p!=0` and

```text
y^2=delta/p.
```

Hence `y`, then `u`, are constant. Equation (8) and the first flux reconstruct
constant `zeta`, so `Z`, nonzero `T`, and `q` are constants. The genuine deck
fixes the constant field but sends `q` to `-q`, contradicting `q!=0`.

It remains to justify division by `y`. At `y=0`, the open chart gives
`u!=0`, while the first flux is

```text
N_1=-1449459uZ.
```

It forces `Z=0`, contradicting `T!=0`. This proves (1).

## 7. Scope and structural lesson

The map and loss ledger is

```text
source:
  the full D-W coefficient plane in the open first-flux chart;

map:
  choose D=delta^2 and pass to (v,p,mu), retaining
  lambda=mu^2=W^2/D^3;

preserved:
  the two exact fluxes, the first-flux wall, the rational trajectory,
  irreducible normalization, and finite branch indices;

choice forgotten:
  delta versus -delta;

restoration sidecar:
  (p,mu)->(-p,-mu), which leaves lambda and the quartic isomorphism class
  unchanged;

hostile controls:
  lambda=0 gives the doubled D-axis quartic, S_6 gives a tangent double
  root, Q_5 collides the two branch divisors, and lambda=-512/385 moves
  one branch exactly onto the excluded wall.                          (34)
```

The useful structural object is not merely the discriminant. It is the
quartic together with its odd-degree squarefree constant fibre. That fibre
rules out every disconnected cover before branch counting begins. The
homogeneous parameter discriminant, the `L_5` resultant, and the wall
evaluation then form three separate sidecars: root collision, cross-factor
collision, and collision with the excluded reconstruction wall.

This theorem closes exactly the `D,W` plane, including its axes. Together
with THM-2429, two of the ten coefficient planes in `mathcal A!=0` are now
empty. The other eight planes, higher-support strata, branches outside the
inherited reduction, split/even short edges, and integral order raising
remain open. Nothing here proves degree twenty two, `JC(2)`, or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_dw_plane_quartic_ramification_thm2437.py
python3 -O 04-computation/jc2_degree22_dw_plane_quartic_ramification_thm2437.py
```

The companion checks the weighted flux reduction, exact resultant (11),
quartic discriminant (21), complete parameter discriminant and leading
coefficient (22)--(23), cross-resultant (24), wall evaluation (25), pairwise
separation of every exceptional locus, the degree-drop octic, both
subresultant classifiers, the simple wall collision, the uniform branch
floor, and the `y=0` hostile boundary. Every truth-bearing check uses explicit
exceptions and remains active under optimized Python.

Normal, optimized, and stored transcripts byte-match. The declared hashes
are over working-tree bytes. **QED.**
