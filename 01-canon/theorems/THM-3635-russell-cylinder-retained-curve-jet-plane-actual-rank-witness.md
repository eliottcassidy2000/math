---
id: THM-3635
title: "Russell-cylinder retained-curve jet plane and actual-rank witness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The retained
  one-jet image of the actual restriction ring is the common-value line plus
  the at-most-two-dimensional plane spanned by the c- and e-tangent columns.  A nonzero
  first-stable determinant forces that plane to contain the all-ones vector
  and gives one exact slope equation.  For the minimal rho=0,u=1 Hermite
  polynomial, an exhaustive actual-ring computation constructs a genuine
  rank-two/rank-two t=0 witness and proves degree-94 minimality for the fixed
  pair U=c,V=e under a symmetric x-degree cutoff.  This is not an all-order
  fold pair or a global Keller claim.
source: root / audit_thm3629 / audit_thm3630 actual-ring continuation, 2026-08-21
audit: >
  PASS -- an independent hostile audit rederived the retained jet plane and
  ordinary-triple conductor, independently reconstructed the finite
  module/SAGBI and normalization certificates, and verified the global
  conductor, exhaustive degree-93/94 boundary, witness hashes, and
  byte-identical normal, optimized, and stored 111-gate transcripts.
depends_on:
  - THM-3634-russell-cylinder-quadratic-fold-first-stable-jet-rank-debt
related:
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3632-russell-cylinder-formal-pair-algebraization-triple-fibre-obstruction
  - THM-3634-russell-cylinder-quadratic-fold-first-stable-jet-rank-debt
script: 04-computation/jc2_russell_cylinder_retained_curve_jet_plane_actual_rank_witness_thm3635.py
output: 05-knowledge/results/jc2_russell_cylinder_retained_curve_jet_plane_actual_rank_witness_thm3635.out
script_sha256: d3f56bdbb83720568c651a265a1559a76aebfd96672826bfa7c9ddc4f9ad0002
output_sha256: 93f4c541df72878444ebd84d147d43f97156bbe93fcd46d96afd781ae5d4ed53
hash_basis: raw LF bytes
---

# THM-3635 -- retained-curve jet plane and actual-rank witness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem has
two levels.  First, it identifies the retained one-jet image for every fixed
collision polynomial.  Second, for the smallest polynomial on the
`rho=0,u=1` Hermite stratum, it computes the whole actual restriction ring and
exhibits a constant-determinant first-stable shadow inside that ring.

The second result proves that THM-3634's remaining rank-two/rank-two cell is
genuinely nonempty at `t=0`.  It does **not** construct target functions whose
Jacobian stays constant away from `t=0`.

All rings, completions, and closed points are over `C`.  The exact companion
works over `Q` and then base-changes.

## 0. Compiler, restriction ring, and notation

Use the exponent-two Danielewski compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3),                    c^2e=b(b+4).                (1)
```

Let `Q in C[x]` satisfy

```text
Q(-1)=Q(1)=-3,                    Q(0)=-3/4.             (2)
```

Restrict `(1)` to `q=Q(x)` and write

```text
S_Q=gamma_Q(R)=C[b(x),c(x),e(x)] subset C[x],
z=e+3.                                                       (3)
```

The three points `x=-1,0,1` all map to `(b,c,e)=(0,0,-3)`.  Put

```text
v_-=Q'(-1),                 u=Q'(0),                 v_+=Q'(1). (4)
```

For `H in S_Q`, let

```text
dH=(H'(-1),H'(0),H'(1)) in C^3.                          (5)
```

## 1. The exact retained one-jet plane

Direct differentiation of `(1)` at the retained fibre gives

```text
db=(0,0,0),

tau_c=dc=(2(v_-+6), 3, 2(6-v_+)),
tau_e=de=(-2(v_-+9), 4u, 2(9-v_+)).                     (6)
```

Define

```text
T_Q=span_C{tau_c,tau_e} subset C^3.                      (7)
```

Then the full retained one-jet image is exactly

```text
{ ((h,h,h),d) : h in C, d in T_Q }.                     (8)
```

Necessity follows from the common target value and the chain rule: at the
common point only the `c`- and `e`-partials survive because `db=0`.
Conversely, every prescribed element of `(8)` is realized by

```text
h+lambda c+mu(e+3) in S_Q.                              (9)
```

In particular, for every actual-ring quadruple `U,V,A,B in S_Q`,

```text
rank[dU dV dA dB] <= 2.                                 (10)
```

Suppose now that the THM-3634 determinant shadow holds:

```text
U'B-AV'=kappa in C*.                                    (11)
```

Let `A_0,B_0` be the common retained values of `A,B`.  Evaluating `(11)` at
the three points gives

```text
kappa (1,1,1)=B_0 dU-A_0 dV in T_Q.                    (12)
```

Thus `(1,1,1) in T_Q`.  Moreover `dim T_Q=2`.  Indeed `tau_c` has middle
entry `3`, so `T_Q` is nonzero.  If it were the one-dimensional line spanned
by `(1,1,1)`, then `tau_c=3(1,1,1)`, forcing

```text
(v_-,v_+)=(-9/2,9/2).                                   (13)
```

But then `tau_e=(-9,4u,9)`, whose first and third entries differ, a
contradiction.

Taking the determinant of `(1,tau_c,tau_e)` yields the exact necessary slope
equation

```text
4u(v_-+v_+)+4v_-v_+-27v_-+27v_+-162=0.                 (14)
```

This is a necessary retained-jet condition, not a sufficient condition for
`(11)` in the global ring.

The sharp hostile quadruple from THM-3634 lies only in
`C+x(x^2-1)C[x]`.  Its three derivative columns `dU,dV,dA` have determinant

```text
-225/2.                                                  (15)
```

Hence that hostile does not even lift through the actual one-jet image for
any `Q`; its advertised actual-ring caveat is essential.

## 2. Ordinary retained triples and the completed local image

Regard the three rows `(tau_c,i,tau_e,i)` as tangent directions in the
`(c,e)` plane.  If all three pairwise determinants are nonzero, the retained
curve image has an ordinary triple point at `(0,0,-3)`.  With branch
parameters `xi_-,xi_0,xi_+`, its completed normalization is

```text
C[[xi_-]] direct_product C[[xi_0]] direct_product C[[xi_+]], (16)
```

and its conductor in the normalization is

```text
xi_-^2 C[[xi_-]] direct_product
xi_0^2 C[[xi_0]] direct_product
xi_+^2 C[[xi_+]].                                      (17)
```

Equivalently, the completed local image consists exactly of triples of
series whose constant terms agree and whose linear-coefficient vector lies
in `T_Q`.  There is no further retained **formal-local** jet condition on
this ordinary stratum.  This does not remove global algebraization debt.

For `rho=0,u=1`, the six Hermite conditions are

```text
Q(-1)=-3, Q(0)=-3/4, Q(1)=-3,
Q'(-1)=-9/2, Q'(0)=1, Q'(1)=9/2.                        (18)
```

Their unique polynomial of degree at most five is

```text
Q_1=x^5+(9/2)x^4-2x^3-(27/4)x^2+x-3/4
   =(4x^5+18x^4-8x^3-27x^2+4x-3)/4.                   (19)
```

Its leading coefficient is nonzero, so degree five is minimal.  Here

```text
tau_c=(3,3,3),             tau_e=(-9,4,9),
Delta_(-,0)=39,            Delta_(0,+)=15,
Delta_(-,+)=54.                                             (20)
```

Thus the minimal Hermite curve is an ordinary retained triple.

## 3. Exact global structure of the minimal actual ring

From now on fix `Q=Q_1` and abbreviate `S=S_(Q_1)`.  The restricted
polynomials `z,c,b` are monic of degrees `12,15,21`.  Define the further
monic elements

```text
g_35=(2/9)(z^3-cb),
g_44=(2/9)(z^2b-c^3),
g_58=(4/81)(z^5-c^4-9c g_44).                          (21)
```

Their subscripts are their exact `x`-degrees.  The six elements of degrees

```text
12,15,21,35,44,58                                       (22)
```

form a finite SAGBI basis for `S`.  A particularly compact certificate is
the following `C[z]`-module presentation.  In residue order modulo `12`, put

```text
(a_0,...,a_11)=(0,73,50,15,88,65,30,79,44,21,58,35),  (23)

(p_0,...,p_11)=
(1, g_58 c, g_35 c, c, g_44^2, g_44 b,
 c^2, g_44 g_35, g_44, b, g_58, g_35).                 (24)
```

Then `deg p_r=a_r`, and exact reduction of the `36` products
`z p_r,c p_r,b p_r` gives

```text
S = direct_sum_(r=0)^11 C[z] p_r.                       (25)
```

Consequently the degree semigroup is generated by `(22)`, its Apéry set
modulo `12` is `(23)`, and its `41` gaps are

```text
1,2,3,4,5,6,7,8,9,10,11,13,14,16,17,18,19,20,
22,23,25,26,28,29,31,32,34,37,38,40,41,43,46,49,
52,53,55,61,64,67,76.                                  (26)
```

The semigroup conductor is `77`, and

```text
dim_C(C[x]/S)=41.                                       (27)
```

The fraction field of `S` is `C(x)`: the pole orders `12` and `35` of
`z,g_35` are coprime, so the polynomial Luroth generator has degree one.
Since `z` is monic of degree `12`, `x` is integral over `S`.  Hence

```text
normalization(S)=C[x].                                  (28)
```

## 4. Conductor and the finite quotient

Put `L=x(x^2-1)`.  The global conductor is

```text
(S:C[x])=h C[x],
h=L^2 H_76,                                             (29)
```

where `h,H_76` are monic, `deg h=82`, `deg H_76=76`, and `H_76` is
squarefree and coprime to `L`.  The factor `L^2` is the global trace of the
local order-two
conductor `(17)` at the retained triple; `H_76` records the remaining
normalization support and is not needed for the local theorem.

The conductor assertion has a finite exact certificate.  Because `z` is
monic of degree `12`, `C[x]` is generated over `S` by
`1,x,...,x^11`.  In the basis `(25)`, the conditions

```text
f x^j in S,                  1<=j<=11,                  (30)
```

have rank `41` on `S_(<=81)` (dimension `41`), and rank `41` on
`S_(<=82)` (dimension `42`).  Thus there is no nonzero conductor element
through degree `81`, while the degree-`82` kernel is one-dimensional and
has monic generator `h`.  Exact division gives `(29)`.

For compact byte-independent pinning, serialize a rational polynomial by
its complete descending coefficient list, including zeros, with every
coefficient in reduced `p/q` form joined by semicolons.  Then

```text
sha256(h)    =23956683350b2a13a5b07b99ec75ac26840a9933b0163e1696c687cdef1124f2,
sha256(H_76) =08b3d82e2ae493e8ffa079c2913a6c2095ac24c1157e8c78c9b2fb2476b4e3b8.
                                                                    (31)
```

## 5. A genuine actual-ring determinant witness

Take the fixed horizontal pair

```text
U=c,                         V=e.                       (32)
```

The ordinary polynomial derivatives `c',e'` are coprime.  Let `(A_0,B_0)`
be the unique reduced Euclidean pair satisfying

```text
c'B_0-e'A_0=1,       deg A_0<deg c'=14,
                     deg B_0<deg e'=11.                (33)
```

Every polynomial solution is

```text
A=A_0+c'T,                  B=B_0+e'T.                  (34)
```

Let `N=C[x]/S`, a `41`-dimensional vector space by `(27)`.  Membership of
both coefficients in the actual ring is the affine system

```text
[c'T]=-[A_0],                 [e'T]=-[B_0] in N.        (35)
```

On `C[x]_(<82)` with ordered basis `1,x,...,x^81`, its `82 x 82`
coefficient matrix and augmented matrix both have exact rank `81`.  The
unique free column is `x^81`.  The homogeneous kernel is generated by the
monic polynomial

```text
K=(1/27)L(27x^2-4x-36)H_76,             deg K=81.      (36)
```

Setting the free coefficient to zero gives a unique `T_*` of degree at most
`80`.  In fact

```text
deg T_*=80,
deg A=94,                    deg B=91,                  (37)

A(-1)=A(0)=A(1)=0,
B(-1)=B(0)=B(1)=1/3,                                (38)

c'B-Ae'=1.                                             (39)
```

The exact coefficient hashes are

```text
sha256(T_*)=b93aa43aec99bfbc94fc7fefd2dd5c57a815d5eb7e1bca1fc7ee2f48fd3fdd66,
sha256(A)  =c7253d6126f4acd03e38437b630e01d7dc6daeee64e92bf142dec70457ff00b5,
sha256(B)  =1b330351288f562690ee99c04e8a71d07c2ba0cfffd33f4980041302adda70d1.
                                                                    (40)
```

Thus `(U,V)` are independent modulo constants and `(A,B)` are independent
in `S`; `(32),(37)--(39)` give a genuine rank-two/rank-two actual-ring
point in the THM-3634 determinant-shadow cell.

## 6. Exhaustive fixed-pair degree boundary

For `d>=0`, define the complete filtered piece

```text
S_(<=d)={f in S: deg_x f<=d}                            (41)
```

and the exact linear map

```text
Phi_d:S_(<=d) direct_sum S_(<=d) -> C[x],
Phi_d(A,B)=c'B-Ae'.                                    (42)
```

This is the whole actual ring in the stated filtration, not a bounded list
of target monomials.  The module basis `(25)` gives

```text
d=93: dim S_(<=d)=53, matrix 108 x 106,
      rank Phi_d=106, rank[Phi_d|1]=107;                (43)

d=94: dim S_(<=d)=54, matrix 109 x 108,
      rank Phi_d=108, rank[Phi_d|1]=108.                (44)
```

Therefore no solution of `(39)` has

```text
max(deg A,deg B)<=93,                                   (45)
```

and the solution in `(37)` is the unique solution under the symmetric
cutoff `max(deg A,deg B)<=94`.

The quantifiers in this minimality statement are strict:

- `Q` is the one polynomial `Q_1` in `(19)`;
- `U,V` are fixed to `c,e`;
- `A,B` range over the complete actual ring `S`, filtered by source
  `x`-degree;
- the cutoff is symmetric in `A,B`.

No minimum is claimed over arbitrary horizontal pairs `U,V`, asymmetric
degree budgets, other Hermite polynomials, target degree, or all-order
quadratic-fold pairs.

## 7. Stopping boundary

The proved payload is:

```text
formal-local: ordinary retained triples have no jet debt beyond
              common value plus the derivative plane;

global t=0:   the minimal actual restriction ring contains an exact
              rank-two/rank-two determinant witness;

global fold:  OPEN -- no F^#,G^# in R[w] with constant pullback Jacobian
              for all t is constructed or excluded here;

Keller/JC(2): OPEN -- no polynomial automorphism or Keller conclusion
              follows from the first-stable shadow.                  (46)
```

The distinction between `(17)` and `(29)` is load-bearing: the former is a
formal-local conductor at the retained ordinary triple, while the latter is
the global conductor of one fixed polynomial curve and carries additional
normalization support.
