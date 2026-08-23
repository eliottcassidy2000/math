---
id: THM-3881
title: "Cusp-ideal residual transport and its rank-two matrix factorization"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  A
  universal two-sidecar identity transports the discriminant residual under
  every additive square/cube lift.  At the THM-3869 three-cusp square point,
  every polynomial-coefficient addition from the full cusp-value-zero ideal,
  together with every Delta-multiple ambiguity in its mixed lift, contracts
  exactly to a rank-two polynomial pair (T,f).  Its coefficient matrix has
  determinant Delta, the two sidecars become L^2 f and P f^2-T^2, and the
  remaining square question is one explicit quartic norm equation.  The
  arbitrary-degree equation, a Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3872 polynomial-coefficient lane, 2026-08-23
audit: >
  PROVISIONAL EXACT PROOF CANDIDATE.  The companion verifies the universal
  transport identity, the norm form for Delta, the exact cusp ideal and both
  Hilbert--Burch syzygies, the rank-two contraction and its two-sided matrix
  factorization, both defect sidecars, the full transported residual, the
  Delta-lift gauge, and two hostile pure-gauge boundary specializations in 32
  active gates.  Normal and optimized runs must byte-match the frozen output.
  Independent audit must recheck syzygy completeness, the converse pair
  parametrization, the every-lift quantifier, and the Mason--Stothers arm
  argument.
depends_on:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
  - THM-3869-three-cusp-square-residual-cardano-line-ramification
  - THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate
related:
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
script: 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
output: 05-knowledge/results/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.out
script_sha256: 87a2d045db3cfd73b072fbcfe6c1826495d68d77d1037988216bca62f858c2c0
output_sha256: 480e99a5f8913c443e23fbf990614eec5029e01c000d9c1ddf08e94fbb18e3fb
semantic_sha256: 4f369d3ca2bb65fa68654fd8ae6ba2193ad6ea6067c5c823ab67a1a2be342ac2
hash_basis: raw LF bytes
---

# THM-3881 -- the omitted ideal is a rank-two norm problem

**PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.**
Work over an algebraically closed field `k` of characteristic zero.  Write
`D=k[x,y]`.  The theorem has two layers:

1. a universal residual-transport identity valid in every commutative ring;
2. an exact rank-two presentation of all polynomial-coefficient `J` additions
   left open by THM-3872.

The second layer does **not** solve the resulting arbitrary-degree square
equation.  Its point is to replace three unconstrained coefficient polynomials
and a hidden lift ambiguity by one explicit two-coordinate norm object.

## 1. Universal two-sidecar transport

Let `P,Q,A,r` lie in a commutative ring and put

```text
P_r=P+2A+r^2,
Q_r=Q+3rP+3rA+r^3.                                      (1)
```

Direct expansion gives the identity

```text
P_r^3-Q_r^2
=P^3-Q^2
 +2(3A+3P+r^2)(AP-rQ)
 +(8A+6P+3r^2)(A^2-r^2P).                               (2)
```

Consequently, if a non-zero-divisor `Delta` divides all three defects, define

```text
S_0=(P^3-Q^2)/Delta,
C=(AP-rQ)/Delta,
B=(A^2-r^2P)/Delta.                                     (3)
```

Then the transported residual is exactly

```text
S_r=S_0+2(3A+3P+r^2)C+(8A+6P+3r^2)B.                   (4)
```

Thus only two sidecars are needed.  Formula `(2)` is a polynomial identity;
it loses no quotient information and needs no geometric hypothesis.

## 2. The three-cusp square point is already a norm form

Retain the THM-3854 branch polynomial `Delta` and the THM-3869 positive
representative `h_*`.  Set

```text
a=x+1,                   L=9x+4,
K=y^2-15x^2-15x-4,
P=aL^2,                  Q=KL^2.                         (5)
```

Then, identically in `D`,

```text
Delta=a^3L^2-K^2,
P^3-Q^2=Delta L^4.                                      (6)
```

The pullbacks of `P,Q` are `h_*^2,h_*^3`.  The base residual `L^4` is the
square from THM-3869; that theorem proves that its Cardano line `L=0` still
ramifies in the cubic field.

## 3. Exact presentation of the omitted cusp ideal

The ideal of the three cusp image points is

```text
J=(j_1,j_2,j_3)=(ax,ay,y^2+4x).                          (7)
```

The minimum mixed lifts of `j_i h_*` from THM-3872 have the compact forms

```text
A_1=xK,
A_2=yK,
A_3=(15x+4)K+aP.                                        (8)
```

There are two generating syzygies

```text
s_1=(-y,x,0),                 s_2=(4,y,-a).              (9)
```

They act differently on the mixed lifts:

```text
s_1 dot (j_1,j_2,j_3)=0,      s_1 dot (A_1,A_2,A_3)=0,
s_2 dot (j_1,j_2,j_3)=0,      s_2 dot (A_1,A_2,A_3)=-Delta. (10)
```

These are **all** syzygies.  Indeed, reduce a relation
`g_1j_1+g_2j_2+g_3j_3=0` modulo `a`.  Since
`j_3 mod a=y^2-4` is nonzero in `k[y]`, one has `a|g_3`.  After dividing by
`a` and using `j_3=y^2+4x`, the remaining relation is a syzygy of `(x,y)`,
hence a multiple of `(-y,x)`.  This gives exactly `(9)`.  Equivalently, the
signed maximal minors of

```text
[ -y   4 ]
[  x   y ]
[  0  -a ]                                                   (11)
```

are `(-j_1,j_2,-j_3)`.

The second line of `(10)` is the missing structural point: polynomial
coefficient changes in `J` and arbitrary `Delta`-multiple changes of the
mixed lift are not separate directions.  The syzygy `s_2` identifies them.

## 4. Contraction to the pair `(T,f)`

For arbitrary `g_1,g_2,g_3 in D`, put

```text
f=g_3,
T=xg_1+yg_2+(15x+4)g_3,
r=g_1j_1+g_2j_2+g_3j_3,
A=g_1A_1+g_2A_2+g_3A_3.                                (12)
```

Equations `(5),(7),(8)` give

```text
r=aT+Kf,
A=KT+aPf.                                                (13)
```

The pair always satisfies one address condition,

```text
T(0,0)=4f(0,0).                                         (14)
```

Conversely, `(14)` is sufficient.  It says that
`T-(15x+4)f` belongs to `(x,y)`, so write that difference as
`xg_1+yg_2` and recover `(12)`.  Therefore the full coefficient universe is
exactly

```text
E={(T,f) in D^2 : T(0,0)=4f(0,0)}.                       (15)
```

In matrix form, `(13)` is

```text
[r]   [a   K ][T]
[A] = [K  aP ][f].                                      (16)
```

The matrix and its adjugate form the exact factorization

```text
[a   K ][aP -K]             [aP -K][a   K]
[K  aP ][-K  a] =Delta I_2= [-K  a][K  aP].             (17)
```

In particular its determinant is `Delta`, so `(T,f)` is uniquely determined
by `(r,A)` in the domain `D`.

This also proves the promised every-lift quantifier.  For a fixed `r in J`,
any two polynomial lifts of `rh_*` differ by a multiple of `Delta`, because
THM-3854 identifies the kernel of `D -> k[t]` with `(Delta)`.  Adding a
multiple of `s_2` realizes every such change by `(10)`.  On pairs the action is

```text
(T,f) -> (T+Kq,f-aq),
r      -> r,
A      -> A-Delta q.                                    (18)
```

Thus `(15)-(18)` include both arbitrary polynomial coefficients and every
mixed-lift gauge, without duplication or loss.

## 5. The two sidecars are a coordinate and a binary norm

Substitute `(5),(13)` into `(3)`.  The cross terms cancel exactly:

```text
AP-rQ=Delta L^2f,
A^2-r^2P=Delta(Pf^2-T^2).                               (19)
```

Hence

```text
C=L^2f,                     B=Pf^2-T^2.                  (20)
```

The full residual left by THM-3872 is therefore

```text
S(T,f)=L^4
 +2(3A+3P+r^2)L^2f
 +(8A+6P+3r^2)(Pf^2-T^2),                               (21)
```

where `r,A` are the two linear forms `(13)` and `(T,f)` ranges over `(15)`.
This is a quartic polynomial in the two module coordinates.  Since `D` is a
UFD and `k` is algebraically closed, a polynomial square in `Frac(D)` is a
polynomial square up to a square unit.  A translated square/cube pair has
square residual **if and only if**

```text
S(T,f)=G^2                                                  (22)
```

for some `G in D`.  This is the exact surviving polynomial problem.

The constant span closed by THM-3872 is the small slice

```text
f=w,                  T=ux+vy+w(15x+4).                   (23)
```

Thus THM-3872 is a genuine bounded cell of `(22)`, while `(21)` exposes the
previously hidden higher-degree directions.

## 6. A first global obstruction inside the lift gauge

Take the pure `s_2` gauge through the base point:

```text
(T,f)=(Kq,-aq),             r=0,       A=-Delta q.        (24)
```

Its residual is

```text
S_q=L^4-6P^2q+12Delta Pq^2-8Delta^2q^3.                  (25)
```

If `S_q` is a square in `D`, then necessarily

```text
a divides q.                                                (26)
```

To prove this, specialize `a=0`, equivalently `x=-1`.  With
`q_0(y)=q(-1,y)`, equation `(25)` becomes

```text
G(y)^2+8(y^2-4)^4q_0(y)^3=625.                           (27)
```

Suppose `q_0` is nonzero and put `n=deg q_0`.  The two nonconstant summands
in `(27)` are coprime, since their sum is the nonzero constant `625`.  Their
common degree is `8+3n=2deg G`; if `8+3n` is odd this is already impossible.
In the remaining case Mason--Stothers gives

```text
8+3n <= deg rad(G(y)(y^2-4)q_0(y))-1
       <= (8+3n)/2+n+1,                                  (28)
```

which says `3+n/2<=0`, impossible.  Hence `q_0=0`, proving `(26)`.

There is an independent second boundary condition.  At `L=0`, or
`x=-4/9`, one obtains

```text
S_q=-8(y^2-8/27)^4 q(-4/9,y)^3.                          (29)
```

Since `k[y]` is a UFD and `k` is algebraically closed, a square survivor must
have `q(-4/9,y)=0` or `q(-4/9,y)` itself a square.  Conditions `(26),(29)` do
not yet close the case `q=aq_1`; that branch remains part of `(22)`.

## 7. Exact scope

This theorem proves the universal identity, the complete presentation of all
coefficient and lift choices, the matrix factorization, the norm sidecars,
and the pure-gauge arm obstruction.  It does **not** prove that `(22)` has no
higher-degree solutions.  In particular, pairs with nonconstant `f`, pairs
with `q=aq_1`, alternate square/cube lifts outside the displayed additive
grammar, a polynomial-plane Keller atlas, and JC(2) remain **OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
python3 -O 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
```

Both runs must byte-match
`05-knowledge/results/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.out`.
