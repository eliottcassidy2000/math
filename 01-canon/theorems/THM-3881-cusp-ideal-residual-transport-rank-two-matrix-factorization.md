---
id: THM-3881
title: "Cusp-ideal residual transport and its rank-two matrix factorization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A universal two-sidecar identity
  transports the discriminant residual under every additive square/cube lift.
  At the THM-3869 three-cusp square point, every polynomial-coefficient
  addition from the full cusp-value-zero ideal, together with every
  Delta-multiple ambiguity in its mixed lift, contracts exactly to a rank-two
  polynomial pair (T,f).  Its coefficient matrix has determinant Delta and
  its sidecars are L^2 f and P f^2-T^2.  The full theorem proves that
  the full T=0 lane and every affine-slope lane T=hf with deg h<=1 have
  the unique square survivor (T,f)=(0,0).  The general arbitrary-degree
  equation, quadratic and higher slopes, a Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3872 polynomial-coefficient lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS for Sections 1--6 and the original scope
  (jc_quartic_c3_construct, 2026-08-23): syzygy completeness/signs, every-lift
  realization, converse address parametrization, determinant uniqueness,
  both sidecars, transport, Mason arm, L square class, normal/-O replay, and
  frozen hashes all passed.  A focused follow-up audit rederived Section 7's
  two factorizations, checked every constant/nonconstant edge without using
  the address condition, and independently proved the affine-slope extension
  `T=hf`, `deg h<=1`.  Its 52-gate audit and the 223-gate primary companion
  pass normally and under `-O`, byte-match their frozen outputs, and isolate
  characteristic three and `deg h=2` as proof boundaries rather than
  counterexamples.
depends_on:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
  - THM-3869-three-cusp-square-residual-cardano-line-ramification
  - THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate
related:
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
script: 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
output: 05-knowledge/results/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.out
script_sha256: 0445ca182f8bf809b352f6ba31fd91ae1bfc84118e168cc44d7de1cf64e880ad
output_sha256: b486419fc5047e279d56b89e8b22d56bb8840c641d5bc13a324358ff02071757
semantic_sha256: 7e5ecbb745d9ac184eb602dbcd09bc84cc44847926d887148ace4984e34ceb87
affine_script: 04-computation/jc2_cusp_residual_affine_slope_thm3881.py
affine_output: 05-knowledge/results/jc2_cusp_residual_affine_slope_thm3881.out
affine_script_sha256: 7966a6b010a85e2941e5ad7a36a54fcf2c7b3fcc78d69257fe972b636d80790e
affine_output_sha256: 116dbf17d57098d32e88c49947d89a157f08af3c1be264fa5b421afa82f2f420
affine_semantic_sha256: dcfd5de0a7eab5e173868bec8de0ec10b9557a2794ad2cded04d6a069a5d12c7
affine_independent_script: 04-computation/jc2_cusp_residual_affine_slope_independent_audit_thm3881.py
affine_independent_output: 05-knowledge/results/jc2_cusp_residual_affine_slope_independent_audit_thm3881.out
affine_independent_script_sha256: bd8b5fb4e3de70f18e57fed84a3ee6f7ded77511cf989e850e33ebc835b6a182
affine_independent_output_sha256: daf1dc441798e761385c5cc791d0fc3703fc54f03b467b51b3d332b5e5c3cdcd
affine_independent_semantic_sha256: dac48c5176e15c79596f6319609d8bfbc94425a049d1fd2f51a9c49256fdce5c
hash_basis: raw LF bytes
---

# THM-3881 -- the omitted ideal is a rank-two norm problem

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
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

## 7. The full `T=0` lane has one square survivor

There is also an all-degree consequence not visible in the constant slice.
Set `T=0`; then `(13)` becomes

```text
r=Kf,                         A=aPf.                       (30)
```

Exact simplification of `(21)`, followed by `(6)`, gives two equivalent
forms

```text
S(0,f)=L^2 H_f,
H_f=L^2(1+af)^3(1+3af)-Delta f^3(2+3af)
   =L^2(1+2af)^3+K^2f^3(2+3af).                          (31)
```

Suppose `f` is nonconstant of total degree `n>=1`.  In the second form of
`(31)`, the term

```text
3aK^2f^4                                                     (32)
```

has degree exactly `4n+5`.  It is the unique term of that degree: the full
first summand has degree at most `3n+5`, and the remaining term `2K^2f^3`
has degree `3n+4`.  Since `n>=1`, both bounds are strictly below `4n+5`.
Its leading form is

```text
3x(y^2-15x^2)^2 f_n^4,
```

which is nonzero in the domain `D`.  Thus `deg H_f=4n+5` is odd, so `H_f`
and equivalently `S(0,f)=L^2H_f` cannot be squares.

The address condition is not needed at the constant edge.  If `f=c!=0`, the
degree-five homogeneous part of `H_c` is

```text
648c^3x^5+3c^4x(y^2-15x^2)^2.                             (33)
```

Its `xy^4` coefficient is `3c^4`, so it is nonzero and `deg H_c=5`.
Finally `f=0` gives `H_f=L^2` and `S(0,0)=L^4`.  Hence

```text
S(0,f) is a square  iff  f=0.                              (34)
```

This degree proof works over every integral-domain coefficient ring of
characteristic different from three.  No algebraic closure or UFD
factorization is needed.  Characteristic three is a genuine boundary of this
argument: the displayed top coefficients vanish, but no characteristic-three
square counterexample is claimed.

## 8. Every affine-slope lane `T=hf`, `deg h<=1`, closes

Let `h` be affine-linear and set

```text
T=hf,       R=ah+K,       M=Kh+aP,       N=P-h^2.          (35)
```

Then `r=Rf`, `A=Mf`, and exact expansion of `(21)` gives

```text
S(hf,f)=L^4+6PL^2f+6(ML^2+PN)f^2
        +2(R^2L^2+4MN)f^3+3R^2Nf^4.                       (36)
```

Write `h=h_1+h_0`, with `h_1` homogeneous linear, and put
`K_2=y^2-15x^2`.  The load-bearing leading forms are

```text
R_2=K_2+xh_1,          N_3=81x^3,          M_4=81x^4.     (37)
```

The coefficient rows of `f,f^2,f^3,f^4` in `(36)` have degrees at most
`5,6,7,7`.  Since the `y^2` coefficient of `R_2` is one, none of the
load-bearing forms can vanish.  If `deg f=n>=1`, the unique leading form is

```text
243x^3(K_2+xh_1)^2 f_n^4,                                 (38)
```

of odd degree `4n+7`; the lower rows have degrees at most
`3n+7,2n+6,n+5,4`.  If `f=c!=0`, the two tied degree-seven rows have
`x^3y^4` coefficient `243c^4`, so they cannot cancel.  Conversely `f=0`
forces `T=0` and gives `S=L^4`.  Therefore, without any address condition,

```text
deg h<=1 and T=hf:
S(T,f) is a square  iff  (T,f)=(0,0).                      (39)
```

This applies a fortiori to the address module `E`.  At `deg h=2` the `f^4`
coefficient has degree ten rather than seven, so this odd-degree mechanism
stops.  That is a proof-boundary witness, not a square solution.

## 9. Exact scope

This theorem proves the universal identity, the complete presentation of all
coefficient and lift choices, the matrix factorization, the norm sidecars,
the pure-gauge arm obstruction, the complete `T=0` closure, and every
affine-slope lane `T=hf` with `deg h<=1`.  It does **not** prove that `(22)`
has no higher-degree solutions.  In particular, arbitrary pairs with
`T!=0`, quadratic and higher slopes, pairs with `q=aq_1`, alternate
square/cube lifts outside the displayed additive grammar, a polynomial-plane
Keller atlas, and JC(2) remain **OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
python3 -O 04-computation/jc2_cusp_ideal_residual_transport_matrix_factorization_thm3881.py
python3 04-computation/jc2_cusp_residual_affine_slope_thm3881.py
python3 -O 04-computation/jc2_cusp_residual_affine_slope_thm3881.py
python3 04-computation/jc2_cusp_residual_affine_slope_independent_audit_thm3881.py
python3 -O 04-computation/jc2_cusp_residual_affine_slope_independent_audit_thm3881.py
```

Each normal/optimized pair must byte-match the corresponding frozen output in
`05-knowledge/results/`.
