---
id: THM-2942
title: "Macaulay extraneous flag factor and Pluecker mutation"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  Every six-quartic-row mutation of the fixed
  degree-seven 20Q+10C+6F Macaulay chart is the genuine ternary
  resultant times a Pluecker flag coordinate.  For the original chart
  the exact factor is q200^6*c300*K, and the two-row mutation has factor
  q200^5*c300*P_alt.  The 210 coordinates satisfy the Grassmann
  exchange laws.  On factorial four-slot forms, K is independent of the
  second interior slot.  In all 1,140 supports of width at most 20, the
  common factor of K and P_alt is positive on n>=0 and the reduced
  cofactors are coprime.  This separates chart walls from the resultant
  wall at the variable Pluecker-coordinate level; the common selected
  basis scale and the resultant wall remain.  No resultant
  nonvanishing or arbitrary-width SFC(4) is proved.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
related:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
script: 04-computation/gmc_macaulay_extraneous_flag_pluecker_thm2942.py
output: 05-knowledge/results/gmc_macaulay_extraneous_flag_pluecker_thm2942.out
script_sha256: 85ce40de2aa777b8af091dfec934de68beeb42676bea3b48dbf815990a51e0e9
output_sha256: 6819fcedf917ca3a09ff1a4413e86d7942189862420af5bb4d2154d182a3b014
constructor_dependency_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
denominator_dependency_sha256: 83d70a95f0943992d0e4b7027eede431d4dc968b66655e37b43fd0acfc692e47
hash_basis: LF-normalized bytes
---

# THM-2942 -- Macaulay extraneous flag factor and Pluecker mutation

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Universal factorization

Let

```text
S=k[x0,x1,x2],              char(k)=0,                 (1)
```

and let `Q,C,F` be arbitrary homogeneous forms of degrees `2,3,4`.
Write their coefficients as

```text
Q=sum q_ijk x0^i x1^j x2^k,
C=sum c_ijk x0^i x1^j x2^k.                           (2)
```

Use the degree-seven Macaulay row order of THM-2921.  Its fixed
`36`-row chart consists of

```text
20 rows from Q*S5,
10 rows from C*S4,
 6 rows from F*S3.                                    (3)
```

The original six degree-three multipliers, in the inherited order

```text
x2^3, x1*x2^2, x1^2*x2, x1^3,
x0*x2^2, x0*x1*x2, x0*x1^2, x0^2*x2, x0^2*x1, x0^3, (4)
```

have indices

```text
J0={0,1,2,3,4,5}.                                    (5)
```

If `Delta_J0(Q,C,F)` is the selected determinant, then

```text
Delta_J0
 =q200^6*c300*K(Q,C)*Res(Q,C,F),                      (6)
```

where the resultant has the standard normalization and

```text
K
 =c120*q200^2-c210*q110*q200
  -c300*q020*q200+c300*q110^2.                       (7)
```

Thus the determinant has two distinct kinds of wall:

```text
Res(Q,C,F)=0                 genuine common-zero wall;
q200*c300*K(Q,C)=0           selected-flag/chart wall. (8)
```

The distinction in `(8)` is load-bearing.  Vanishing of one selected
Macaulay determinant does not by itself exhibit a common projective
zero.

## 2. Six exact blocks

The determinant in `(6)` has coefficient multidegree

```text
(deg_Q,deg_C,deg_F)=(20,10,6),                        (9)
```

whereas

```text
deg Res(Q,C,F)=(12,8,6).                              (10)
```

If `Q,C,F` have a common projective zero, evaluation at that zero gives
a right kernel of every maximal Macaulay minor.  Irreducibility of the
universal resultant therefore gives

```text
Delta_J0=E_J0(Q,C)*Res(Q,C,F),                        (11)
```

with `E_J0` independent of `F` and of bidegree `(8,2)` in `(Q,C)`.

Set `F=x2^4`.  The six selected quartic rows become unit pivots.  After
deleting their rows and columns, order the remaining rows and columns
by `x2` degree.  Terms containing an additional `x2` strictly raise
that degree, so the residual matrix is block triangular.  Put

```text
q(t)=q200+q110*t+q020*t^2,
c(t)=c300+c210*t+c120*t^2+c030*t^3,
R=Res_t(q(t),c(t)).                                   (12)
```

The six diagonal blocks have sizes and determinants

| `x2` degree | size | determinant |
|---:|---:|---|
| `0` | `8` | `c300*q200^2*R` |
| `1` | `7` | `q200^2*R` |
| `2` | `6` | `q200*R` |
| `3` | `5` | `R` |
| `4` | `3` | `K` |
| `5` | `1` | `q200` |

Their product is

```text
q200^6*c300*K*R^4.                                   (13)
```

The standard specialization identity is

```text
Res(Q,C,x2^4)=R^4.                                   (14)
```

Since `E_J0` is independent of `F`, equations `(11)--(14)` prove `(6)`
on a dense open set and hence as a polynomial identity.

For the pure powers

```text
Q=(u0*x0+u1*x1+u2*x2)^2,
C=(v0*x0+v1*x1+v2*x2)^3,                             (15)
```

equation `(7)` becomes

```text
K=3*v0*u0^2*(u0*v1-u1*v0)^2.                        (16)
```

Together with

```text
Res(Q,C,F)=det[u;v;w]^24
```

when `F=(w dot x)^4`, this recovers exactly the flag factor in
THM-2927.

## 3. All 210 quartic charts are one Pluecker object

Let

```text
A=S/(Q,C).                                            (17)
```

On the dense selected flag subchart where `Q,C` are a regular sequence
and the first thirty rows in `(3)` are independent,

```text
dim A_3=dim A_7=6.                                   (18)
```

those thirty rows are a selected basis of the degree-seven ideal piece
`(Q,C)_7`.  Multiplication by `F` therefore descends to

```text
mu_F:A_3 -> A_7.                                     (19)
```

Its determinant is the norm of `F` on the `(2,3)` complete
intersection, hence the ternary resultant, up to the fixed basis
factor already calibrated in `(6)`.

The ten degree-three monomials in `(4)` map to ten vectors in the
six-dimensional space `A_3`.  For every six-subset

```text
J subset {0,...,9},                  |J|=6,            (20)
```

let `p_J(Q,C)` be their signed Pluecker coordinate.  Equivalently,
form the `4`-by-`10` relation matrix whose rows are

```text
x2*Q, x1*Q, x0*Q, C                                  (21)
```

in the monomial order `(4)`; `p_J` is the correspondingly signed
complementary `4`-minor.

With the complementary-minor sign in `(21)`, exterior multiplication
in `(19)` gives

```text
Delta_J
 =-q200^5*c300*p_J(Q,C)*Res(Q,C,F).                  (22)
```

The original coordinate is

```text
p_J0=-q200*K.                                        (23)
```

Thus `(22)` specializes to `(6)`.  In particular all `210` selected
Macaulay determinants obey every Grassmann-Pluecker exchange relation
after removal of their common resultant and common flag scale.

This gives a lawful chart atlas: changing quartic rows changes the
Pluecker coordinate, not the common-zero divisor.

## 4. The two-row mutation

Replace the quartic multiplier indices `1,2` by `6,7`:

```text
J1={0,3,4,5,6,7}.                                    (24)
```

Its complementary minor is

```text
p_J1=-P_alt,                                         (25)
```

where

```text
P_alt
 = c012*q020*q200^2
  -c021*q011*q200^2
  -c210*q002*q020*q200
  +c210*q011^2*q200
  +c300*q002*q020*q110
  -c300*q011^2*q110.                                 (26)
```

Therefore

```text
Delta_J1
 =q200^5*c300*P_alt*Res(Q,C,F).                      (27)
```

Equations `(6)` and `(27)` show exactly what a two-chart proof must
establish.  At a point where the common scalar `q200*c300` is nonzero,
simultaneous vanishing of the two charts can arise either from the
genuine resultant wall or from the simultaneous flag wall

```text
K=P_alt=0.                                           (28)
```

Separating the flag wall does not separate the resultant wall.

## 5. A one-cut positivity operator

The following independent lemma is useful for chart cofactors.

Let

```text
p(x)=sum_(j=0)^d a_j x^j,             B>0,            (29)
```

and suppose the coefficients have one sign cut: for some `c`,

```text
a_j<=0  for j<c,               a_j>=0 for j>=c.       (30)
```

If `p(B)>0`, then

```text
p^(m)(x)>0                  (x>=B, 1<=m<=d),           (31)
Delta^m p(B)>0             (1<=m<=d).                 (32)
```

Indeed,

```text
B^m p^(m)(B)=sum_j (j)_m a_j B^j,                    (33)
```

and the weights `(j)_m` are nondecreasing across the sign cut.  Comparing
them with the cut weight and using the positive top coefficient gives
strict positivity.  Taylor expansion at `B` proves `(31)`, and

```text
Delta^m p(B)
 =int_[0,1]^m p^(m)(B+t1+...+tm) dt1...dtm           (34)
```

proves `(32)`.

If instead `B` is a positive integral root and the hypotheses
`(29)--(30)` hold, the constant Gregory-Newton coefficient vanishes
but every coefficient of order `1,...,d` is positive.  All higher
orders vanish for the degree reason.  This is a conditional positivity
gate; the theorem does not assert that every factorial flag cofactor
has one sign cut.

## 6. Factorial four-slot consequences

Return to the factorial functional `L(s^j)=j!` and support

```text
(n,n+a,n+b,n+M),             0<a<b<M.                 (35)
```

Use THM-2921's endpoint-elimination coordinates `(x0,x1,x2)`.  The
coefficients appearing in `K` are precisely those of

```text
Q(x0,x1,0),                         C(x0,x1,0).        (36)
```

They depend on `(M,a,n)` but not on the second interior offset `b`.
Consequently

```text
K_(M,a,b)(n)=K_(M,a)(n).                              (37)
```

This is the exact reason the first chart can have a stable
first-gap-controlled sign transition while the full resultant still
depends on both interior slots.

The companion also gives the following finite exact separation result.
For every one of the

```text
1,140 supports with 3<=M<=20 and 0<a<b<M,             (38)
```

let `G=gcd(K,P_alt)` in `Z[n]`, normalized to have positive leading
coefficient.  Then

```text
G(0)>0,
every coefficient of G is nonnegative,
gcd(K/G,P_alt/G)=1.                                  (39)
```

Hence `G(n)>0` for every `n>=0`, and the two reduced flag cofactors
have no common complex root.  In particular the two flag coordinates
cannot vanish simultaneously at a nonnegative integral depth in this
finite support bank: factorial positivity gives

```text
q200(n)>0                                                   (39a)
```

after the positive THM-2925 scaling, so `p_J0=-q200*K` and
`p_J1=-P_alt` have a common zero only if `K,P_alt` do.

This is a statement about the two variable Pluecker coordinates.  The
common selected-basis scale `c300` in `(22)` may still vanish, and the
resultant may still vanish.

A seductive closed formula for `deg G` is false.  The unique hostile
inside `(38)` is

```text
(M,a,b)=(11,7,8):
G=16*n^2+320*n+1584=16*(n+9)*(n+11),                 (40)
```

of degree `2`, whereas the naive formula predicts `1`.  No degree law
is claimed.

## 7. What this does and does not prove

The theorem proves three reusable facts:

1. the fixed determinant is the resultant times an explicit flag
   coordinate;
2. all quartic-row changes form one exact Pluecker atlas; and
3. the original/two-row-mutated variable Pluecker walls are separated
   throughout the finite bank `(38)`.

It does **not** prove that `Res(Q,C,F)` is nonzero.  Therefore it does
not by itself close any new width, prove arbitrary-width SFC(4), or
prove an arbitrary-radial GMC statement.  A finite-width application
must still show that at least one full determinant is nonzero, not
merely that at least one variable Pluecker coordinate is nonzero; the
common selected-basis scale must also be retained.

The quotient interpretation is also chart-relative.  On the boundary
where the first thirty rows cease to be a basis of `(Q,C)_7`, the
displayed quotient bases do not apply.  The identities `(6)`, `(22)`
and `(27)`, first proved on the dense selected flag subchart, remain
global polynomial identities.  Repeated support offsets likewise
collapse the response flags and are outside `(35)`.

## 8. Exact verification

The companion:

- derives all six symbolic blocks in Section 2 and the binary
  resultant specialization;
- checks the original and mutated complementary minors symbolically;
- checks all `210` quartic-row charts over `F_1000003`, finding `207`
  nonzero sample coordinates with one common scale and verifying a
  three-term Pluecker exchange on both sides;
- checks the pure-power specialization `(16)`;
- checks `649` positive-base and `78` positive-root instances of the
  one-cut derivative/difference lemma;
- checks the first-gap independence `(37)` in `480` exact cases; and
- checks `(39)` on all `1,140` supports in `(38)`, including the unique
  hostile `(40)`.

The finite-field census audits the exterior-algebra identity; the
characteristic-zero proof is Sections 1--4.  The finite factorial scan
is explicitly limited to `(38)`.

Run

```text
python 04-computation/gmc_macaulay_extraneous_flag_pluecker_thm2942.py
python -O 04-computation/gmc_macaulay_extraneous_flag_pluecker_thm2942.py
```

Normal and optimized output must byte-match the stored transcript and
the LF-normalized hashes above.

**QED, pending independent hostile audit.**
