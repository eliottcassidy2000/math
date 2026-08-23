---
id: THM-3821
title: "The first RZ2 extension on the cubic pseudo-plane has no Darboux pair"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The complete
  first rz^2 extension of the THM-3814 canonical nodal profile admits no
  Darboux pair.  A hypothetical
  pair must enter one of two exact Kummer anatomies: the generic odd ladder
  7,5,3,1 with a Riccati square payment, or the P=0 skip ladder 7,3,1 with a
  linear payment.  The next r^2z bucket excludes every nonzero Kummer root;
  the remaining confluent monomial towers fail by the origin jets and the
  quadratic-arm root ratio.  At the intermediate top step S=0 is impossible,
  while T=0 remains allowed through mu=0 before the later contradiction.
  This is a complete no-go for the stated ansatz, not for arbitrary profiles
  on the surface and not a planar-JC theorem.
source: jc_zero_debt_lift / cubic-pseudoplane rz2 odd-ladder lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDITS PASS AFTER REPAIR (root / jc-cohn-boundary and
  root post-promotion, 2026-08-23).  The audits rederived the Poisson signs,
  monic normal form, every Kummer valuation family, both terminal laws, and
  the constant, linear, nonzero-root, and monomial closures.  A post-promotion
  audit found one endpoint seam: the claimed exact
  origin valuation for R fails when ord_0(v)=0 because the Pg block can tie
  the forcing.  The repaired proof uses the full r^2z origin bucket to force
  R(0)=0; exact equality remains valid when ord_0(v)>0.  It also repaired the
  one-sided wording because T=mu*S permits mu=0, and retained arbitrary g(e)
  through the confluent-root evaluation.  Independently, the all-degree audit
  expanded the full degenerate r^2z polynomial, obtained the impossible
  half-integral resonance m=5/2, and checked the generic root payment, origin
  jets, odd opposite-root seam, and root-of-unity contradiction.  The repaired
  deterministic companion checks all of these in 72 active gates without
  inactive Python asserts; normal and optimized runs byte-match the frozen
  output and raw hashes recorded below.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3814-nodal-rz-kummer-profile-degree-gate
related:
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_rz2_odd_ladder_thm3821.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_rz2_odd_ladder_thm3821.out
script_sha256: a4567cd4bdd37df9a2db260d9c00e1275211ef1735690591a7511f458a269e18
output_sha256: 22eaf1a707edb91974c83c9f2c11a53f551468a1d3b805c2715ae79558910852
semantic_sha256: f73ecf90a4c8fc5a77ce7471c54f8a6eccc0a4c542d7b95ca7fcd5b63fc3c4d7
hash_basis: raw LF bytes
---

# THM-3821 -- the first RZ2 layer is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be an
algebraically closed field of characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles `f,g,h,kappa,p,q,S,T in k[e]`, define

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T.              (2)
```

There are no such profiles for which

```text
{A,C}=1.                                               (3)
```

In particular, this complete first `rz^2` extension of the THM-3814 profile
contains no Darboux pair.

Assume the contrary.  Sections 1--5 force `mu in k`, nonzero `alpha,beta in
k`, and a polynomial `v in k[e]` such that, on writing

```text
K=kappa-mu f,
P=q-mu p,
Q=h-mu g,
D=3e^2-2mu e-1,                                       (4)
```

one has

```text
T=mu S,
S=alpha e^4v^7,
K=beta e^2v^4,
D f=2eK-1/12.                                         (5)
```

Exactly one of the following intermediate branches occurs.

### Generic branch: `P!=0`

There are `delta in k*`, `eta in k`, and a nonzero `U in k[e]` such that

```text
P=delta e^3v^5,
p=e^2v^3U,
Q=(5delta/(7alpha)) e v U+(eta alpha/beta)e^2v^3.      (6)
```

The profiles satisfy the terminal Riccati law `(32)` below.  In particular,
at every nonzero root `rho` of `v`,

```text
28alpha beta f(rho)=5delta U(rho)^2.                   (7)
```

### Degenerate branch: `P=0`

There are `theta in k*` and a nonzero `U in k[e]` such that

```text
Q=theta e^2v^3,
p=e v U.                                               (8)
```

The profiles satisfy the terminal linear-payment law `(38)` below, and at
every nonzero root `rho` of `v`,

```text
4beta f(rho)+3theta U(rho)=0.                          (9)
```

Sections 6--8 exclude both branches for every polynomial `v`, completing the
contradiction.

## 1. The descending canonical buckets

As in THM-3814, the monic presentation

```text
B=k[r,e][z]/(z^3-r^2e-r)                              (10)
```

makes `(1,z,z^2)` a free basis and gives a unique normal form.  The arm
coefficient is unchanged:

```text
[z]({A,C}-1)=(36e^2f-24e kappa-12f+1)/2.              (11)
```

Thus `(3)` requires

```text
(3e^2-1)f-2e kappa=-1/12,              f(0)=1/12.     (12)
```

The next three descending coefficients are

```text
[r^5]=-21e^2(-ST'+TS'),                               (13)

[r^4]=-3e(-4efT'+4e kappa S'-7eS kappa'
           +7eTf'+2fT-2kappa S-12ST'+12TS'),          (14)

[r^3z^2]=-3(-5epT'+5eqS'-7eSq'+7eTp'-pT+qS).          (15)
```

The exact companion also freezes the longer `r^3z` and `r^3` buckets from
which Sections 3--5 descend.

## 2. The new `S` profile is necessarily nonzero

Suppose `S=0`.  If `T=0`, equation `(2)` is the complete profile excluded
by THM-3814.  If `T!=0`, equation `(14)` becomes

```text
4efT'-7eTf'-2Tf=0.                                    (16)
```

In `k(e)`, this says `(T^4/(e^2f^7))'=0`.  Unique factorization therefore
gives

```text
T=alpha e^4w^7,
f=beta e^2w^4                                         (17)
```

with nonzero constants.  This contradicts `f(0)=1/12`.  Hence `S!=0`.

Equation `(13)` and the characteristic-zero constant-field fact now give
`T=mu S`.  With `K` as in `(4)`, equation `(12)` becomes the last identity
in `(5)`.  This relation does not force `mu!=0`; the boundary `T=0` remains
available when `mu=0`.  The polynomial `K` cannot vanish: otherwise the
nonconstant quadratic `D` times `f` would be `-1/12`.

Substitution in `(14)` gives

```text
4eK S'-7eS K'-2KS=0,
(S^4/(e^2K^7))'=0.                                    (18)
```

At an irreducible factor other than `e`, the valuation equation is
`4a=7b`; at `e` it is `4a=2+7b`.  Their nonnegative solutions are

```text
(a,b)=(7j,4j),
(a,b)=(4+7j,2+4j),                 j>=0,              (19)
```

which proves the first two towers in `(5)`.

## 3. The generic branch integrates to 7-5-3-1

With `P=q-mu p`, equation `(15)` becomes

```text
7eS P'-5eP S'-PS=0.                                   (20)
```

If `P!=0`, then `(P^7/(eS^5))'=0`.  Substituting `(5)` and applying the
same UFD argument gives

```text
P=delta e^3v^5                         (delta in k*).  (21)
```

The `r^3z` bucket, written using `Q=h-mu g`, is

```text
5e(pP'-Pp')+7eS Q'-3eQ S'-2QS=0.                     (22)
```

Put `u=p/P` and `w=Q/S` temporarily in `k(e)`.  Equation `(18)` transforms
`(22)` into

```text
(Kw)'=(5P^2K/(7S^2))u'.                               (23)
```

The coefficient `P^2K/S^2` is constant by `(5),(21)`.  Integration gives

```text
Q=5Pp/(7S)+eta S/K.                                   (24)
```

The second term is polynomial, while `P/S` is a nonzero constant times
`1/(ev^2)`.  Since `Q` is polynomial, `(24)` forces

```text
p=ev^2R                                                (25)
```

for some `R in k[e]`, and then

```text
Q=(5delta/(7alpha))R+(eta alpha/beta)e^2v^3.          (26)
```

The pure `r^3` equation after the differences `(4)` is

```text
0=-4e^2Kf'+4e^2fK'-6eKS'+12eSK'+2KS
  -5ePg'+3egP'-3eQp'+5epQ'+Pg-Qp.                    (27)
```

Let `rho!=0` be a root of `v` of multiplicity `m`, and let
`a=ord_rho(R)`.  The arm law in `(5)` gives

```text
D(rho)f(rho)=-1/12,        hence D(rho)f(rho)!=0.     (28)
```

After `(25)--(26)`, the `f,K` forcing in `(27)` has exact order `4m-1`.
The `R`-square block has order `2a+2m-1`, with leading coefficient a
nonzero scalar times `a-3m`.  All other terms occur later.  If `a<m`, the
`R` block is uniquely earliest and `a-3m!=0`; if `a>m`, the forcing is
uniquely earliest.  Thus

```text
ord_rho(R)=m.                                         (29)
```

At the origin, write `d=ord_0(v)` and `a=ord_0(R)`.  Here `f(0)=1/12`.
When `d>0`, the two corresponding orders are `3+4d` and `1+2d+2a`;
the possible coefficient resonance `a=3d+2` lies strictly above their
unique matching value, so `a=d+1`.  When `d=0`, a later `Pg` term can tie
the forcing at order three, and that order comparison alone is not exact.
The next full canonical bucket instead gives

```text
[r^2z]({A,C}-1)|_(e=0)=-5delta R(0)/(14alpha),

ord_0(R)>=d+1,       with equality when d>0.           (30)
```

Thus `R(0)=0` in the remaining `d=0` case.

Because `k` is algebraically closed, `(29)--(30)` imply

```text
R=evU                                                  (31)
```

for a polynomial `U`.  Once Sections 6--7 exclude constant `v`, either a
nonzero root and `(29)` or a positive origin order and the exact part of
`(30)` forces `R`, hence `U`, to be nonzero.  Equations
`(25)--(26),(31)` are exactly `(6)`.

## 4. The generic terminal square payment

Substitution of `(5)--(6)` in `(27)` gives the following exact law:

```text
0=
 14alpha^2 beta^2 e^3v^7(3ev'+v)
 +21alpha^2 eta e v^2(2eUv'-evU'+Uv)
 +2beta(2ev'+v)(28alpha beta f-5delta U^2)
 +ev(-28alpha beta^2 f'+10beta delta UU')
 +35alpha beta delta v((3ev'+2v)g-evg').              (32)
```

At a nonzero root `rho` of multiplicity `m`, all terms except the third
have order at least `m`.  The coefficient of order `m-1` is

```text
4beta rho m u(rho)
 (28alpha beta f(rho)-5delta U(rho)^2),                (33)
```

where `v=(e-rho)^m u` and `u(rho)!=0`.  This proves the square payment `(7)`.
Notice that `(28)` makes the payment nonvacuous: in particular `U(rho)!=0`.

## 5. The P=0 branch has a 7-3-1 skip ladder

Now suppose `P=0`.  Equation `(22)` reduces to

```text
7eS Q'-3eQ S'-2QS=0.                                  (34)
```

The possibility `Q=0` is empty.  Indeed `(27)` then factors as

```text
2beta e^3v^3
 (3alpha e^4v^7v'+alpha e^3v^8
  +8efv'-2evf'+4fv)=0.                                (35)
```

At a nonzero root of `v`, the `8efv'` term is uniquely earliest and
nonzero by `(28)`.  If `v` has no nonzero root, algebraic closure gives
`v=ce^d`; at the origin the leading coefficient in parentheses is
`(8d+4)f(0)c`, again nonzero.  Thus `Q!=0`.

Equation `(34)` gives `(Q^7/(e^2S^3))'=0`, hence

```text
Q=theta e^2v^3                           (theta in k*). (36)
```

Compare the `Q,p` block in `(27)` with the same `f,K` forcing.  At a
nonzero root of multiplicity `m`, its order is
`ord(p)+3m-1`; equality with `4m-1` forces `ord(p)=m`.  The only coefficient
resonance is `ord(p)=5m`, strictly above the match.  At the origin, equality
of `ord(p)+2+3d` and `3+4d` forces `ord_0(p)=d+1`; the resonance
`3+5d` again lies above it.  Hence

```text
p=evU                                                   (37)
```

with nonzero polynomial `U`, proving `(8)`.

Substitution in `(27)` leaves the exact terminal law

```text
0=
 2beta(3alpha e^4v^7v'+alpha e^3v^8
       +8efv'-2evf'+4fv)
 -3theta(-4eUv'+evU'-2Uv).                            (38)
```

At a nonzero root `rho`, its order-`m-1` coefficient is

```text
4rho m u(rho)(4beta f(rho)+3theta U(rho)),             (39)
```

which proves the linear payment `(9)`.

## 6. The common tower parameter is nonconstant

It remains to exclude a constant `v`.  Absorb its nonzero value into the
three tower constants.  The arm equation in `(5)` is then

```text
(3e^2-2mu e-1)f=2beta e^3-1/12.                       (40)
```

Exact division by the displayed quadratic forces

```text
4mu^2+3=0,                  beta=-mu/4,
f=1/12-mu e/6.                                         (41)
```

In the generic branch, substitute `(41)` and `(6)` into the next untouched
coefficient.  Its origin value is

```text
[e^0]([r^2z^2]({A,C}-1)/e^3)=5delta/2 !=0,            (42)
```

contradicting `(3)`.

In the degenerate branch, put `U_0=U(0)`.  The same coefficient and then
the terminal `r^3` coefficient give

```text
[e^0]([r^2z^2]/e^3)=-6mu U_0,
[e^0]([r^3]/e^3)=18theta U_0-mu/2.                    (43)
```

The first identity forces `U_0=0`, since `(41)` makes `mu!=0`; the second
then leaves `-mu/2!=0`.  Thus constant `v` is impossible in both branches,
as asserted in the theorem statement.

## 7. The linear tower is also empty

Suppose `deg v=1`.  Absorb its leading coefficient into the tower constants
and write

```text
v=e-t.                                                  (44)
```

First let `t!=0`.  In the generic branch, `(28)` and the square payment
`(7)` make both `f(t)` and `U(t)` nonzero.  But direct evaluation of the
next normal coefficient gives

```text
[r^2z]({A,C}-1)|e=t
 =60delta t^2U(t)f(t)/(7alpha),                         (45)
```

which cannot vanish.

In the degenerate branch, the linear payment `(9)` again makes `U(t)!=0`.
The same coefficient is now

```text
[r^2z]|e=t=t^2(-2mu+3t)U(t).                           (46)
```

Thus `t=2mu/3`.  Yet division of the arm numerator
`2beta e^3(e-t)^4-1/12` by `D` has linear remainder coefficient
`2beta/27` at precisely this address, contradicting `beta!=0`.

It remains to recompute the confluent root `t=0`; none of the preceding
steps divides by `t` here.  In the generic branch, the arm linear remainder
requires

```text
64mu^6+240mu^4+216mu^2+27=0.                          (47)
```

The arm constant remainder already excludes `mu=0`.  The origin `r^2`
bucket first gives `mu g(0)=0`, hence `g(0)=0`.  The origin `z^2` bucket
then requires, after cancelling the nonzero `2beta mu`,

```text
16mu^4+44mu^2+21=0.                                   (48)
```

As polynomials in `x=mu^2`, the resultant of the left sides of
`(47),(48)` is

```text
-4534272 !=0.                                          (49)
```

So the generic confluent branch is empty.

Finally take the degenerate branch at `t=0` and put `U_0=U(0)`.  The first
two relevant leading coefficients are

```text
[e^2]([r^2z])=-3mu U_0,

[e^7]([r^3])
 =2beta(256beta mu^5+768beta mu^3+432beta mu)/81
   +54theta U_0.                                      (50)
```

The first forces `U_0=0`.  The arm constant remainder is

```text
256beta mu^5+768beta mu^3+432beta mu-243=0,           (51)
```

so the second coefficient in `(50)` becomes `6beta!=0`.  This closes the
last linear case.  Together with Section 6 it proves `deg v>=2`.

## 8. The next bucket closes every remaining degree

The preceding constant and linear computations are useful independent hostile
controls, but the untouched `r^2z` bucket gives an all-degree closure.

### 8.1 Generic branch

After `(5)--(6)`, direct reduction of that bucket modulo `v` is

```text
[r^2z] = (60delta/(7alpha)) e^2 f U v'        modulo (v). (52)
```

Let `rho!=0` be a root of `v` of multiplicity `m`, and write
`v=(e-rho)^m u` with `u(rho)!=0`.  Its coefficient of order `m-1` is

```text
(60delta/(7alpha)) rho^2 f(rho)U(rho)m u(rho).        (53)
```

The arm identity `(28)` makes `f(rho)` nonzero, and the square payment `(7)`
then makes `U(rho)` nonzero.  Thus `(53)` cannot vanish.  The generic branch
has no nonzero root of `v`.

Since `k` is algebraically closed, `v=c e^d`.  Section 6 excludes `d=0`, so
take `d>=1` and absorb `c` into the nonzero tower constants.  Put

```text
N=3+4d,
D=3e^2-2mu e-1,
D f=2beta e^N-1/12.                                  (54)
```

If `mu=0`, the two roots of `D` are opposite.  Polynomial divisibility in
`(54)` would require their `N`th powers to agree, impossible because `N` is
odd.

Suppose `mu!=0`.  Every tower correction beyond `f,g` has sufficiently high
origin order that the next two canonical buckets give

```text
[r^2]|e=0=mu g(0),
[z^2]|e=0=-3f(0)-mu f'(0)-9g(0).                     (55)
```

The constant and linear jets of `(54)` are

```text
f(0)=1/12,                  f'(0)=-mu/6.              (56)
```

Equations `(55)--(56)` force `g(0)=0` and then

```text
mu^2=3/2.                                             (57)
```

Let `a,b` be the roots of `D`.  They are distinct and nonzero, with
`a+b=2mu/3` and `ab=-1/3`.  Evaluating `(54)` at both roots forces
`a^N=b^N`, so `q=a/b` is a root of unity.  On the other hand,

```text
q+q^(-1)=(a^2+b^2)/(ab)=-(4mu^2/3+2)=-4.             (58)
```

Hence `q^2+4q+1=0`.  Under a complex embedding its two roots are real, while
the only real roots of unity are `+1,-1`; neither satisfies this quadratic.
This contradiction closes the generic branch.

### 8.2 Degenerate branch

In the branch `(8)`, reduction of the same bucket modulo `v` is

```text
[r^2z] = e^2(3e-2mu)U v'                    modulo (v). (59)
```

At a nonzero root `rho` the linear payment `(9)` and `(28)` make `U(rho)`
nonzero.  The order-`m-1` coefficient in `(59)` first forces

```text
3rho-2mu=0.                                           (60)
```

This apparent address is not a survivor.  Expand the full bucket one order
farther at `(60)`.  Its coefficient of `(e-rho)^m` is exactly

```text
(3/2)rho^2 U(rho)u(rho)(2m-5).                        (61)
```

All nonlinear `v` terms begin at order at least `3m-1>m`.  Since `m` is a
positive integer, `2m-5` is nonzero, contradicting `(61)`.  Thus the
degenerate branch also has no nonzero root of `v`.

Again write `v=c e^d`.  Section 6 handles `d=0`; absorb `c` and take `d>=1`.
The coefficient of `e^d` in the terminal law `(38)` is

```text
(2/3)(2d+1)(beta+9theta U(0))=0,                      (62)
```

so `U(0)=-beta/(9theta)!=0`.  If `mu!=0`, the coefficient of `e^(d+1)` in
the `r^2z` bucket is

```text
-mu(2d+1)U(0),                                        (63)
```

which is nonzero.  If `mu=0`, the opposite-root/odd-exponent contradiction
from `(54)` applies unchanged.  This closes the degenerate branch and proves
the asserted nonexistence.

## 9. Scope and exact controls

The theorem identifies the smallest new information that the `rz^2` layer
can carry.  THM-3814 failed because a fourth-power forcing arrived one
valuation step before every fifth-power correction.  Here the successive
top buckets force the debt-paying ladder

```text
v^7 -> v^5 -> v^3 -> v                               (64)
```

in the generic branch, with a precisely typed `e^4,e^3,e^2,e` companion.
The terminal square `(7)` is the new sidecar absent from the `rz` profile,
but `(52)` shows why it cannot pay: the next bucket differentiates the common
root one order too early.  In the degenerate branch, the apparent address
`rho=2mu/3` is destroyed one order later by the half-integral resonance
`m=5/2`.  The final obstruction is therefore not merely a bounded-degree
failure; it is the incompatibility between integral divisor multiplicities
and the odd Kummer ladder.

This theorem closes exactly the eight-profile ansatz `(2)`.  It does not
exclude arbitrary polynomial pairs on `B`, a second `rz^j` layer, a different
arm immersion, or a planar-JC counterexample by another mechanism.

The companion named in the metadata recomputes the full arbitrary-function
bracket, freezes all six source buckets, checks every Kummer valuation
family, verifies the integrated relation `(24)`, reconstructs both terminal
laws from the actual `r^3` coefficient, checks the local payment
coefficients `(33),(39)`, and replays both constant-tower contradictions
`(42)--(43)` and every separated linear-root gate `(45)--(51)`.  Its final
all-degree block independently checks `(52),(59)`, the addressed coefficient
`(61)`, the origin jets `(55)--(57)`, the quadratic root-ratio law `(58)`,
and the monomial payment conflict `(62)--(63)`.  It uses no finite-field or
bounded-degree inference.  **QED.**
