---
id: THM-3821
title: "Cubic pseudo-plane RZ2 profiles enter an odd ladder"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  Any Darboux pair in the first rz^2 extension
  of THM-3814 must enter one of two exact Kummer anatomies.  The generic
  branch has the odd exponent ladder 7,5,3,1 and ends in a Riccati square
  payment; the P=0 branch has the skip ladder 7,3,1 and ends in a linear
  payment.  All one-sided and zero-terminal top branches are empty, and the
  common Kummer parameter must be nonconstant.  This is a necessary
  top-branch classification, not an existence theorem, a complete rz^2
  no-go, or a planar-JC counterexample claim.
source: jc_zero_debt_lift / cubic-pseudoplane rz2 odd-ladder lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The deterministic companion has 46 active
  gates checking the Poisson Casimir and signs, unique monic reduction, six
  descending source buckets, the asymmetric 7/4 contradiction, all 7/4,
  7/5, and 7/3 valuation families, the integrated polynomial relation, the
  generic and degenerate local valuation transfers, both terminal laws, and
  their square/linear nonzero-root payments, and the independent generic
  and degenerate constant-tower contradictions.  Normal/-O/frozen/hash/docs
  replay and independent proof rederivation are required before promotion.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3814-nodal-rz-kummer-profile-degree-gate
related:
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_cubic_pseudoplane_rz2_odd_ladder_thm3821.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_rz2_odd_ladder_thm3821.out
script_sha256: a43c03b0aeacbc9422cbd9f5248f87c3cccadd52ad5aea473ef7234b99373d1e
output_sha256: 7d2d1ce8e91e76df9bc7ed90df3b09e69e5d9054d5e3f7d2848b9f41fd07cf69
semantic_sha256: e32d280719d0be2daf6d0d5315ffb3cc8633894a3a26578bb85d0c9a83fd1783
hash_basis: raw LF bytes
---

# THM-3821 -- the RZ2 layer exposes the 7-5-3-1 ladder

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Let `k` be an algebraically closed field of
characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles `f,g,h,kappa,p,q,S,T in k[e]`, define

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T.              (2)
```

Suppose

```text
{A,C}=1.                                               (3)
```

Then there are `mu in k`, nonzero `alpha,beta in k`, and a nonconstant
`v in k[e]` such that, on writing

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

Moreover exactly one of the following necessary branches occurs.

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

These conclusions are necessary only.  The theorem neither constructs nor
excludes solutions of the remaining terminal laws.

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

## 2. The top profile is necessarily two-sided

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
in `(5)`.  The polynomial `K` cannot vanish: otherwise the nonconstant
quadratic `D` times `f` would be `-1/12`.

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
The two corresponding orders are `3+4d` and `1+2d+2a`; the possible
coefficient resonance `a=3d+2` lies strictly above their unique matching
value.  Therefore

```text
ord_0(R)=d+1.                                         (30)
```

Because `k` is algebraically closed, `(29)--(30)` imply

```text
R=evU                                                  (31)
```

for a polynomial `U`.  The zero polynomial is excluded by the same
earliest-order forcing.  Equations `(25)--(26),(31)` are exactly `(6)`.

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

## 7. Scope and exact controls

The theorem identifies the smallest new information that the `rz^2` layer
can carry.  THM-3814 failed because a fourth-power forcing arrived one
valuation step before every fifth-power correction.  Here the successive
top buckets force the debt-paying ladder

```text
v^7 -> v^5 -> v^3 -> v                               (44)
```

in the generic branch, with a precisely typed `e^4,e^3,e^2,e` companion.
The terminal square `(7)` is the new sidecar absent from the `rz` profile.
The degenerate branch deletes the `v^5` rung but retains a linear payment.

Neither `(32)` nor `(38)` has been solved together with the lower canonical
buckets.  Therefore this result is not an existence theorem and does not
close `(2)`.  It says exactly where a positive construction must live and
which local component spectra it must synchronize.

The companion named in the metadata recomputes the full arbitrary-function
bracket, freezes all six source buckets, checks every Kummer valuation
family, verifies the integrated relation `(24)`, reconstructs both terminal
laws from the actual `r^3` coefficient, checks the local payment
coefficients `(33),(39)`, and replays both constant-tower contradictions
`(42)--(43)`.  It uses no finite-field or bounded-degree inference.  **QED.**
