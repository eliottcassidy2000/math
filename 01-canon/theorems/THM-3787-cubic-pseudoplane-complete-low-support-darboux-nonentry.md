---
id: THM-3787
title: "Cubic pseudo-plane complete low-support Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  THM-3785 cubic Russell pseudo-plane, every Darboux pair with two Euler
  weights in one output and three in the other is impossible, including the
  common-step arithmetic-progression cells.  Every common-step three-by-three
  AP cell whose scalar lies in the central three-term bracket bucket is also
  impossible.  The proof is uniform in all weights and closes the negative
  divisor, weight-zero, repeated-degree, and A=1 radical seams.  General
  two-by-four and noncentral or non-AP three-by-three supports remain open;
  no planar Jacobian counterexample is claimed.
source: jc_zero_debt_lift / pseudo-plane Euler-support census, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The sign/zero boundary,
  every nominal-degree cancellation direction, the finite-radical extension
  and every place above infinity, and both a=1 seams were rederived from the
  displayed equations.  One overstrong sentence about constant profiles at
  positive weights divisible by three was repaired; the scalar-bucket leading
  coefficient is nonzero there as well, so the proof and theorem are
  unchanged.  The deterministic companion verifies the
  complete two-by-three collision grammar, every bounded sign/zero seam as a
  hostile control, both common-step resonance identities and their exceptional
  leading-cancellation alternatives, the A=1 first-bucket ODE, the two
  resonant scalar factorizations, the three-by-three adjacent integrations,
  radical elimination, and the a=1 Wronskian formula.  Normal and optimized
  runs byte-match the frozen transcript.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3781-common-step-three-by-three-danielewski-darboux-nonentry
  - THM-3786-quadratic-etale-surface-irregular-two-by-three-support-no-go
script: 04-computation/jc2_cubic_pseudoplane_low_support_thm3787.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_low_support_thm3787.out
script_sha256: d42f6593cbb9b30ca2708716f2f4de6c2e0e57f8ec6e5426a7e506b884a71f16
output_sha256: 73c929260d4c5fa6cdcc9b5aad67124727c3ebb99c993fdaa89454abd40291b5
semantic_sha256: 3ad9952bfea2c2d9a51723cbcff23ecbf7d1d166d7a09c568f4cc675344f6fe7
hash_basis: raw LF bytes
---

# THM-3787 -- the cubic pseudo-plane has no low-support Darboux entrance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
an algebraically closed field of characteristic zero, let `c in k*`, and let

```text
B=k[R,Z,E]/(R^2E-Z^3+c^3R)                           (1)
```

carry the Poisson bracket from THM-3785:

```text
{R,Z}=3R^2,        {R,E}=9Z^2,        {Z,E}=3c^3+6RE. (2)
```

Give `(R,Z,E)` weights `(3,1,-3)`.  The bracket has weight shift `+2`.
This theorem proves the following two uniform statements.

1. If `F,G in B` have exactly two and three active Euler weights,
   respectively, then `{F,G}!=1`.  The same holds after swapping the outputs.
2. If both supports are three-term arithmetic progressions with common step
   and the scalar bracket is carried by their central three-contribution
   bucket, then `{F,G}!=1`.

Here an active weight-zero scalar may always be subtracted and is not counted
as support.  A nonconstant weight-zero profile is retained and treated below.
The theorem does **not** close `2 by 4`, non-AP `3 by 3`, or AP `3 by 3`
cells with a noncentral scalar bucket.

## 1. Inherited homogeneous calculus

Use THM-3785's source chart

```text
w=c+x^3y,              Delta=w^3-c^3,
R=x^3,                 Z=xw,                 E=x^(-3)Delta.          (3)
```

For `u in Z`, write `rho(u) in {0,1,2}` for its residue modulo three and

```text
m(u)=max(0,ceil(-u/3)).                                (4)
```

The exact weight module is

```text
B_u=x^u w^rho(u) Delta^m(u) k[w^3].                   (5)
```

In particular every negative-weight profile is nonconstant and vanishes at
all three roots of `Delta`.  For homogeneous elements

```text
A=x^a f(w),                  C=x^b g(w),               (6)
```

put

```text
[a,f;b,g]=a f g'-b f'g.                               (7)
```

Then

```text
{A,C}=x^(a+b+2)[a,f;b,g].                             (8)
```

We use three immediate consequences repeatedly.

- A singleton contribution bucket cannot equal a nonzero scalar: that would
  be the homogeneous Darboux pair excluded by THM-3785.
- If two nonzero-weight homogeneous terms commute, their weights have the
  same sign.  If one weight is zero and the other is nonzero, the zero-weight
  profile is scalar and can be removed.  Indeed `(7)=0` is a logarithmic
  derivative identity; opposite signs would make a product of nonconstant
  polynomials constant.
- If nonzero same-sign weights `u,v` commute, unique factorization gives

```text
f=lambda p^(|u|/g),       g=mu p^(|v|/g),
g=gcd(|u|,|v|).                                           (9)
```

The letters `lambda,mu` in `(9)` are nonzero field constants.  These facts
also type every weight-zero boundary below; no use is made of the false claim
that all weight-zero elements are scalar.

## 2. The complete irregular two-by-three collision grammar

Write the exact supports as

```text
supp(F)={a,a+r},             supp(G)={b,b+s,b+t},
r>0,                         0<s<t.                   (10)
```

After subtracting the bottom contribution weight, the six bracket addresses
are

```text
0, s, t, r, r+s, r+t.                                 (11)
```

Thus a collision occurs only when

```text
r=s,                         r=t,              or     t=r+s.          (12)
```

There are two collisions exactly when `r=s` and `t=2r`, the common-step AP.
Without a collision the scalar bucket is a forbidden singleton.  Suppose now
that exactly one relation in `(12)` holds.  The scalar must lie there.

### 2.1 The lower collision `r=s`, `t!=2r`

The scalar condition is `a+b=-r-2`.  The singleton bottom endpoint makes
`a,b<0`.  Since `t>r` and the non-AP condition excludes the only case `t=2`
with `r=1`, one has `t>=3`.  The singleton top endpoint has total weight
`t-2>0`, so both `a+r` and `b+t` are positive.  But the singleton at address
`t` asks the negative weight `a` to commute with the positive weight `b+t`,
contradicting Section 1.

### 2.2 The upper collision `r=t`

Again `a+b=-r-2` and `a,b<0`.  The singleton at `s` gives `b+s<0`.  If
`r>=3`, the top endpoint has total weight `r-2>0`, so `a+r>0`; the singleton
at `r+s` then has opposite signs.  If `r=2`, necessarily `s=1`.  The top
weights sum to zero, hence both are zero; the singleton at `r+s` forces the
weight-zero `F` profile commuting with weight `-1` to be scalar and removable.
The case `r=1` has no positive `s<r`.

### 2.3 The diagonal collision `t=r+s`, `r!=s`

The scalar condition is `a+b=-r-s-2`, so `a,b<0`.  The singleton addresses
`s` and `r` force `b+s<0` and `a+r<0`.  If `r>=3`, the top endpoint has
positive total weight `r-2`, contradicting `a+r<0`.  If `r=2`, its total is
zero and the singleton at `r` forces the resulting weight-zero profile to be
scalar and removable.  If `r=1`, the top endpoint weights are commuting
integers with sum `-1`; they cannot have the same strict sign, and a zero
weight would again be scalar and removable.  This is also impossible.

This proves that every irregular `2 by 3` and `3 by 2` cell is empty.

## 3. The lower scalar bucket in a common-step two-by-three AP

It remains to take

```text
supp(F)={a,a+d},             supp(G)={b,b+d,b+2d}.    (13)
```

First put the scalar in the lower double bucket, so

```text
a+b+d=-2.                                             (14)
```

For `d=1`, the top weights are both zero and each of the two scalar summands
contains a negative profile as an undifferentiated factor.  The scalar bucket
therefore vanishes at `Delta=0`.  Hence assume `d>=2`.  Endpoint sign
synchronization lets us write

```text
A=-a,        B=d+2-A,        C=d-A,
D=d+A-2,    T=A-2,          1<=A<=d-1.               (15)
```

The low endpoint weights are `(-A,-B)`, the high endpoint weights are
`(C,D)`, and the middle `G` weight is `T`.

The case `A=1` is the seam at which the ordinary leading-degree argument
loses a coefficient.  Endpoint commutation gives profiles proportional to

```text
p, p^(d+1)                         at the low endpoint,
q, q                               at the high endpoint.              (16)
```

If `M` is the middle profile, the nonscalar double bucket is exactly

```text
(d-1)q M'+q'M=q_0[(d-1)q p'+q'p].                    (17)
```

Thus `H=M-q_0p` satisfies

```text
(d-1)qH'+q'H=0.                                      (18)
```

If `H!=0`, then `H^(d-1)q` is constant.  Polynomiality makes `q,H`
constant, while `H` is a difference of two weight `-1` profiles and is
divisible by `Delta`; hence `H=0`.  The scalar bucket is now either zero or
a nonzero multiple of

```text
p^d[(d-1)q p'+q'p],                                  (19)
```

whose leading coefficient has positive degree because `Delta|p`.

Now let `A>=2`.  Put

```text
delta=gcd(A,B),                  epsilon=gcd(C,D),
f=p^(A/delta),   g=lambda p^(B/delta),
P=q^(C/epsilon), Q=mu q^(D/epsilon).                 (20)
```

Here `f,g` are the low endpoint profiles and `P,Q` the high endpoint
profiles.  Write `kappa=deg p>0`, `eta=deg q>=0`, and `n=deg M`.
When `T=0`, a scalar `M` is removable, so then an active middle component has
`n>0`.  If `T>0`, the profile may have `n=0` when the positive weight is
divisible by three, but the scalar bracket still keeps its leading term: its
leading coefficient contains the nonzero term `-T A kappa/delta`.
The scalar bucket and the nonscalar double bucket force the two nominal
degree equalities

```text
A kappa/delta+n=C eta/epsilon+B kappa/delta,
A kappa/delta+D eta/epsilon=C eta/epsilon+n.          (21)
```

There is no hidden leading-cancellation branch.  If the second equality in
`(21)` failed because the `P,M` bracket lost its leading term, then

```text
C n=T C eta/epsilon.                                 (22)
```

Combining `(22)` with the first equality gives `2A=d+2`; but the required
strict nominal-degree inequality has the opposite sign, since then
`3A-d-2=A>0`.  Thus `(21)` is exhaustive.  Subtracting its two equations
gives

```text
(3A-d-2)(kappa/delta+eta/epsilon)=0,                  (23)
```

so necessarily

```text
d=3A-2.                                               (24)
```

At this resonance `(A,B)=(A,2A)` and `(C,D)=(C,2C)` with `C=2A-2`.
After rescaling endpoint terms, the nonscalar bucket integrates to

```text
M=q_0pq+H,                  CqH'=Tq'H.                (25)
```

The scalar bucket factors exactly as

```text
(Cq p'+A p q') (c_0p-c_1H/q).                        (26)
```

At infinity the first factor has degree `kappa+eta-1` and nonzero leading
coefficient `C kappa+A eta`.  If `c_0!=0`, the `p` term dominates because
`deg H=(T/C)eta<eta`; if `c_0=0`, the total degree is

```text
kappa-1+(T/C)eta>0                                   (27)
```

because the negative endpoint gives `Delta|p`.  Thus `(26)` is never one.

## 4. The upper scalar bucket in a common-step two-by-three AP

Put the scalar in the upper double bucket instead:

```text
a+b+2d=-2.                                            (28)
```

For `d=1`, all profiles in both scalar summands have negative weight, so the
scalar vanishes at `Delta=0`.  For `d=2`, each scalar summand still contains
a negative profile as an undifferentiated factor and vanishes there, including
the two nonconstant weight-zero endpoint profiles.  For `d=3`, endpoint signs
would require two positive integers summing to one.  Hence `d>=4`, and

```text
A=-a,      B=2d+2-A,      C=d+2-A,
D=d-A,    T=A-2,          3<=A<=d-1.                 (29)
```

The endpoint profiles have the same form `(20)`, now for pairs `(-A,-B)`
and `(D,T)`, while the middle `G` profile `M` has weight `-C`.
The scalar degree balance is exact.  The other bucket can lose its nominal
leading term only if

```text
A n=C A kappa/delta.                                 (30)
```

If its nominal degree were strictly larger than the fixed term, `(30)` and
the scalar balance would give `2A=d+2`, which makes the supposedly positive
difference `3A-2d-2` equal `2-A<0`.  Thus no exceptional branch exists, and

```text
(3A-2d-2)(kappa/delta+eta/epsilon)=0.                 (31)
```

Consequently `3A=2d+2`.  Write

```text
A=2C_0+2,                    d=3C_0+2,       C_0>=1. (32)
```

The low endpoint profiles are `p,p^2`, the high endpoint profiles are
`q,q^2`, and the other bucket integrates to

```text
M=q_0pq+H,                  A pH'=(C_0+2)p'H.         (33)
```

The scalar bucket factors as

```text
(C_0q p'+A p q') (c_0q+c_1H/p).                      (34)
```

If `c_0!=0`, its infinity degree is positive.  If `c_0=0`, equation `(33)`
gives `H^A=lambda p^(C_0+2)`.  Every root of `p` is therefore a root of
`H`; since `Delta|p`, one has `deg H>=3`, and the degree of `(34)` is
`eta+deg H-1>0`.  This closes the second scalar placement and completes all
`2 by 3` and `3 by 2` supports.

## 5. Scalar-centred common-step three-by-three AP supports

Finally suppose both supports are common-step APs and the scalar is carried
by the central three-contribution bucket.  After swapping the outputs if
needed, their weights are

```text
F:  -(d+a),  -a,  d-a,
G:  -(d-a+2), a-2, d+a-2,            1<=a<=d-1.      (35)
```

The case `d=1` has two weight-zero top endpoints, and every central scalar
summand vanishes at `Delta=0`.  Hence `d>=2`.  The range in `(35)` follows
from endpoint sign synchronization; output swap sends `a` to `2-a` and
selects the displayed representative.

Let `K,N` be endpoint radicals in a finite algebraic extension of `k(w)`:

```text
f_-=A_0 K^(d+a),       g_-=B_0 K^(d-a+2),
f_+=L_0 N^(d-a),       g_+=M_0 N^(d+a-2).             (36)
```

Concretely, if the low and high common UFD bases are `p,q`, take
`K^u=p`, `N^v=q` for the corresponding endpoint gcds.  Because `k` is
algebraically closed, the constant field of this finite extension is still
`k`.  Let `alpha,beta` be the two middle profiles.  The lower and upper
two-contribution buckets integrate exactly to

```text
alpha/K^a=lambda K^(a-2) beta+C_0,
beta/N^(a-2)=mu alpha N^a+D_0.                        (37)
```

Eliminating `beta` gives, with `nu=lambda mu!=0` and `T=KN`,

```text
alpha(1-nu T^(2a-2))
   =K^a(lambda D_0 T^(a-2)+C_0).                     (38)
```

Suppose `a>=2`, and choose a normalized place above infinity.  Write
`delta_K>0` and `delta_N>=0` for the pole orders of `K,N`; the strict first
inequality follows because the negative endpoint profile is divisible by
`Delta`.  If the `D_0` term dominates in `(38)`, the pole order of `alpha`
is

```text
-a delta_N<=0.                                       (39)
```

If that term is absent or cancels, its pole order is at most

```text
(2-a)delta_K-(2a-2)delta_N<=0.                       (40)
```

But `alpha` is a nonzero weight `-a` polynomial profile and therefore has
strictly positive pole order.  This contradiction closes every `a>=2`.

It remains to treat `a=1`, where `(38)` deliberately loses its large
denominator.  Put `n=d+1`, `m=d-1`.  The endpoint profiles in `(36)` are
`K^n,N^m`, and `(37)` becomes

```text
alpha=lambda beta+C_0K,
beta=mu alpha+D_0/N.                                  (41)
```

If `nu=1`, consistency and the positive pole order of `KN` force
`C_0=D_0=0`; the endpoint scalar coefficient and the middle Wronskian then
both vanish.  If `nu!=1`, direct substitution factors the complete central
scalar bucket as

```text
(KN)'[
  L_0B_0(1-nu) n m K^(n-1)N^(m-1)
  + C_0D_0/((1-nu)N^2)
].                                                    (42)
```

The first term in the brackets has strictly larger pole order than the
second by

```text
(n-1)delta_K+(m+1)delta_N>0.                          (43)
```

Moreover, after multiplication by `(KN)'`, that dominant term is the
endpoint polynomial bracket and has degree
`deg(K^n)+deg(N^m)-1>=2`, since the negative profile `K^n` is divisible by
`Delta`.  It cannot cancel to a scalar.  This closes the exceptional radical
seam and proves the scalar-centred `3 by 3` statement.

## 6. Exact remaining support frontier

THM-3785 closed homogeneous and all `2 by 2` supports.  Sections 2--4 now
close every `2 by 3` and `3 by 2` support, and Section 5 closes the most
collision-rich central common-AP `3 by 3` cell.  The first unexcluded sparse
shapes are therefore

```text
2 by 4,
3 by 3 with a noncentral scalar bucket,
3 by 3 outside the common-step AP grammar.             (44)
```

These are support-shape statements only.  They neither construct a Darboux
pair on `(1)` nor prove that none exists.  By THM-3785, any actual survivor
would still pull back through the cubic atlas to a planar Keller map of field
degree divisible by three and at least nine, by THM-3794.

The deterministic companion named in the metadata checks the bucket census,
all small zero-weight controls, the resonance and exceptional-cancellation
identities, both exact scalar factorizations, and the finite-radical
integration/Wronskian packet.  **QED.**
