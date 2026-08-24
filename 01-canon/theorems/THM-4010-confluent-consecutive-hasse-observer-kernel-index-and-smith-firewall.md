---
id: THM-4010
title: "Confluent consecutive Hasse observer kernel, index, and Smith firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every integer
  translate of m+1 consecutive nodes and every jet order k, equality of the
  first k Hasse jets is equivalent to congruence modulo the kth power of the
  monic node polynomial. Consequently an off-node value is determined modulo
  the optimal kth power of the falling-factorial value. On degree below
  N=(m+1)k, the jet-image lattice has index
  product_(j=1)^m (j!)^(k^2), is translation-invariant, and has mod-p rank
  k*min(m+1,p). The tempting Smith list `(j!)^k` repeated k times is false:
  its rank-minimal hostile is (m,k)=(3,2), where the actual final factors are
  12,108 rather than 36,36. The positive higher-p-power Smith exponents remain
  OPEN; no cross-frontier arithmetic or dynamical consequence is asserted.
source: root + independent no-import audit, 2026-08-24
audit: >
  PASS (independent standard-library reconstruction, 2026-08-24). The audit
  rebuilds Hasse matrices, Bareiss determinants, an exact unimodular Smith
  reducer, exhaustive determinantal divisors through rank nine, translation
  and mutation controls, polynomial division, endpoint and values-only
  hostiles, mod-p ranks, and the arbitrary-two-node k=2 formula. It computes
  28 Smith forms through rank 28 and 516,244 exact gates. Primary and audit
  normal/optimized streams match their frozen outputs.
depends_on:
  - THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas
related:
  - THM-4001-cyclotomic-factorial-all-arity-decoder
script: 04-computation/confluent_consecutive_hasse_observer_thm4010.py
output: 05-knowledge/results/confluent_consecutive_hasse_observer_thm4010.out
script_sha256: f28dfba947010cb1b82891b1c1b6981fbc2f0555865840eddb75693d29c10888
output_sha256: 2a97b20bea13267d3a1d324975f921e9399e59861f2206fd9c3a72ce1108cc7e
semantic_sha256: ac4767298585cf7aaed247899cec0e675837b55d57d279d0295681484ee7c0fa
independent_audit_script: 04-computation/confluent_consecutive_hasse_observer_thm4010_independent_audit.py
independent_audit_output: 05-knowledge/results/confluent_consecutive_hasse_observer_thm4010_independent_audit.out
independent_audit_script_sha256: b122225779dc95260c3e6681b932b4ce6a2a48da0304910ad0f23386f8d9bcc1
independent_audit_output_sha256: 42e6daa3eab5a9a6126aa5f7897a311794e2887a0c0f9b062bd5c0bf0a988236
independent_audit_semantic_sha256: 23c5b3fbbaeb75be24ffb5d2221abaff12697bade93a5b0d199c9f11c8d97af4
hash_basis: raw LF bytes
---

# THM-4010 -- the confluent consecutive Hasse observer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Fix

```text
a in Z,                 m>=0,                 k>=1,
x_j=a+j (0<=j<=m),      N=(m+1)k,
F_(a,m)(X)=product_(j=0)^m (X-x_j).                    (1)
```

Write `D^[r]` for the Hasse derivative, normalized by

```text
P(x+T)=sum_(r>=0) D^[r]P(x) T^r,                      (2)
```

and define the integral jet observer

```text
J_(a,m,k)(P)=(D^[r]P(x_j))_(0<=j<=m,0<=r<k).          (3)
```

This is the higher-jet frontier left by [THM-4000](THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas.md).

## 1. Exact kernel, remainder, and optimal new-value modulus

For all `P,Q in Z[X]`,

```text
J_(a,m,k)(P)=J_(a,m,k)(Q)
iff P-Q in F_(a,m)(X)^k Z[X].                          (4)
```

Indeed, Hasse Taylor expansion at `x_j` says that the first `k` jets of
`H=P-Q` vanish exactly when `(X-x_j)^k` divides `H`. The distinct linear
factors are coprime over `Q[X]`, hence their product `F^k` divides `H` there;
monicity puts the quotient back in `Z[X]`. The reverse direction is immediate.

Monic division therefore gives a unique `H_P in Z[X]` of degree `<N` with the
same observed jets as `P`. If `B` is an integer outside the node set, then

```text
P(B)==H_P(B) (mod |F_(a,m)(B)|^k).                    (5)
```

This modulus is optimal, not merely sufficient: the exact ideal of possible
value differences is

```text
{H(B):J_(a,m,k)(H)=0}=F_(a,m)(B)^k Z,                 (6)
```

and `H=F^k` attains its generator. At a sampled node the order-zero jet gives
exact equality; no modulus zero is used. Retaining values alone when `k>1`
is insufficient: `F` vanishes at every node but has a nonzero first jet and
does not belong to `(F^k)`.

## 2. Exact lattice index and translation invariance

Restrict `(3)` to `Z[X]_<N`. In the monomial basis its square matrix is

```text
M_((j,r),q)=binom(q,r)(a+j)^(q-r).                     (7)
```

For `0<=j<=m, 0<=r<k`, use the monic Hermite--Newton basis

```text
G_(j,r)(X)=product_(i<j)(X-x_i)^k (X-x_j)^r.          (8)
```

Its degrees are `jk+r=0,...,N-1`, so the change from monomials is integral
unitriangular. The jet matrix in this basis is triangular, with diagonal

```text
D^[r]G_(j,r)(x_j)=product_(i<j)(x_j-x_i)^k=(j!)^k.    (9)
```

Consequently the jet image is a full-rank sublattice of `Z^N` of index

```text
|det M|=product_(j=1)^m (j!)^(k^2).                   (10)
```

Translation `X -> X+a` is an integral binomial unitriangular coefficient
change, so the complete Smith form, not only the determinant, is independent
of `a`. Hasse normalization is load-bearing: ordinary derivative row `r`
multiplies the determinant by `r!`; the total ratio is
`product_(r=0)^(k-1)(r!)^(m+1)`, first nontrivially `8` at `(m,k)=(2,3)`.

## 3. Uniform Smith information that survives

Let `d_1|...|d_N` be the Smith invariant factors and put
`alpha_(p,i)=v_p(d_i)`. For every prime `p`,

```text
rank_(F_p)(M mod p)=k min(m+1,p),                      (11)
sum_i alpha_(p,i)=k^2 sum_(j=1)^m v_p(j!).             (12)
```

For `(11)`, consecutive nodes collapse modulo `p` to `min(m+1,p)` distinct
classes, giving the upper bound. One representative jet block per class is
the Hasse-confluent CRT map at distinct field points, hence has full rank and
gives equality. Equation `(12)` is the valuation of `(10)`. Thus no prime
`p>m` occurs, and exactly `k min(m+1,p)` invariant factors have zero
`p`-valuation.

It follows that the number of global unit factors is

```text
N      if m<=1,
2k     if m>=2.                                        (13)
```

For the second line, the adjacent-node `2k`-minor is unimodular, while `(11)`
at `p=2` forbids another unit.

There is one further closed boundary case. At arbitrary nodes `0,d`, `d>0`,
with two Hasse jets at each node,

```text
SNF=diag(1,1,d*gcd(d,2),d^3/gcd(d,2)).                (14)
```

Clearing the first node leaves the block

```text
[d^2 d^3; 2d 3d^2],                                   (15)
```

whose entry gcd is `d*gcd(d,2)` and determinant is `d^4`.

## 4. The repeated-factorial Smith guess is false

The triangular diagonal `(9)` determines the determinant, but it is not in
general the Smith diagonal. The rank-minimal failure is `(m,k)=(3,2)`, `N=8`:

```text
naive   (1,1,1,1,4,4,36,36),
actual  (1,1,1,1,4,4,12,108).                         (16)
```

The exact determinantal divisors are

```text
(Delta_0,...,Delta_8)=(1,1,1,1,1,4,16,192,20736).    (17)
```

A direct hostile uses rows `(0,1,2,3,4,5,7)` and columns
`(0,1,2,3,4,5,6)` of `(7)`: its determinant is `2112`. The naive law would
make every rank-seven minor divisible by `576`, but `576` does not divide
`2112`. All smaller ranks are ordinary (`k=1`), unimodular (`m<=1`), or the
sole mixed corner `(m,k)=(2,2)`, whose actual Smith form is
`(1,1,1,1,4,4)`; hence the hostile is rank-minimal. Another small failure is

```text
(m,k)=(2,3): actual (1,1,1,1,1,1,2,16,16),
             naive  (1,1,1,1,1,1,8,8,8).             (18)
```

## 5. Exact finite audit and open boundary

The primary companion computes the 24 pairs `0<=m<=5,1<=k<=4`. The no-import
standard-library audit expands this to all 28 pairs

```text
0<=m<=6,                 1<=k<=4,                     (19)
```

through rank 28. It independently checks Bareiss determinants, exact
unimodular Smith reduction, every determinantal divisor in all 16 cases of
rank at most nine, three translations, three unimodular mutations, mod-`p`
ranks for `p=2,3,5,7`, twelve instances of `(14)`, polynomial division,
sampled and off-node endpoints, and the values-only hostile. Normal and
optimized streams match the frozen outputs across `516,244` exact gates.

What remains **OPEN** is a closed formula for the ordered positive exponents

```text
alpha_(p,kp+1),...,alpha_(p,N)                         (20)
```

when `p<=m`, equivalently the higher-`p^e` filtration of the consecutive
Hasse-jet image. Equations `(10)`--`(12)` determine their number and total
sum, not their distribution; `(16)` is the first obstruction. This theorem
does not transfer to Rule 30 phase dynamics, LRC(14), the planar Jacobian
problem, factorial-coefficient conjectures, class ranks, or Hopf problems.
