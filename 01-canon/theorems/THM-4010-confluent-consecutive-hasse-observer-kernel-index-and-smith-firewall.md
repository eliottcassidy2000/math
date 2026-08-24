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
  k*min(m+1,p). Writing n=m+1, when n>p the first positive p-Smith
  exponent is exactly 1+v_p(k); the exponent-one multiplicity is zero if
  p divides k and otherwise min(p,n-p). The complete p-primary partition is
  the residue-cluster multiset union, and for k=2, p<n<=2p it is explicitly
  `0^(2p),1^(n-p),3^(n-p)` for odd p and `0^4,2^(2n-4)` for p=2.
  The tempting Smith list `(j!)^k` repeated k times is false:
  its rank-minimal hostile is (m,k)=(3,2), where the actual final factors are
  12,108 rather than 36,36. Later positive p-power layers remain OPEN; no
  cross-frontier arithmetic or dynamical consequence is asserted.
source: root + independent no-import audit, 2026-08-24
audit: >
  PASS (independent standard-library reconstruction, 2026-08-24). The audit
  rebuilds Hasse matrices, Bareiss determinants, an exact unimodular Smith
  reducer, exhaustive determinantal divisors through rank nine, translation
  and mutation controls, polynomial division, endpoint and values-only
  hostiles, mod-p ranks, and the arbitrary-two-node k=2 formula. It computes
  28 Smith forms through rank 28 and 516,244 exact gates. Primary and audit
  normal/optimized streams match their frozen outputs. A second independent
  modular-DVR referee verifies the first positive layer, CRT cluster union,
  complete k=2 pair band, and both minimal hostiles on 180 atlas cases through
  rank 80 and 9,501,593 gates. It caught and repaired a cancellation-unsafe
  attainment sentence; the canonical proof uses the literal residual entry
  `pk` after clearing the representative-node identity block.
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
primewise_script: 04-computation/confluent_sampler_primewise_extension_thm4010.py
primewise_output: 05-knowledge/results/confluent_sampler_primewise_extension_thm4010.out
primewise_script_sha256: 79f2b53a945815e6729202661c2ae5b045363d99175bef380f0652d52a4c6f49
primewise_output_sha256: 10cb170be545de5b5a946ee24786b9ca7982aac0c76c9cd2869f07c03148db5e
primewise_semantic_sha256: a9e30379316861df868043a67f3306f8a343e474aec052d22889be3eb122a989
primewise_independent_script: 04-computation/confluent_sampler_primewise_extension_thm4010_independent_audit.py
primewise_independent_output: 05-knowledge/results/confluent_sampler_primewise_extension_thm4010_independent_audit.out
primewise_independent_script_sha256: c7c2d6a6c865958d83d1bf49f9f1dd61050496d946d5017e6d2361d06bce996a
primewise_independent_output_sha256: 0ee98c1861df3c9e0af95c5d2249e42f57d88f65be714f7054c09bb22fbb98c7
primewise_independent_semantic_sha256: 0f8dc0d0d56792c2526b4395230393963e77df72c9b102d2b6086fae3ffd0f26
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

## 4. First positive `p`-adic layer and residue clusters

Put `n=m+1`.  If `n<=p`, equation `(11)` gives full rank modulo `p`, so
every `p`-Smith exponent is zero.  Suppose `n>p` and order the exponents as

```text
alpha_(p,1)<=...<=alpha_(p,nk).                        (14a)
```

Then

```text
alpha_(p,1)=...=alpha_(p,kp)=0,
alpha_(p,kp+1)=1+v_p(k).                               (14b)
```

The first line is `(11)`.  For the genuinely `p`-adic step, first group the
nodes by residue modulo `p`.  The monic factors

```text
product_(j congruent c mod p)(X-j)^k                  (14c)
```

have pairwise unit resultants over `Z_(p)`.  Chinese remainder on both the
polynomial quotient and its jet target therefore makes the global
`p`-exponent partition the sorted multiset union of the residue-cluster
partitions.

In one cluster, translate the representative node to zero and clear its
`k`-jet identity block.  Write the remaining nodes as `ph`, `h>=1`.  On the
mod-`p` kernel, the residual entry from coefficient degree `k+s` to jet
`r=k-t`, where `s>=0` and `1<=t<=k`, is

```text
binom(k+s,k-t)(ph)^(s+t).                              (14d)
```

With `L=s+t`, the `L` consecutive numerator factors in
`binom(k+s,L)` contain `k`.  Hence

```text
v_p(p^L binom(k+s,L))
 >=L+v_p(k)-v_p(L!)
 >=1+v_p(k).                                          (14e)
```

This lower bound is attained without any summand-cancellation argument: the
literal cleared-matrix entry `(h,r,q)=(1,k-1,k)` is exactly `pk`.  Since
`n>p` supplies the duplicate node `p`, `(14b)` follows.

The exponent-one multiplicity is also exact:

```text
#{i:alpha_(p,i)=1}=0                         if p divides k,
                     min(p,n-p)              if p does not divide k. (14f)
```

When `p` is prime to `k`, divide the duplicate-node block by `p` and reduce
modulo `p`.  Terms with `L>=2` vanish; at an active residue `c` the remaining
map is, in jet `k-1`,

```text
Q -> k h G'(c)^k Q(c),
G(X)=product_(d=0)^(p-1)(X-d).                         (14g)
```

Exactly `min(p,n-p)` residue clusters contain a duplicate.  Evaluation of
`Q` at those distinct residues has full row rank, since its available
dimension is `k(n-p)`.  This proves `(14f)`; if `p|k`, `(14b)` already
places the first positive exponent above one.

The cluster union closes one complete band.  For `k=2` and `p<n<=2p`, every
cluster has size one or two, so `(14)` at gap `d=p` gives

```text
p odd:  0^(2p), 1^(n-p), 3^(n-p),
p=2:    0^4,    2^(2n-4).                              (14h)
```

The factor `k` is load-bearing: `(p,m,k)=(2,2,2)` has partition
`0^4,2^2`, not an exponent-one factor.  Nor is there one layer-one factor per
extra node once a residue cluster has size three.  The first genuinely
confluent hostile with `k>=2` and `p` prime to `k` is

```text
(p,m,k)=(3,6,2): 0^6,1^3,3^3,4^2,                    (14i)
```

where the four extra nodes produce only three exponent-one factors.

## 5. The repeated-factorial Smith guess is false

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

## 6. Exact finite audit and open boundary

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

The primewise extension separately checks `226` complete prime partitions for
`m<=30`, `k<=8`, and rank through `66`, totaling `16,614,629` gates.  Its
independent `Z/p^E` reducer and subset-DP minor engine checks `180` atlas cases
through rank `80`, the complete pair band through `p=17`, `164` boundaries
with `n<=p`, and `9,501,593` gates.  Both normal/optimized pairs match their
frozen outputs.

What remains **OPEN** is a closed formula for the later positive exponents

```text
alpha_(p,kp+2),...,alpha_(p,N)                         (20)
```

when `p<=m`, equivalently the higher Bockstein filtration internal to clusters
of size at least three.  A finite atlas suggests, but does not prove, a stable
`k=2` local pattern for cluster size below `p`.  Equations `(10)`--`(14f)`
fix the zero layer, first positive layer, exponent-one count, and total mass,
not the later distribution; `(16)` is the first global obstruction. This theorem
does not transfer to Rule 30 phase dynamics, LRC(14), the planar Jacobian
problem, factorial-coefficient conjectures, class ranks, or Hopf problems.
