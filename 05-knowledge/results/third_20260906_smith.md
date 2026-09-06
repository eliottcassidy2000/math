# Every mixed (m,2,1) Smith form: three residual factors and a capped projective residue

**Status: PROVED in the stated scope; FINITE-EXACT controls passed;
[independent full-matrix and proof audit accepted](third_20260906_smith_root_audit.md).** This classifies the entire infinite multiplicity family
`(m,2,1)` for every integer `m>=2`, at every prime, on the declared integral
homogeneous coefficient module. It includes a complete global integer formula
in affine coordinates and an intrinsic projective classification. No theorem
ID, external novelty, or priority claim is made.

The [uniform-unit consumer](third_20260906_smith_density.md) gives exact
mean recovery gains and full-kernel sizes from the complete residue ladder.

## Inheritance, first hostile, and live board

The immediate mechanism is the
[complete `(2,2,1)` theorem](second_20260906_smith.md): after clearing one
complete bank, the residual rank is three, so three determinantal invariants
suffice. The new
[arbitrary projective-jet transport](continuing1_20260906_jets_projective.md)
supplies lawful complete-bank covariance, the weighted bracket determinant,
and residue splitting, with exact denominator attainment. Its generic inverse
formula is inherited rather than re-claimed.

The canonical hostile is the old dyadic equal-metric pair with exponents
`(0,0,4,7,9)` versus `(0,0,4,8,8)`. The corrected near miss is treating a
raw chart-dependent numerator as an intrinsic invariant: the previous theorem
had to mask its inactive digits. The least-used sidecar here is the **ideal
relation between different residual minors**, before taking any valuations.

The six-concept board is: arbitrary complete Hasse multiplicity; fixed residual
rank three; the full second determinantal ideal; the coefficient integers
`m,m+1,m+2`; reference changes in a bracket ratio; and capped finite-precision
kernels. Increasing `m` preserves the first three objects, so this operation
can be solved uniformly instead of adding isolated small Smith banks.

The first `m=3` hostile has four different spectra at one ternary metric.
Swapping the unequal multiplicity labels changes the spectrum, even though
its labelled distance pattern survives. A second hostile changes an infinite
linear cancellation to valuation one under a lawful projective chart. These
probes led to the precise cap below. No incorrect exploratory statement was
promoted. Targeted searches recovered the old `(2,2,1)` law and the incoming
general precision formula, but not this all-`m` full-partition formula.

## 1. All integer invariant factors in one line of gcds

Let `m>=2` and let `a,b` be nonzero distinct integers. On integer polynomials
of degree at most `m+2`, observe Hasse orders `0,...,m-1` at zero, orders
zero and one at `a`, and order zero at `b`. Put

```text
G_1 = gcd(a^m, m*a^(m-1), b^m),
G_2 = gcd(a^(2m), a^m*b^m*(b-a),
                    a^(m-1)*b^m*(m*b-(m+1)*a)),
D   = |a|^(2m)*|b|^m*|b-a|^2.
```

All gcds are positive. The complete integer Smith list is

```text
1 repeated m times,       G_1, G_2/G_1, D/G_2.       (1)
```

The displayed factors satisfy divisibility in this order. Formula `(1)` is
valid for every `m`; no search over higher reciprocal coefficients or over
more than three residual rows is required.

### The mechanism: the apparent quadratic packets are redundant

The first `m` Hasse rows at zero form an identity matrix in degrees
`0,...,m-1`. Clearing those columns by integral row operations leaves

```text
R_m = [ a^m          a^(m+1)        a^(m+2)       ]
      [ m*a^(m-1)    (m+1)*a^m      (m+2)*a^(m+1)]
      [ b^m          b^(m+1)        b^(m+2)       ].          (2)
```

Its entry ideal has gcd `G_1`, and its determinant has absolute value `D`.
For the second ideal, set

```text
L=m*b-(m+1)*a,
U=a^(2m),        V=a^m*b^m*(b-a),
W=a^(m-1)*b^m*L.
```

The three same-node value/derivative minors are `U,2aU,a^2U`. The three
value/value minors are `V,(a+b)V,abV`. The three derivative/value minors are

```text
W,             (a+b)W+V,             abW+bV.        (3)
```

For the latter identities use

```text
m*b^2-(m+2)*a^2 = (a+b)L+a(b-a),
(m+1)b-(m+2)a   = L+(b-a).
```

Thus **all nine** second minors generate exactly `(U,V,W)`, including both
ideal containments because `U,V,W` themselves are minors. This proves `G_2`
and hence `(1)` by determinantal divisors. In particular the apparent
`m+2` term supplies no additional prime wall or residue coordinate.

As an exact precision corollary, for a prime `p` put
`A=v_p(a),B=v_p(b),C=v_p(b-a)`. Then the largest exponent simplifies to

```text
L_p=max(mB+2C, mA+C,
              (m+1)A+2C-v_p(mb-(m+1)a)),           (4)
```

where a zero numerator contributes minus infinity. This follows by subtracting
the minimum of the three second-minor valuations from `2mA+mB+2C`.
The inherited arbitrary-jet inverse formula has many reciprocal coefficients;
here their maximum equals the three elementary candidates in `(4)`. This
collapse is a consequence of the complete residual ideal, not a truncation
of the inverse theorem.

A further consequence is that the worst precision is always attained at a
**value observation in one of the two small banks**. The simple point has
inverse denominator exponent `mB+2C`. The double point's order-zero and
order-one reciprocal coefficients have denominator exponents `mA+C` and
`(m+1)A+2C-v_p(mb-(m+1)a)`. Their maximum is the exact denominator of its
value column, by the primitive-content attainment theorem in
[complete projective higher jets](continuing1_20260906_jets_projective.md).
Their combined maximum is exactly `(4)`. Complete local chart changes preserve
the denominator of a whole bank, and the value column attains that bank's
maximum in every unimodular frame. Thus this observation-owner conclusion
also holds in the projective observer below.

If `L_p>0`, for every `N>=1` there is an integral source perturbation invisible
modulo `p^(N+L_p-1)` in the data but visible modulo `p^N` in its coefficients,
whose only changed datum is the value at the double or simple point. All
jets at the growing `m`-bank remain unchanged. Multiply that cardinal value
column by its exact integer denominator and then by `p^(N-1)`; denominator
attainment proves the required coefficient and data valuations. This is a
uniform-precision ownership theorem, not a claim that the heavy bank can be
deleted or that arbitrary unequal data-precision budgets are unchanged.

## 2. Full projective classification at every prime

Let `v_0,v_1,v_2` be primitive integer vectors in pairwise distinct rational
directions. Give them multiplicities `m,2,1` in that order. On homogeneous
binary forms of degree `m+2`, observe all coefficients

```text
[T^r] F(v_i+T w_i),             0<=r<m_i,
```

where `det(v_i,w_i)=1`. These are complete Hasse banks; ordinary higher
derivatives cannot be substituted without factorial factors.

Fix a prime `p`, set

```text
A=v_p(det(v_0,v_1)),  B=v_p(det(v_0,v_2)),
C=v_p(det(v_1,v_2)),  mu=v_p(m).
```

There are `m` initial zero exponents, followed by precisely three exponents
classified below. The first bank remains tied to its multiplicity `m`.
Permuting distinct multiplicities while retaining the same geometric points
is not an allowed symmetry of the observer.

**Case I: `A<=B`.** Set `h=min(A,mu)`. The three remaining exponents are

```text
((m-1)A+h,      (m+1)A-h,      mB+2C).              (5)
```

This includes a closest pair consisting of the double and simple nodes,
a closest pair consisting of the `m`-fold and simple nodes, and every
all-equal depth pattern.

**Case II: `A>B` and `d=A-B` is different from `mu`.** Put `e=B=C`, and set

```text
rho=(m-1)e+min(e,(m-1)d+mu),
theta=min(d,mu),
sigma=2me+(m-1)d+theta.
```

The three exponents are

```text
(rho,     sigma-rho,     (m+2)e+(m+1)d-theta).       (6)
```

In particular every prime not dividing `m` is metric-only throughout the
family. Divisibility of `m+1` or `m+2` does not create a separate exception.

**Case III: `A>B` and `d=A-B=mu>0`.** Put `e=B=C` and

```text
K=min(e,m*d).
```

If `K=0`, put `kappa=0`. Otherwise choose any primitive local reference
vector `w` with all `det(v_i,w)` units at `p`, and define

```text
tau_w = det(v_0,v_1)*det(v_2,w)
         / [p^d det(v_0,v_2)*det(v_1,w)] in Z_p^*,
kappa = min(K, v_p(m/p^d-(m+1)*tau_w)).              (7)
```

At a zero argument, take the capped valuation to be `K`. A suitable reference
exists because `e>=1` puts all three directions in the same residue class;
a unimodular tangent to `v_0` suffices. The capped valuation is independent
of that reference. The three exponents are

```text
((m-1)e+K,
 (m+1)e+m*d-K+kappa,
 (m+2)e+m*d-kappa).                                 (8)
```

These are the only potentially nonmetric factors. The first residual factor
is metric; the last two shift oppositely as `kappa` varies. All displayed
lists are sorted. For example in `(8)`, their final difference is
`e+K-2kappa>=e-K>=0`; the other difference is nonnegative as well.

This is the complete local partition. Multiplying the prime powers in each
position recovers the global integer Smith factors. Primes outside the
pairwise brackets contribute only zero exponents.

## 3. Exhaustive valuation proof

The first and second residual determinantal valuations from Section 1 are

```text
rho=min(mA,(m-1)A+mu,mB),
sigma=min(2mA, mA+mB+C,
                  (m-1)A+mB+v_p(mb-(m+1)a)).        (9)
```

The total valuation is `T=2mA+mB+2C`, and the residual list is
`(rho,sigma-rho,T-sigma)`.

If `A<=B`, the linear numerator has valuation at least `A`. Its candidate
in `(9)` is therefore at least `mA+mB>=2mA`; the value/value candidate is
also at least `2mA`. Hence `sigma=2mA`, and the first ideal gives exactly
`rho=(m-1)A+min(A,mu)`. This proves `(5)`.

If `A>B`, ultrametricity gives `B=C=e`, and write

```text
a=p^(e+d)u,       b=p^e v,       u,v in Z_p^*.
```

The second ideal becomes

```text
sigma=2me+(m-1)d
        +min((m+1)d,e+d,v_p(mv-(m+1)p^d u)).         (10)
```

If `p` does not divide `m`, the last valuation is zero. If `p` divides `m`,
then `m+1` is a unit. When `d!=mu`, the two terms inside that numerator have
distinct valuations, so the valuation is exactly `min(d,mu)`. This is no
larger than either other cap in `(10)`. The first ideal gives `rho` in
Case II, proving `(6)` for all nonresonant cases, including `mu=0`.

When `d=mu`, let `c=m/p^d` and write

```text
v_p(mv-(m+1)p^d u)=d+v_p(cv-(m+1)u).
```

The competing first two minors in `(10)` then cap the extra valuation at
`min(md,e)=K`. Since `v` is a unit, the capped coordinate is
`min(K,v_p(c-(m+1)u/v))`. Moreover the first ideal reduces to
`rho=(m-1)e+K`, and

```text
sigma=2me+md+kappa,       T=(3m+2)e+2md.
```

Their differences give `(8)`. The two caps have distinct mechanisms:
`e` comes from the value/value minor, while `md` comes from the same-node
value/derivative minor. Keeping only the most cancelled linear minor would
lose both competitors and yield a false precision law.

## 4. Intrinsic cap, complete residue ladder, and real failure boundaries

For two unit-separated references `w,w'`, Plucker's identity implies that

```text
det(v_2,w)det(v_1,w') / [det(v_1,w)det(v_2,w')]
 =1 modulo p^e,
```

because the cross difference contains `det(v_1,v_2)` of valuation `e`.
Both normalized cross ratios in `(7)` are units, so
`tau_w=tau_(w') modulo p^e`. Since `m+1` is a unit in Case III and `K<=e`,
their linear expressions have the same valuation capped at `K`. Thus `kappa`
is intrinsically defined. Unit rescalings of the three representatives and
projective `GL_2(Z_p)` changes cancel directly in the exact bracket ratio.

In a chart with `w` at infinity, translate `v_0` to zero. The ratio becomes
`tau=a/(p^d b)=u/v`, so it is precisely the coordinate used by the valuation
proof. No exact rational cross-ratio independence is claimed. It is the
specific capped valuation consumed by the Smith module that is independent.

The entire ladder is attained for each fixed `m,p,e` in Case III:

- if `K=0`, only `kappa=0` occurs;
- at an odd prime and `K>=1`, every integer `0,...,K` occurs;
- at two and `K>=1`, exactly `1,...,K` occur.

Indeed `c=m/p^d` and `m+1` are units. Take `v=m+1`. For `1<=r<K`,
use `u=c+p^r`, giving exact numerator valuation `r`; for `r=K`, use
`u=c`, giving an exact zero numerator. At an odd prime, choose a nonzero
residue `z` such that `c+z` is a unit and take `u=c+z`, giving valuation
zero. Such a choice exists because there are at least three residues. At
two, both units are odd, so valuation zero is impossible; all positive
values and the cap are attained as above.

Thus there are exactly `K+1` different partitions at odd primes, and exactly
`K` at two when `K>=1`. Every different `kappa` gives a different ordered
factor list. At two the case `K=1` is metric-only despite the resonance;
active unit dependence starts at `e>=2`. The old `(2,2,1)` bit is recovered:
`m=2,p=2,d=1,K=min(e,2)` gives two states only once `e>=2`. At odd primes,
unit dependence already starts at `e=1` when `p` divides `m`.

For fixed multiplicity, the ladder saturates at `m*v_p(m)` as `e` grows.
Across multiplicities it is unbounded: choosing `m=p` gives a ladder of
`p+1` states at odd `p` once `e>=p`. This identifies the missing information
and its exact finite residue budget for the entire family.

The raw valuation is not a projective invariant. At `m=3,p=3,e=d=1`, take
`a=9,b=12`; then `3b-4a=0`. The lawful chart `x -> x/(1+x)` changes the
normalized ratio to `tau=13/40`, so `v_3(1-4tau)=1`. The capped value is still
`K=1` and the full Smith module is unchanged. The first invalid implication
would be to preserve an individual infinite cancellation rather than the
capped determinantal ideal. The strongest survivor is precisely `(7)`.

## 5. Projective scope and the multiplicity-label hostile

The incoming
[complete projective higher-jet theorem](continuing1_20260906_jets_projective.md)
proves that the full observer is invariant under integral unit chart changes,
complete local Hasse-coordinate changes, and unit representative changes.
It also proves the determinant

```text
|det E|=|det(v_0,v_1)|^(2m)
         |det(v_0,v_2)|^m |det(v_1,v_2)|^2.          (11)
```

There is a common affine chart with unit denominators unless `p=2` and the
three directions occupy all three projective residue classes. Choose the
pole in an unoccupied class. In the exceptional three-class case `(11)` is
a unit, so every factor is a unit, already included in `(5)`.
Otherwise normalize in that chart, translate the `m`-fold node to zero,
and apply `(1)` and the valuation argument over `Z_p`. These integral ideal
identities work over `Z_p` directly; alternatively integer approximation above
the known determinant precision transfers every determinantal divisor.
This proves that Sections 2--4 classify primitive directions including infinity.

For the first new odd-prime example, fix `m=3,p=3,e=3,d=1`, doubled node
`a=81u`, and simple node `b=108`, with the threefold node at zero. The four
choices `u=2,4,10,1` give respectively

```text
(0,0,0,9,12,18),
(0,0,0,9,13,17),
(0,0,0,9,14,16),
(0,0,0,9,15,15).                                   (12)
```

All have the same multiplicity-labelled bracket depths `(4,3,3)` and the
same determinant valuation 39. These are full six-dimensional observers,
not selected minor statistics.

At `u=1`, swapping which of the positions `0` and `81` carries multiplicity
three changes the last list to the first. The two positions form the same
closest pair, so the associated weighted metric still has the same distances
under that swap. The affine translation taking the new heavy node to zero
gives `a=-81,b=27`; its unit coordinate has `kappa=0`. Thus unequal
multiplicity labels and their position are load-bearing. When `m=2` the
first two multiplicities coincide, and the previous endpoint-swap invariance
is recovered automatically from the common observer; no such symmetry is
asserted for unequal banks.

## 6. Exact full-kernel hierarchy

Within one resonant weighted metric, put

```text
P=(m+1)e+md-K,       Q=(m+2)e+md.
```

The varying factor pair is `(P+kappa,Q-kappa)`. For adjacent attainable
states `kappa` and `kappa+1`, write `H_kappa(N)` for the kernel of the full
observer modulo `p^N`. Then for every integer `N>=1`,

```text
|H_(kappa+1)(N)| / |H_kappa(N)|
 =p  if P+kappa+1 <= N <= Q-kappa-1,
 =1  otherwise.                                    (13)
```

The interval is nonempty whenever both states are attainable, since its
number of integer levels is `e+K-2kappa-1>=e-K+1>=1`.
All other factors agree; summing their capped exponents proves `(13)`.
Larger `kappa` therefore increases or preserves the complete kernel at every
precision while lowering the sharp uniform recovery loss by exactly one
per adjacent step. This is a full-observer consequence because all factors
have been classified, not a claim from an isolated pair in an unknown
partition.

## 7. Computation, connection, and next boundary

The [standalone verifier](../../04-computation/third_20260906_smith.py) and
[matching output](third_20260906_smith.out) compare all nine second minors,
the global gcd formula, the case classification at four primes, independent
DVR elimination, literal homogeneous full Hasse matrices, projective and
tangent changes, every attainable residue state in the stated finite bank,
and exact kernel ratios. The first `m=3` label and chart hostiles are retained.
Finite banks challenge the analytic proof; the universal quantifiers come
from the minor identities and complete valuation cases.

```text
python3 -B 04-computation/third_20260906_smith.py
python3 -B -O 04-computation/third_20260906_smith.py
```

Normal and optimized output agree byte for byte, with **284,133 active gates**:
7,220 integer residual rows and 28,880 prime partitions, 4,889 attained
capped-residue states, 200 literal full projective matrices, 120 reference
changes, and 24 literal full-inverse small-bank value-attainment controls.
The residual bank is `m=2,...,20`, distinct nonzero `a,b` in `[-10,10]`, and
primes `2,3,5,7`. The residue bank uses `m=2,...,48`, the same primes dividing
`m`, and depths `e=0,1,2,3,m*v_p(m),m*v_p(m)+2`, retaining every attainable
state. These finite bounds do not replace the all-`m` analytic proof.

Raw LF-byte SHA-256:

```text
source f5a3deb05f8a5068687d6b050c5824b45a54c10765a5b7cfa709ceb474db2fcc
output b4e08a8d0220742ae8f6692f1644c8b3ecd06cd26555304456b428044c7eaaad
```

The source-to-target map clears a **complete** heavy bank, preserving the
entire integral cokernel up to unit factors. It retains a three-row residual
observer for every multiplicity. Passing from minors to ideals keeps all
invariant factors; selecting one numerator alone loses its competitors.
The projective quotient forgets raw cancellation beyond the needed cap, and
`kappa` restores exactly the coordinate that the three-factor consumer uses.

The next change of residual size is `(m,2,2)` or `(m,3,1)`. Their remaining
rank is four, so the three invariants used here no longer determine the full
partition. The present theorem supplies no classification for those families
and no general LRC or Laurent consequence. A next round must preserve another
intermediate ideal rather than treating increasing multiplicity alone as the
obstruction.
