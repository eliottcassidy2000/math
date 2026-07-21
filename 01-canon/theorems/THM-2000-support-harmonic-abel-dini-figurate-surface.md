---
id: THM-2000
title: SUPPORT-HARMONIC ABEL--DINI LAW AND THE TWO-AXIS FIGURATE MASS SURFACE
status: PROVED analytic/exact.  A sequence interpreted literally as a subset of the harmonic numbers is governed by its value support, not its indexing multiplicity.  The support/multiset collision tax, Abel--Stieltjes identity, multiplicative-block and partial-sum Dini laws, full Bertrand near-linear boundary and regular-variation tail, Kakeya eventual-overlap/all-strict criteria (including the exact simplex interval decomposition), master-figurate beta-integral and finite digamma surface, divisibility-resonance pi-squared ridge, simplex/power equal-mass ladder, polygonal digamma axis, five Faulhaber values, Gauss triangular-theta product, ladder ratio/sum/product trichotomy, and reciprocal tournament-series reversal are proved.  Exact and high-precision referees are optimization-safe.  Numerical atlas values carry no unproved irrationality or transcendence claim; global H-spectrum divergence remains open
source: codex-2026-07-21 reciprocal-sequence continuation and audit of concurrent THM-1985/1990
depends_on: [THM-462, THM-1127, THM-1360, THM-1990]
related: [THM-841, THM-853, THM-854-noholes-completeness-and-rank2-polygonal-law.md, THM-874, THM-1370-h-spectrum-omits-7-21-all-n.md, THM-1985, THM-2005, THM-2010, THM-2016, MISTAKE-209, MISTAKE-210, MISTAKE-219, MISTAKE-220]
external: "S. Kakeya, On the partial sums of an infinite series, Tohoku Sci. Rep. 3 (1914), 159--164 (achievement-set tail criterion); Applegate--Pol--Sloane, The Toothpick Sequence and Other Sequences from Cellular Automata, arXiv:1004.3036 (dyadic formula for A139250)"
script: 04-computation/support_harmonic_abel_dini_figurate_surface_thm2000.py
output: 05-knowledge/results/support_harmonic_abel_dini_figurate_surface_thm2000.out
lean: 04-computation/lean/TournamentH7/TournamentH7/SupportHarmonicFigurate.lean
script_sha256: 71ea088a8efcb64503d2f223b17fa19cc34216c5e767043efde0348bfb696403
output_sha256: 56717ff8713ea18cf2c9201caca2300810b558af02413204bf5a6fea4e8c9a5a
lean_sha256: c9550106a5c5aa8d8ae18b173b1aa5bf26c552cfa11d615d03f987cb04a1170e
---

# THM-2000 -- support-harmonic Abel--Dini law and the figurate mass surface

## 1. A sequence is a support, not a multiset

Let `a=(a_n)` be a sequence of positive integers and let

```text
nu_a(m)=#{n:a_n=m},                    A=support(a)={m:nu_a(m)>0}.       (1)
```

There are two different reciprocal functionals:

```text
sigma_set(a)  =sum_(m in A)1/m,
sigma_multi(a)=sum_n1/a_n=sum_m nu_a(m)/m.                            (2)
```

The user's instruction -- “the reciprocal of an integer sequence is simply a
subset of the harmonic numbers” -- selects `sigma_set`.  On every finite
truncation one has the exact **collision-tax decomposition**

```text
sigma_multi(a)=sigma_set(a)+sum_m (nu_a(m)-1)_+/m.                    (3)
```

Monotone convergence extends (3) to infinite values in `[0,infinity]`.
Consequently `sigma_set` is invariant under reindexing, finite or infinite
repetition, and deletion of repeated copies.  `sigma_multi` is not.

This repairs several immediate examples.

* `0!,1!,2!,...` has termwise mass `e` but support mass `e-1`.
* Fibonacci and Catalan each repeat the initial value `1`; their support masses
  are their customary termwise reciprocal constants minus one.
* The labeled-tournament values

  ```text
  2^C(n,2):       1,2,8,64,...
  ```

  and the switching values

  ```text
  2^C(n-1,2):     1,1,2,8,64,...
  ```

  have the **same support**.  Their termwise sums differ by one, but their
  subset-harmonic masses are equal.  The old `+1` is exactly the collision tax
  of the repeated `1`.

This distinction matters especially for census sequences, which often begin
with several copies of `1`.  A numerical “constant” must state whether it is a
support value or a term-multiset value.

## 2. Abel--Stieltjes turns the sequence into logarithmic occupancy

For an arbitrary set `A subset N`, define its counting function

```text
A(x)=#{m in A:m<=x}.                                                     (4)
```

Partial summation gives, for every real `x>=1`,

```text
sum_(m in A,m<=x)1/m
  =A(x)/x+integral_1^x A(t)/t^2 dt.                                    (5)
```

This is the faithful Abel transform of a harmonic subset.  It is insensitive
to indexing and keeps exactly the information relevant to reciprocal mass.

Equation (5) gives the exact convergence criterion

```text
sigma_set(A)<infinity
  iff integral_1^infinity A(t)/t^2 dt<infinity.                         (6)
```

Necessity follows from (5), since the integral is the partial mass minus the
nonnegative term `A(x)/x` and is therefore bounded by the partial mass.  For
sufficiency, the only extra point is `A(x)/x ->0`; that follows from the tail
integral because

```text
integral_x^(2x) A(t)/t^2 dt >=A(x)/(2x).                                (7)
```

For a completely discrete Dini form, fix `b>1` and put

```text
B_k=#(A intersect [b^k,b^(k+1))).                                      (8)
```

Every reciprocal in that block lies between `b^(-k-1)` and `b^(-k)`, so

```text
B_k/b^(k+1)
 <=sum_(m in A intersect [b^k,b^(k+1)))1/m
 <=B_k/b^k.                                                            (9)
```

Hence

```text
sigma_set(A)<infinity  iff  sum_k B_k/b^k<infinity.                    (10)
```

Thus reciprocal mass is **summed relative occupancy on logarithmic scale**.
Natural density zero is nowhere near sufficient: the primes occupy only about
`1/k` of the `k`-th exponential block, but those occupancies still have a
divergent harmonic sum.  Triangular numbers occupy exponentially less of each
successive block and converge.

### 2.1 Generating, heat, and Dirichlet transforms of one support measure

Put one unit atom at every occupied integer.  Its ordinary support generator
and heat trace are

```text
G_A(x)=sum_(a in A)x^(a-1),
Theta_A(t)=sum_(a in A)e^(-at).
```

Tonelli gives the exact transform dictionary

```text
sigma(A)=integral_0^1 G_A(x)dx
        =integral_0^infinity Theta_A(t)dt,                            (HT1)
D_A(s)=sum_(a in A)a^(-s)
      =1/Gamma(s) integral_0^infinity t^(s-1)Theta_A(t)dt.           (HT2)
```

Here (HT2) is an identity in the extended nonnegative reals for real `s>0`;
for complex `s` it holds in every common half-plane of absolute convergence.

Thus Abel counting, the support Dirichlet series, and theta/heat identities are
Stieltjes, Mellin, and Laplace views of the same counting measure.  In
particular the full profile is a lattice valuation:

```text
D_(A union B)(s)+D_(A intersect B)(s)=D_A(s)+D_B(s).                 (VAL)
```

For real `s>=0`, (VAL) is understood in `[0,infinity]`; for complex `s`, use
a common half-plane of absolute convergence.

This is the pointwise indicator identity integrated against `m^(-s)`.  It also
shows why one scalar cannot be a complete fingerprint: mass-preserving
Egyptian refinements can move through many different supports on the same
valuation level set (the automatic/Egyptian companion atlas is THM-2005).

If `A(x)~C x^alpha`, Karamata's theorem gives

```text
Theta_A(t)~C Gamma(alpha+1)t^(-alpha),             t downarrow 0.
```

For a pure power law, reciprocal convergence is exactly local integrability of
this heat singularity (`alpha<1`).  The quarter-square trace splits into square
and oblong lattice traces; the lifted `G_1` trace splits into shifted square and
oblong Green functions, explaining the `coth/tanh` constants below.  Gauss's
triangular-number theta product, applied in Section 8 to tournament supports,
is another exact heat-trace evaluation, not an unrelated special-function
coincidence.

## 3. The actual boundary is Bertrand, not “linear versus superlinear”

Let `b_1<b_2<...` be the strictly increasing enumeration of the support.  Then
`b_n=O(n)` forces support divergence, and a lower bound
`b_n >= c n^(1+epsilon)` forces support convergence.  This must be stated for
the support enumeration: an indexed sequence may repeat a sparse value (for
example a power of two) arbitrarily many times without changing its support.
The converse asserted in THM-1990 is false even without repetitions: the
entire slowly varying region between these estimates is the classical
Bertrand boundary.

For example,

```text
b_n asymptotic n(log n)^alpha
```

gives

```text
sum_n1/b_n converges  iff alpha>1.                                    (11)
```

In particular `n log n` is superlinear and still divergent.  More generally,
with `log_r` the `r`-fold iterated logarithm,

```text
b_n asymptotic
 n log n log_2 n ... log_(r-1)n (log_r n)^alpha                       (12)
```

has convergent reciprocal series exactly when `alpha>1`.  Equivalently, for a
general product of iterated-log powers, the first exponent which differs from
one decides the series: greater than one converges, less than one diverges.
This is the integral test after the successive substitutions
`u_1=log x,u_2=log u_1,...`.

The primes, with `p_n asymptotic n log n`, are the canonical zero-density
divergent support.  Any claimed “polynomial versus `#P`” interpretation must
therefore retain (5) or (10); a coarse growth label is not a convergence iff.

### 3.1 The Abel--Dini partial-sum edge

The other Dini law used in THM-1985 also has a short exact proof.  Let
`x_n>0`, put `S_n=sum_(j<=n)x_j`, and suppose `S_n -> infinity`.  Then

```text
sum_n x_n/S_n=infinity,
sum_n x_n/S_n^(1+epsilon)<infinity            for every epsilon>0.   (D)
```

For convergence, start at `n=2` (so `S_(n-1)>0`).  Monotonicity of
`t^(-1-epsilon)` gives

```text
x_n/S_n^(1+epsilon)
 <=integral_(S_(n-1))^(S_n)t^(-1-epsilon)dt,
```

and the integrals telescope from `S_1`; the omitted first term is finite.  For divergence put `y_n=x_n/S_(n-1)` and
ignore the first term.  If `y_n>=1` infinitely often, then
`x_n/S_n=y_n/(1+y_n)>=1/2` infinitely often.  Otherwise `y_n<1` eventually,
and

```text
log(1+y_n)<=y_n<=2y_n/(1+y_n)=2x_n/S_n.
```

But `sum log(1+y_n)=log(S_n/S_N) -> infinity`.  This proves both halves,
including the exact exponent-one boundary, without importing an asymptotic
growth class.

### 3.2 Regular variation gives the tail constant

The Abel viewpoint also quantifies convergence.  If the increasing support
enumeration satisfies

```text
b_n asymptotic c n^p,                    c>0, p>1,
```

then asymptotic inversion and the `p`-series tail give

```text
A(x) asymptotic c^(-1/p)x^(1/p),
sum_(b_n>x)1/b_n asymptotic
  c^(-1/p)/(p-1) x^(1/p-1),
sum_(n>N)1/b_n asymptotic
  1/[c(p-1)] N^(1-p).                                         (RV)
```

This attaches a geometric meaning to several constants in the atlas:

```text
triangular numbers:         value-tail ~sqrt(2/x),
quarter-squares:            value-tail ~2/sqrt(x),
Farey endpoint totals:      value-tail ~(pi/sqrt6)/sqrt(x),
first balance side S_1:     value-tail ~1/(2x^(2/3)),
d-simplex numbers:          value-tail
  ~(d!)^(1/d)/(d-1) x^(1/d-1).
```

At degree two, the coefficient of `x^(-1/2)` in the reciprocal tail is the
same coefficient that appears in `A(x)~C sqrt(x)`.  This Abel self-duality
explains why quarter-square and Farey supports have parallel tails despite
very different arithmetic definitions.

### 3.3 The mass is a needle-width budget

Fix `0<c<=1` and attach to each `a in A` an arc of length `c/a` on the unit
circle.  For arbitrary deterministic placements, the first Borel--Cantelli
lemma on the circle gives

```text
sigma(A)<infinity  =>  almost every point lies in only finitely many arcs.
```

If instead the arc rotations are independent and uniform, then for each fixed
point the hit events are independent with probabilities `c/a`.  The second
Borel--Cantelli lemma and Fubini give

```text
sigma(A)=infinity  =>  almost every point is hit infinitely often,
                       for almost every random placement.            (BC)
```

This is a precise Kakeya-needle analogy, not a deterministic Kakeya theorem.
It says that Abel--Bertrand mass is exactly the random repeated-cover threshold;
arithmetic correlations and direction constraints are the extra difficulty in
LRC-style deterministic covering.  Triangular, quarter-square, and Farey
supports have finite total width, while primes sit on the divergent Bertrand
edge.

### 3.4 Kakeya meets Kakeya: the achievement set of reciprocal sub-supports

There is a second natural object after choosing a convergent support
`A={b_1<b_2<...}`: choose **another subset** of its reciprocal atoms and retain
every possible mass,

```text
E(A)={sum_n epsilon_n/b_n : epsilon_n in {0,1}}.                    (K1)
```

Put `x_n=1/b_n` and `R_n=sum_(j>n)x_j`.  Kakeya's achievement-set criterion
says that `x_n<=R_n` eventually makes `E(A)` a finite union of closed
intervals; if `x_n>R_n` for every `n`, it is a Cantor set.  The first statement
has a direct greedy proof.  On an all-overlap tail, start with any residual
`0<=y_n<=x_n+R_n`; choose `epsilon_n=1` when `y_n>=x_n` and zero otherwise.
In the first case `0<=y_n-x_n<=R_n`; in the second,
`0<=y_n<x_n<=R_n`.  The residual therefore stays inside the next tail and
tends to zero.  Prefixing the finitely many earlier atoms gives a finite union
of intervals.  In the all-strict case, the two cylinders at every level are
separated; the binary coding map is continuous and injective, hence its image
is a Cantor set.

Consequently every regularly varying support `b_n~c n^p`, `p>1`, lies on the
interval side, because (RV) gives

```text
x_n/R_n~(p-1)/n -> 0.                                               (K2)
```

The simplex ladder is completely explicit.  For an integer `k>=2` and
`x_n=1/C(n,k)`, `n>=k`, its exact remainder is

```text
R_n=k/[(k-1)C(n,k-1)],
x_n/R_n=(k-1)/(n-k+1).                                              (K3)
```

Thus precisely the first `k-2` atoms, `n=k,...,2k-3`, dominate their
remainders, while every later atom overlaps its tail.  Hence

```text
E_k = union_(epsilon in {0,1}^{k-2})
 [ sum_(n=k)^(2k-3) epsilon_n/C(n,k),
   sum_(n=k)^(2k-3) epsilon_n/C(n,k)+L_k ],                         (K4)

L_k=k/[(k-1)C(2k-3,k-1)].                                         (K5)
```

The `2^(k-2)` intervals in (K4) are pairwise disjoint.  In particular

```text
E_2=[0,2],
E_3=[0,1/2] union [1,3/2],
E_4=[0,2/15] union [1/5,1/3]
    union [1,17/15] union [6/5,4/3].                               (K6)
```

Now compare the equal-mass power support `{k^j:j>=0}`.  Its atom/tail ratio is
constant:

```text
k^(-j) / sum_(r>j)k^(-r)=k-1.                                     (K7)
```

For `k=2`, its achievement set is again `[0,2]`.  For every `k>=3`, it is the
base-`k` digit set with digits `{0,1}`, a self-similar Cantor set of Hausdorff
dimension `log_k 2`.  Therefore simplex numbers and powers of `k` share the
same scalar mass `k/(k-1)` but, for `k>=3`, have opposite achievement-set
type: the simplex set has finitely many interval components and nonempty
interior, whereas the power set is perfect, nowhere dense, and totally
disconnected.  The scalar quotient forgets even the presence of intervals;
the Dirichlet/block profile remembers why.

## 4. The master figurate array has a reciprocal mass surface

THM-1360 and the S109 figurate atlas use the two-axis master array

```text
N(s,d,m)
 =(s-2)C(m+d-2,d)+C(m+d-2,d-1),       s>=3,d>=2,m>=1.                 (13)
```

The slice `s=3` is the `d`-simplex column; the slice `d=2` is the
`s`-gonal column.  Write

```text
M(s,d)=sum_(m>=1)1/N(s,d,m).                                           (14)
```

Put `a=s-2`.  The binomial ratio

```text
C(m+d-2,d)=(m-1)C(m+d-2,d-1)/d                                      (15)
```

also gives the explicit factorization

```text
N(s,d,m)=[m(m+1)...(m+d-2)/d!] [d+(s-2)(m-1)],
```

which exposes positivity and degree on both axes.  Reciprocating it gives

```text
1/N(s,d,m)
 = [1/C(m+d-2,d-1)] d/[d+a(m-1)].                                    (16)
```

Now

```text
1/C(m+d-2,d-1)=(d-1)B(m,d-1),
d/[d+a(m-1)]=d integral_0^1 u^[d+a(m-1)-1]du.                         (17)
```

Insert the beta integral, sum the geometric series in `t u^a`, and use
Tonelli (all terms are nonnegative).  This proves the exact two-dimensional
representation

```text
M(s,d)=d(d-1) integral_0^1 integral_0^1
       (1-t)^(d-2)u^(d-1)/(1-tu^(s-2)) dt du.                          (18)
```

Equivalently,

```text
M(s,d)=d(d-1) sum_(k>=0) B(k+1,d-1)/[d+(s-2)k].                       (19)
```

The common first term is `N(s,d,1)=1`.  Every later denominator increases
strictly with either `s` or `d`: shape adds a positive simplex layer, while
dimension uses the coning recurrence

```text
N(s,d+1,m)=sum_(i=1)^m N(s,d,i).                                      (20)
```

Thus `M(s,d)` decreases strictly along both axes.  It tends to one along
either axis.  For `d -> infinity`, compare with the simplex slice below; for
`s -> infinity`, use dominated convergence under the same slice.

### 4.1 Simplex numbers and powers have the same mass

At `s=3`,

```text
N(3,d,m)=C(m+d-1,d).                                                   (21)
```

The telescoping identity from THM-1990 gives

```text
M(3,d)=sum_(m>=1)1/C(m+d-1,d)=d/(d-1).                                (22)
```

But also

```text
sum_(r>=0)1/d^r=d/(d-1).                                              (23)
```

Therefore:

> **Equal-mass duality.**  For every `d>=2`, the polynomially sparse
> `d`-simplex numbers and the exponentially sparse powers of `d` have exactly
> the same support-harmonic mass.

Triangular numbers and powers of two are only the first instance.  The whole
ladder is

```text
d=2: 2,  d=3:3/2,  d=4:4/3, ... ->1.                                 (24)
```

Equal total mass does not imply equal logarithmic occupancy: (10) sees the
difference immediately.  The beta transform in (18), not a density match, is
what equates them.

### 4.2 The polygonal axis is a digamma ladder

At `d=2`, the `m`-th `s`-gonal number is

```text
P_s(m)=m[(s-2)m-(s-4)]/2.                                             (25)
```

For `s!=4`, partial fractions give

```text
1/P_s(m)=2/(s-4)[1/(m-(s-4)/(s-2))-1/m].                              (26)
```

Using

```text
sum_(m>=1)[1/(m-r)-1/m]=psi(1)-psi(1-r),                              (27)
```

one obtains

```text
M(s,2)=2/(s-4)[psi(1)-psi(2/(s-2))],             s!=4,               (28)
M(4,2)=zeta(2)=pi^2/6.                                                (29)
```

The first values are

```text
triangular s=3:  2,
square     s=4:  pi^2/6,
pentagonal s=5:  3 log 3-pi/sqrt(3),
hexagonal  s=6:  2 log 2.                                            (30)
```

Termwise shape growth proves that these masses strictly decrease to one.
Thus the rational simplex ladder and analytic polygonal ladder are the two
boundary axes of one positive beta-integral surface.

### 4.3 Every interior point has a finite digamma form

The surface is more rigid than its two boundary axes suggest.  Put

```text
a=s-2,                 beta=d/a-1,
H_j=sum_(k=1)^j1/k,    H_j^(2)=sum_(k=1)^j1/k^2.
```

The factorization above says

```text
1/N(s,d,m)=d!/[a (product_(j=0)^(d-2)(m+j)) (m+beta)].               (PF1)
```

Suppose first that `beta` is not one of `0,1,...,d-2`.  Define

```text
C_beta=d!/[a product_(j=0)^(d-2)(j-beta)],
C_j=(-1)^j d!/[a(beta-j)j!(d-2-j)!].                                (PF2)
```

Residues at the simple poles give

```text
1/N=C_beta/(m+beta)+sum_(j=0)^(d-2)C_j/(m+j),
C_beta+sum_j C_j=0.                                                  (PF3)
```

The cancellation makes the harmonic divergences disappear, yielding the
finite exact form

```text
M(s,d)=-C_beta[psi(d/a)+gamma]-sum_j C_j H_j.                        (PF4)
```

Since `d/a` is positive rational, recurrence followed by Gauss's digamma
formula reduces (PF4) algorithmically to a finite real cyclotomic expression:
a rational, an algebraic multiple of `pi`, and an algebraic-linear combination
of logarithms of positive algebraic sine values.  This is an exact description,
not a blanket irrationality or transcendence claim; cancellations can occur.

There is a sharper **divisibility resonance**.  A pole collision occurs exactly
when

```text
a>=2 and a divides d,              r=d/a-1 in {0,...,d-2}.           (PF5)
```

Writing

```text
B=(-1)^r d!/[a r!(d-2-r)!],
A_r=B(H_r-H_(d-2-r)),
A_j=(-1)^j d!/[a(r-j)j!(d-2-j)!]       (j!=r),                       (PF6)
```

the two simple poles merge into

```text
1/N=B/(m+r)^2+sum_j A_j/(m+j),       sum_j A_j=0,
M(s,d)=B[pi^2/6-H_r^(2)]-sum_j A_jH_j.                              (PF7)
```

Thus every resonant point is a rational plus the nonzero rational
`B*pi^2/6`, and is therefore transcendental.  The first exact rows expose the
ridge:

```text
M(4,2)=pi^2/6,
M(4,3)=18-24log2,
M(4,4)=21-2pi^2,
M(4,6)=15pi^2-1175/8,
M(5,3)=pi^2/3-2,
M(5,6)=1205/18-(20/3)pi^2.                                         (PF8)
```

The extra shape pole is linear, so it can collide with at most one of the
consecutive rising-factor poles; no triple or higher resonance is possible.
The square-pyramidal identity `F_2=M(4,3)=18-24log2` lies immediately off the
even-`d` resonance ray of `s=4`, explaining why it carries `log2` while its
neighbors `M(4,2),M(4,4),M(4,6),...` carry `pi^2`.

The simplex axis `a=1` is deliberately just outside the collision set:
`beta=d-1`, and its digamma expression cancels all the way down to the rational
`d/(d-1)`.  In this sense the surface has a rational boundary and a family of
transcendental `pi^2` ridges indexed by the divisor relation `s-2 | d`.

## 5. The Faulhaber axis: five exact reciprocal masses

For the Rosetta/power-sum triangle let

```text
F_p(n)=sum_(k=1)^n k^p.                                                (31)
```

`F_0(n)=n` is the harmonic divergence.  For every `p>=1`, Faulhaber's leading
term makes `F_p(n)` of order `n^(p+1)`, so the support mass converges.  The first
five convergent values have exact forms:

```text
sum_n1/F_1(n)=2,                                                       (32)
sum_n1/F_2(n)=18-24 log 2,                                            (33)
sum_n1/F_3(n)=4pi^2/3-12.                                             (34)
sum_n1/F_4(n)=(-270+480log2)/7
 -(90/7)[D(1-r_4)+D(2+r_4)],
   r_4=(sqrt21-3)/6,             D(x)=psi(x)+gamma,                  (F4)
sum_n1/F_5(n)=60-4pi^2
 -8sqrt3 pi cot(pi(sqrt3-1)/2).                                     (F5)
```

Equation (32) is triangular telescoping.  For (33), use

```text
1/F_2(n)=6/[n(n+1)(2n+1)]
        =6/n+6/(n+1)-24/(2n+1)                                      (35)
```

and take the cancelling harmonic limit.  For (34), Faulhaber gives
`F_3(n)=T_n^2`, while

```text
1/[n^2(n+1)^2]=1/n^2+1/(n+1)^2-2/[n(n+1)].                           (36)
```

The resulting sum is `4(2zeta(2)-3)`.  For `p=4`, the fully rational first
split is

```text
1/F_4(n)=270(2n+1)/[7(3n^2+3n-1)]
         +480/[7(2n+1)]-30/n-30/(n+1).
```

Factoring the remaining quadratic with `3r_4(r_4+1)=1` gives the digamma pair
in (F4).  For `p=5`,

```text
F_5(n)=n^2(n+1)^2(2n^2+2n-1)/12,
1/F_5(n)=48/(2n^2+2n-1)-12/n^2-12/(n+1)^2.
```

Putting `r=(sqrt3-1)/2`, so `r(r+1)=1/2`, splits the quadratic term into
`8sqrt3[1/(n-r)-1/(n+1+r)]`.  Digamma reflection and recurrence then give
(F5).  The `p=4` residues add rather than subtract across its quadratic roots,
so an algebraic-argument digamma pair remains; this is still an exact form,
not a finite-prefix label.

Notice that `F_2(n)=N(4,3,n)`:
the square-pyramidal Faulhaber column meets the master figurate surface at
`(s,d)=(4,3)`, and the numerical rows agree.

The referee records a high-precision value for `p=6` but makes no closed-form
or arithmetic-type claim for it.  Equations (F4)--(F5) likewise assert exact
forms only; no blanket arithmetic classification is inferred from them.

The first triangular balance tower supplies a second semantic test.  At row
`n>=1`, center `C=n(n+1)` makes the two interval sides

```text
{C-n,...,C}            and            {C+1,...,C+n}
```

have the same value

```text
S_1(n)=n(n+1)(2n+1)/2=3F_2(n)=3,15,42,90,... .                     (BAL1)
```

Since `S_1(n+1)-S_1(n)=3(n+1)^2`, these common values are distinct.  Partial
fractions give

```text
1/S_1(n)=2/n+2/(n+1)-8/(2n+1),
sum_(n>=1)1/S_1(n)=6-8log2.                                         (BAL2)
```

This is exactly one third of the square-pyramidal mass (33).  Under support
semantics the equal left/right value is counted once.  Counting the two labeled
sides would double (BAL2), while taking the reciprocal of the whole balanced
row `2S_1(n)` would halve it.  The carrier choice is part of the theorem, not a
cosmetic indexing convention.

### 5.1 The true tournament `c_3` maximum is a parity-spliced Faulhaber column

The earlier reciprocal atlases called `C(n,3)` the maximum number of cyclic
triples.  It is only the number of triple slots.  The true maximum is

```text
M_3(n)=n(n^2-1)/24,                    n odd,
M_3(n)=n(n^2-4)/24,                    n even.                       (C3-1)
```

At odd order `n=2k+1`, this is exactly the square-pyramidal number

```text
M_3(2k+1)=k(k+1)(2k+1)/6=F_2(k).
```

At even order `n=2k`, it is `k(k^2-1)/3`.  Moreover the full interlaced
sequence is strictly increasing and its increments repeat triangular numbers:

```text
M_3(n+1)-M_3(n)=T_floor(n/2),                  n>=1,                 (C3-2)
```

where `M_3(1)=M_3(2)=0`; the two low-order cases are direct.

Thus support and indexed semantics coincide here.  The odd and even parts give

```text
sum_(k>=1)1/M_3(2k+1)=18-24log2,
sum_(k>=2)1/M_3(2k)=3/4,
sum_(n>=3)1/M_3(n)=75/4-24log2.                                (C3-3)
```

The odd-order identity is literally the `p=2` Faulhaber mass; the even part
telescopes from `3/[k(k^2-1)]`.  This replaces the false `3/2` tournament claim
while preserving `3/2` as the correct reciprocal mass of all triple slots.

## 6. Ratios, partial sums, and partial products occupy three regimes

The old Walsh dimension-ladder ratios are `2k-1`.  Starting at `k=2`, apply
three accumulation functors:

```text
ratio:          r_k=2k-1,
partial sum:    s_k=sum_(j=2)^k(2j-1)=k^2-1,
partial product:p_k=product_(j=2)^k(2j-1)=(2k-1)!!.                   (37)
```

Their support-harmonic behavior is completely different:

```text
sum_(k>=2)1/r_k       =infinity,                                     (38)
sum_(k>=2)1/s_k       =3/4,                                          (39)
sum_(k>=2)1/p_k
 =e^(1/2)sqrt(pi/2) erf(1/sqrt(2))-1.                                (40)
```

For (39),

```text
1/(k^2-1)=1/2[1/(k-1)-1/(k+1)].                                     (41)
```

For (40), put

```text
F(x)=sum_(k>=0)x^k/(2k-1)!!.
```

Coefficient comparison gives

```text
2xF'(x)=(1+x)F(x)-1,                 F(0)=1,
F(x)=1+e^(x/2)sqrt(pi x/2) erf(sqrt(x/2)).
```

Evaluating at `x=1` and deleting the `k=0,1` terms proves (40); the referee
also checks the value to 45 digits.  The same local ladder therefore produces a divergent harmonic
edge, a rational telescoping mass, and an error-function mass depending only
on whether one retains ratios, integrates them additively, or integrates them
multiplicatively.  “The operation, not the count, is the self-similar object”
becomes an exact analytic statement here.

## 7. Sequence-atlas consequences

The support criterion sorts several established repo families without forcing
spurious closed forms.

| family | support growth / exact fact | reciprocal status |
|---|---|---|
| Moser `A000127` rows | quartic | convergent; numeric support mass `2.0174822491...` |
| polygonal diagonal `G` | quartic quasi-polynomial | convergent; support mass `2.3827162754...` |
| quarter-squares `A002620` | union of squares and oblongs | exact support mass `zeta(2)+1` |
| Fibonacci | exponential, one repeated `1` | standard reciprocal-Fibonacci constant minus `1` = `2.3598856662...` |
| Farey endpoint total `2 sum_(d<=n)phi(d)-1` | asymptotic `(6/pi^2)n^2` | convergent; no toothpick recurrence required |
| toothpick `A139250` | dyadic quadratic lower bound | convergent |
| unlabeled tournament census `A000568` | quadratic-exponential | convergent; arithmetic nature unnamed |
| self-line/orbit sequences | only finite prefixes known in the cited atlas | numeric prefixes only; no all-`n` inference |
| primes | `p_n asymptotic n log n` | divergent Bertrand boundary |

The newly pulled THM-2010 invariant sequences are, at present, exact **finite
prefixes** for `n=3,...,6`, not asymptotic sequences with established
reciprocal constants.  Their collision taxes already matter: the prefixes
`|R|=(2,2,6,8)` and `|disc|=(1,2,2,5)` each have tax `1/2`, while
`max|R|=(3,3,15,15)` has tax `1/3+1/15=2/5`.  The prefixes for `|specA|`,
`|H|`, cycle-vector count, arborescence-invariant count, maximum total
arborescences, metagraph edges, and WL colors have no collisions in that
range.  The referee records every exact prefix mass, but no four-term trend is
promoted to convergence, density, or arithmetic type.

The quarter-square row is a useful warning against reading only the leading
term.  Its positive support splits disjointly as

```text
{floor(n^2/4):n>=2}={k^2:k>=1} disjoint_union {k(k+1):k>=1}.
```

No positive square equals an oblong number, because `k(k+1)` lies strictly
between the consecutive squares `k^2` and `(k+1)^2`.  Therefore

```text
sum_(m in support(A002620))1/m
 =sum_(k>=1)1/k^2+sum_(k>=1)1/[k(k+1)]
 =zeta(2)+1.                                                         (QS)
```

For the Farey endpoint totals, the classical estimate
`sum_(d<=n)phi(d)=3n^2/pi^2+O(n log n)` sharpens mere convergence to

```text
sum_(k>N)1/[2sum_(d<=k)phi(d)-1]
 =pi^2/(6N)+O(log N/N^2).                                           (F)
```

So their reciprocal constant is approached at a controlled `1/N` rate; a
finite prefix should never be presented as the constant itself.

### 7.1 Polynomial and quasipolynomial sequences have finite polygamma forms

There is a general exact closure principle behind the numerical Moser and `G`
rows.  Let `P in Q[x]` have degree at least two and no zero on `n>=n_0`.
Over its algebraic roots, write

```text
1/P(z)=sum_rho sum_(k=1)^(m_rho)c_(rho,k)/(z-rho)^k.                 (POLY1)
```

The simple-pole residues sum to zero because `1/P(z)=O(z^-2)`.  Therefore

```text
sum_(n>=n_0)1/P(n)
 =-sum_rho c_(rho,1)psi(n_0-rho)
  +sum_rho sum_(k>=2)c_(rho,k)(-1)^k
     psi^(k-1)(n_0-rho)/(k-1)!.                                    (POLY2)
```

This follows by cancelling the simple harmonic divergences and using the
Hurwitz-zeta/polygamma identity for repeated poles.  Splitting into residue
classes proves the analogous statement for every rational quasipolynomial.
For support rather than term-multiset mass, one must still remove collisions;
an eventually injective sequence requires only its finite collision tax.

Consequently `A000127` has a finite quartic-root digamma form.  The parity
quasipolynomial `G` has two such quartic forms, with its repeated initial `1`
removed.  Calling their displayed decimals “numeric support masses” describes
the convenient evaluation, not the exact analytic class; no arithmetic type is
being claimed.

### 7.2 Run support becomes a monotone reciprocal-mass filtration

The S109 at-most-`j`-runs filtration supplies a more structural comparison.
Its row sums are

```text
R_j(r)=sum_(i=0)^j C(r+1,2i),                    r>=0,               (RUN1)
```

and increase pointwise to `2^r`.  The first row has the exact hyperbolic mass

```text
R_1(r)=1+C(r+1,2),
sum_(r>=0)1/R_1(r)=2pi/sqrt7 tanh(pi sqrt7/2).                       (RUN2)
```

`R_2(r)=A000127(r+1)`.  Since `R_j>=R_1` for `j>=1`, dominated convergence
turns the combinatorial inclusion of run modes into the mass chain

```text
2pi/sqrt7 tanh(pi sqrt7/2)
 > sigma(A000127)=2.0174822491... > ... downarrow 2
 =sum_(r>=0)1/2^r.                                                   (RUN3)
```

The diagonal sums `G_j(n)` count tilings with at most `j` domino runs.  Here

```text
G_1(n)=1+floor(n^2/4),
support(G_1)={1} union {k^2+1:k>=1} union {k(k+1)+1:k>=1},
sigma(G_1)=pi/2 coth(pi)+pi/sqrt3 tanh(pi sqrt3/2)-1/2
          =2.8748213280... .                                        (RUN4)
```

The initial duplicate `1` is collapsed.  The branches are disjoint because
the underlying quarter-square sequence is strictly increasing from `n=2`.
Again `G_j>=G_1`, and appending a square gives strict increase after the shared
initial collision, so dominated convergence yields

```text
sigma(G_1)>sigma(G_2=G)=2.3827162754...>...
 downarrow sigma(Fibonacci support)=2.3598856662... .                 (RUN5)
```

Thus filling missing Vandermonde/run modes makes the counts larger and their
reciprocal masses smaller.  Brown completeness, support semantics, and Abel
occupancy all see the same filtration.

There is a particularly useful opposite LRC example.  THM-1127 proves for its
literal toothpick count `N(K)` that

```text
N(K+194040)=N(K)+121304,                 K>9.                         (42)
```

Choosing any phase representative `K_0>9` and iterating (42) shows that the support contains the infinite
arithmetic progression

```text
N(K_0)+m*121304,                 m>=0.                                (43)
```

The reciprocal sum over (43) diverges.  Hence this LRC toothpick count has
divergent support mass even though the geometric toothpick population and the
Farey endpoint total are polynomially sparse enough to converge.  “Toothpick
self-similarity” is not one analytic universality class; its affine phase law
decides the harmonic behavior.

For completeness, the convergence row for the classical toothpick population
does not rely on a numerical fit.  The Applegate--Pol--Sloane dyadic formula

```text
T(2^j)=(2^(2j+1)+1)/3
```

and monotonicity imply, whenever `2^j<=n<2^(j+1)`, that
`T(n)>=T(2^j)>n^2/6`.  Hence `sum_n1/T(n)` is bounded by a constant multiple
of `sum_n1/n^2`.

The Hamiltonian-value spectrum requires the same honesty.  The later file
`THM-1370-h-spectrum-omits-7-21-all-n.md` proves
that `7,21` never occur and that every other odd value through `609` does occur.
It labels

```text
support(H)=odd positive integers minus {7,21}                          (44)
```

as a **conjecture**, not a theorem.  Therefore positive density and divergence
of the global H-value support remain open.  Removing two values from the odds
would preserve divergence, but that conclusion is conditional on (44).

## 8. Gauss closes the tournament triangular theta

Let `q` satisfy `|q|<1`.  The common support of labeled-tournament and
switching values gives

```text
Psi(q)=sum_(r>=0)q^[r(r+1)/2].                                        (45)
```

Gauss's triangular-number product identity is

```text
Psi(q)=(q^2;q^2)_infinity^2/(q;q)_infinity
      =(q;q)_infinity(-q;q)_infinity^2.                              (46)
```

For `0<q<1`, using the positive real branches, this is equivalently

```text
Psi(q)=theta_2(0,sqrt(q))/(2q^(1/8)).
```

At `q=1/2`,

```text
Psi(1/2)=1.641632560655153866293843...                                (47)
```

This is exactly the labeled-tournament support-harmonic mass.  The switching
**termwise** sum is one larger because it repeats `1`; its support mass is
again (47).  Thus the “unnamed partial theta” target in THM-1990 is already a
classical product/theta value.  The signed pentagonal Euler product is a
different, complementary specialization.

## 9. Reciprocalization reverses the tournament analytic failure

Let

```text
L_n=2^C(n,2),                    V_n=A000568(n).                       (48)
```

Orbit counting gives the elementary bounds

```text
L_n/n! <=V_n<=L_n.                                                     (49)
```

They immediately correct several historical analytic-series claims.

1. For every fixed complex `s`, the terms `V_n/n^s` do not tend to zero.
   Therefore the index-Dirichlet series

   ```text
   sum_n V_n/n^s                                                       (50)
   ```

   diverges everywhere.

2. Consecutive absolute terms of the labeled EGF satisfy

   ```text
   |L_(n+1)z^(n+1)/(n+1)!| / |L_n z^n/n!|
     =2^n|z|/(n+1).                                                    (51)
   ```

   For every `z!=0` this tends to infinity.  The series has radius zero,
   not positive radius or an entire continuation.

The user's reciprocal transform produces the correct analytic objects.  Put

```text
Z_V(s)=sum_n V_n^(-s).                                                 (52)
```

If `sigma=Re(s)>0`, (49) gives

```text
|V_n^(-s)|<= (n!)^sigma 2^[-sigma C(n,2)],                            (53)
```

whose ratio tends to zero.  At `sigma=0` the absolute terms are one; for
`sigma<0` they grow.  Hence `Z_V` has exact abscissa of absolute convergence

```text
Re(s)=0.                                                              (54)
```

Similarly,

```text
sum_n z^n/V_n                                                         (55)
```

is entire, because `V_n^(1/n)` tends to infinity.  Support-deduplication only
changes the finitely repeated initial values and does not affect (54)--(55).
The normalized Burnside correction is another legitimate decaying sequence;
the old growing-coefficient series is not.

## 10. Carrier and Tournament Analysis audit

The faithful vertices are occupied integer values, or equivalently the
multiplicative blocks (8), not sequence indices.  The support quotient
preserves membership, counting function, reciprocal mass, and convergence; a
collision-tax sidecar remembers what it deliberately destroys.

Several carriers were challenged: indices, values, logarithmic blocks,
figurate-axis cells, Dirichlet exponents, recurrence states, and proof
obligations.  Ordering finitely many families by computed support mass gives a
transitive tournament (no cycles, singleton SCCs, one Hamiltonian path after
exact-value/name tie-breaking), but it is only atlas telemetry.  A finite mass
order cannot determine tail convergence, equality of closed forms, or
arithmetic type.  The proof carrier is the counting function plus its block
occupancies.

## 11. Scope and reproducibility

The analytic identities above are proofs, not fits.  The script independently
checks the finite rational identities, collision taxes, Abel formula, block
sandwich, master-array recurrences and monotonicities, polygonal special
values, Faulhaber decompositions, Gauss product, and ladder triad.  Numerical
atlas rows are labelled as partial values where appropriate.

No unnamed reciprocal constant is asserted irrational or transcendental.
No finite sequence prefix is extrapolated to an all-`n` law.  In particular,
this theorem does not prove the conjectural global H-spectrum completeness.
