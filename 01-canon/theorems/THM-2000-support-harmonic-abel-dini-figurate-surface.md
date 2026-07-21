---
id: THM-2000
title: SUPPORT-HARMONIC ABEL--DINI LAW AND THE TWO-AXIS FIGURATE MASS SURFACE
status: PROVED analytic/exact.  A sequence interpreted literally as a subset of the harmonic numbers is governed by its value support, not its indexing multiplicity.  The support/multiset collision tax, Abel--Stieltjes identity, multiplicative-block Dini criterion, full Bertrand near-linear boundary, master-figurate beta-integral surface, simplex/power equal-mass ladder, polygonal digamma axis, three Faulhaber values, Gauss triangular-theta product, ladder ratio/sum/product trichotomy, and reciprocal tournament-series reversal are proved.  Exact and high-precision referees are optimization-safe.  Numerical atlas values carry no unproved irrationality or transcendence claim; global H-spectrum divergence remains open
source: codex-2026-07-21 reciprocal-sequence continuation and audit of concurrent THM-1985/1990
depends_on: [THM-1127, THM-1360, THM-1990]
related: [THM-841, THM-853, THM-874, THM-1370, THM-1985, THM-2005, MISTAKE-209, MISTAKE-210]
external:
  - "Lawrence Downey, Boon W. Ong, and James A. Sellers, Beyond the Basel Problem: Sums of Reciprocals of Figurate Numbers, College Mathematics Journal 39 (2008), 391--394, JSTOR stable 27646686, https://www.jstor.org/stable/27646686"
  - "Archived preprint: https://web.archive.org/web/20130529032918/http://www.math.psu.edu/sellersj/downey_ong_sellers_cmj_preprint.pdf"
  - "Applegate--Pol--Sloane, The Toothpick Sequence and Other Sequences from Cellular Automata, arXiv:1004.3036 (dyadic formula for A139250)"
script: 04-computation/support_harmonic_abel_dini_figurate_surface_thm2000.py
output: 05-knowledge/results/support_harmonic_abel_dini_figurate_surface_thm2000.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/SupportHarmonicFigurate.lean
script_sha256: 27f7fe4ff60a3d78efac586eb11b138dda5dd5262a5177046829e06c032b6502
output_sha256: fbc4930f5861be1f276ea7c1e45944583f2217d2375f7c068ce3c3e53fabcdc2
formalization_sha256: 801a41c4d9bdb9b418427007a69b7b1a75a90ba39ff8f772b4388fe3c25c279f
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
a_n asymptotic n(log n)^alpha
```

gives

```text
sum_n1/a_n converges  iff alpha>1.                                    (11)
```

In particular `n log n` is superlinear and still divergent.  More generally,
with `log_r` the `r`-fold iterated logarithm,

```text
a_n asymptotic
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

For convergence, monotonicity of `t^(-1-epsilon)` gives

```text
x_n/S_n^(1+epsilon)
 <=integral_(S_(n-1))^(S_n)t^(-1-epsilon)dt,
```

and the integrals telescope.  For divergence put `y_n=x_n/S_(n-1)` and
ignore the first term.  If `y_n>=1` infinitely often, then
`x_n/S_n=y_n/(1+y_n)>=1/2` infinitely often.  Otherwise `y_n<1` eventually,
and

```text
log(1+y_n)<=y_n<=2y_n/(1+y_n)=2x_n/S_n.
```

But `sum log(1+y_n)=log(S_n/S_N) -> infinity`.  This proves both halves,
including the exact exponent-one boundary, without importing an asymptotic
growth class.

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

## 5. The Faulhaber axis: three exact reciprocal masses

For the Rosetta/power-sum triangle let

```text
F_p(n)=sum_(k=1)^n k^p.                                                (31)
```

`F_0(n)=n` is the harmonic divergence.  For every `p>=1`, Faulhaber's leading
term makes `F_p(n)` of order `n^(p+1)`, so the support mass converges.  The first
three convergent values are exact:

```text
sum_n1/F_1(n)=2,                                                       (32)
sum_n1/F_2(n)=18-24 log 2,                                            (33)
sum_n1/F_3(n)=4pi^2/3-12.                                             (34)
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

The resulting sum is `4(2zeta(2)-3)`.  Notice that `F_2(n)=N(4,3,n)`:
the square-pyramidal Faulhaber column meets the master figurate surface at
`(s,d)=(4,3)`, and the numerical rows agree.

The referee records high-precision values for `p=4,5,6` but makes no closed-
form or arithmetic-type claim for them.

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
| Fibonacci | exponential, one repeated `1` | standard reciprocal-Fibonacci constant minus `1` = `2.3598856662...` |
| Farey endpoint total `2 sum_(d<=n)phi(d)-1` | asymptotic `(6/pi^2)n^2` | convergent; no toothpick recurrence required |
| toothpick `A139250` | dyadic quadratic lower bound | convergent |
| unlabeled tournament census `A000568` | quadratic-exponential | convergent; arithmetic nature unnamed |
| self-line/orbit sequences | only finite prefixes known in the cited atlas | numeric prefixes only; no all-`n` inference |
| primes | `p_n asymptotic n log n` | divergent Bertrand boundary |

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

The Hamiltonian-value spectrum requires the same honesty.  THM-1370 proves
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

The user's reciprocal transform produces the correct analytic objects.  First
retain the indexing and put

```text
Z_V(s)=sum_n V_n^(-s).                                                 (52)
```

If `sigma=Re(s)>0`, (49) gives

```text
|V_n^(-s)|<= (n!)^sigma 2^[-sigma C(n,2)],                            (53)
```

whose ratio tends to zero.  At `sigma=0` the absolute terms are one; for
`sigma<0` they grow.  Hence this indexed `Z_V` has exact abscissa of absolute
convergence

```text
Re(s)=0.                                                              (54)
```

Similarly,

```text
sum_n z^n/V_n                                                         (55)
```

is entire, because `V_n^(1/n)` tends to infinity.  For the support profile,
`V_n -> infinity` makes the deduplicated support infinite, while it is a
subseries of the indexed profile for `sigma>0`.  Its absolute terms are again
one at `sigma=0`; hence its abscissa is also exactly zero.  This argument does
not assume that all collisions are confined to the initial values.  The
normalized Burnside correction is another legitimate decaying sequence; the
old growing-coefficient series is not.

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

The analytic identities above are proofs, not fits.  The stored output was
reproduced byte-for-byte under ordinary and optimized Python.  The script
independently checks the finite rational identities, collision taxes, Abel
formula, block sandwich, master-array recurrences and monotonicities,
polygonal special values, Faulhaber decompositions, Gauss product, and ladder
triad.  Numerical atlas rows are labelled as partial values where appropriate.

The sorry-free Lean module
`TournamentH7/SupportHarmonicFigurate.lean` certifies the finite algebraic
kernel: master and ordinary-polygonal factorizations, reciprocal
decompositions, the maximum-`c_3` denominator algebra, and finite block
sandwiches.  Infinite summation, beta-integral interchange, and special-value
evaluation remain paper proofs rather than postulated axioms.

No unnamed reciprocal constant is asserted irrational or transcendental.
No finite sequence prefix is extrapolated to an all-`n` law.  In particular,
this theorem does not prove the conjectural global H-spectrum completeness.
