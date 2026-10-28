# Aliquot source classes: term count is not polynomial degree

Date: 2026-10-28. Status: **PROVED** elementary classifications and bounds;
**FINITE-EXACT** independent numerical checks. No literature-priority claim.
The user clarified that “cubic” means **exactly three proper-divisor terms**.

## Results first

Write `s(n)=sum_{d|n, d<n} d`, including the divisor **1** once. Then

* exactly three proper-divisor terms means **n=p³ or n=pq**, with distinct
  primes p,q in the second case;
* exactly five means **n=p⁵ or n=p²q**, with p≠q;
* after giving prime powers their own first class, the user's partition becomes
  **prime powers / distinct-prime products / all other mixed-prime sources**;
  the third class always has at least five proper-divisor terms;
* **20 belongs to the five-term class:** `s(20)=1+2+4+5+10=22`.
  Its additional equality `22=20+1+1` is true but is not a sum of its distinct
  proper divisors. It is the independent condition `s(20)=20+2`;
* among **all** sources with exactly five proper divisors, the abundant ones
  are precisely **12,18,20**, the only perfect one is **28**, and every other
  such source is deficient.

“Three-term” and “five-term” are preferable to cubic and quintic: polynomial
degree is a different invariant. No assertion about fixed-excess solutions of
all other shapes is needed for these conclusions.

## 1. Inheritance, scope and live concepts

Session82 [ALIQUOT_FIBERS.md](https://github.com/eliottcassidy2000/math/blob/6192c0df22ad13989070f4c6e8d3eb0ad7000a3f/05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/ALIQUOT_FIBERS.md)
already supplies the factorized sigma
formula, `s(pq)=1+p+q` for **distinct** primes, finite complete inverse fibers,
and the bound `n≤(t−1)²` for `s(n)=t≥2`. Its
[FIBER_PORTS.md](https://github.com/eliottcassidy2000/math/blob/6192c0df22ad13989070f4c6e8d3eb0ad7000a3f/05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/FIBER_PORTS.md)
emphasizes that one target can hide different named sources. Session81
[ALIQUOT_OBSERVERS.md](https://github.com/eliottcassidy2000/math/blob/8269c6a93b9856f1221c45f57992998af9facc80/05-knowledge/results/session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md)
already uses valuation profiles to repair
such information loss. These are inherited mechanisms, not new claims.

This lane adds a precise classification of the user's requested *term-count*
classes, its low-layer abundance consequence, and hostile controls against
confusing source shape, polynomial degree and evaluated target. The statements
are elementary derivations, not claims of previously unknown number theory.

Live board: divisor lattice / exponent vector / summand count / polynomial
degree / evaluated fiber / abundance. Closest mechanism: multiplicative divisor
sums. Canonical hostile: prime squares cannot be substituted into the distinct
semiprime formula. Corrected near miss: `20+1+1` counts a repeated unit and the
improper divisor 20. Least-used sidecar: the complete exponent multiset.

No code or theorem from another lane is a dependency of the checker here.

## 2. Canonical summands and formal polynomial

For `n≥2`, let `n=product_i p_i^a_i`, with distinct primes and positive integer
exponents. Its proper divisors correspond bijectively to vectors

`0≤b_i≤a_i`, excluding `b=(a_1,...,a_r)`.

Unique factorization makes all the evaluated divisor monomials distinct. Thus

```
P_a(x_1,...,x_r) = product_i(1+x_i+...+x_i^a_i) - product_i x_i^a_i,
s(n) = P_a(p_1,...,p_r).
```

All coefficients of the remaining monomials are one. This is a canonical
expanded divisor-sum representation; it is not a minimum-length arithmetic
circuit and not a claim that alternative numerical decompositions are unique.

With `ω(n)=r`, `Ω(n)=sum a_i`, and `τ(n)=product(a_i+1)`:

| Coordinate | Exact meaning |
|---|---|
| `ω(n)` | number of distinct prime variables |
| `Ω(n)` | degree of the deleted top monomial |
| `D(n)=Ω(n)−1` | total degree of `P_a` |
| `T(n)=τ(n)−1` | number of proper-divisor terms |
| exponent multiset `{a_i}` | full divisor-lattice shape, up to prime relabeling |
| `s(n)` | evaluation at the actual primes |

The degree assertion follows because the unique degree-Ω monomial is deleted
and lowering any positive exponent once leaves a degree-(Ω−1) monomial.
For a prime, `P_(1)=1` has degree zero and one term.

Elementary global comparisons are

`Ω(n)≤T(n)≤2^Ω(n)−1`, and `2^ω(n)−1≤T(n)`.

The first lower bound follows by expanding `product(1+a_i)`; equality holds
exactly for prime powers. The upper bound follows from `a+1≤2^a` for positive
integers a; equality holds exactly when n is squarefree. The last lower bound
uses `a_i+1≥2`. These are inequalities about the canonical expansion, not
operation counts after factorization.

Even `(ω,Ω,τ)` does not recover the exponent shape: multisets `(1,5,5)` and
`(2,2,7)` both give `(ω,Ω,τ)=(3,11,72)`. For example 38880=`2^5*3^5*5`
and 28800=`2^7*3^2*5^2` have 71 proper-divisor terms and polynomial degree10,
but different divisor-lattice shapes. An exhaustive partition check through
Ω=11 shows that Ω=11 is the first such collision **when all three coordinates
are retained**. That minimum is a small FINITE-EXACT certificate, not needed
by any classification proof below.

## 3. Exact term-count classification and the repaired partition

**Theorem.** For each integer `k≥1`, sources with exactly k proper divisors
have precisely the exponent multisets obtained by writing `k+1` as an
unordered product of integers at least2, then subtracting1 from every factor.

**Proof.** `T(n)=k` iff `product(a_i+1)=k+1`. Unique factorization of n gives
the forward implication; assigning distinct primes to any such exponent
multiset gives the converse. The multiplicative partition of k+1 need not be
unique. Prime labels and their assignments remain free. ∎

The first layers are:

| Proper-divisor terms k | Complete source shapes; distinct letters are distinct primes |
|---:|---|
|1|p|
|2|p²|
|3|p³, pq|
|4|p⁴|
|5|p⁵, p²q|
|6|p⁶|
|7|p⁷, p³q, pqr|
|8|p⁸, p²q²|
|9|p⁹, p⁴q|
|10|p¹⁰|
|11|p¹¹, p⁵q, p³q², p²qr|

In particular a mixed-prime source with k terms exists **iff k+1 is
composite**: at least two factors `(a_i+1)` force compositeness; conversely
any `k+1=uv` with u,v≥2 gives `p^(u−1)q^(v−1)`.

For **every n≥2**, the disjoint, exhaustive three-way partition is:

1. **Prime-power class:** `ω(n)=1`, including primes.
2. **Distinct-semiprime class:** `ω(n)=Ω(n)=2`, so `n=pq` with p≠q.
3. **Higher mixed class:** `ω(n)≥2` and `Ω(n)≥3`.

Class1 has priority: its prime cubes also have three terms, and its fifth
powers also have five terms. Without that priority the proposed labels overlap.
For class3, either an exponent is at least2 and another is positive, giving
`τ≥3*2=6`, or at least three exponents are positive, giving `τ≥2³=8`.
Thus class3 has at least five terms. These three classes do not classify
positive n=1, which should be separately recorded as the unit with `s(1)=0`
and no proper divisors. Target1 has exactly all prime sources, and is infinite;
for targets t≥2 every source is composite.

## 4. Requested examples, including the missing unit

| Source | Factorization | Full proper-divisor sum | Terms | Degree | Class |
|---:|---|---|---:|---:|---|
|6|2*3|1+2+3=6|3|1|distinct semiprime|
|10|2*5|1+2+5=8|3|1|distinct semiprime|
|15|3*5|1+3+5=9|3|1|distinct semiprime|
|14|2*7|1+2+7=10|3|1|distinct semiprime|
|12|2²*3|1+2+3+4+6=16|5|2|higher mixed|
|18|2*3²|1+2+3+6+9=21|5|2|higher mixed|
|20|2²*5|1+2+4+5+10=22|5|2|higher mixed|
|28|2²*7|1+2+4+7+14=28|5|2|higher mixed|
|24|2³*3|1+2+3+4+6+8+12=36|7|3|higher mixed|
|36|2²*3²|1+2+3+4+6+9+12+18=55|8|3|higher mixed|

For p²q, the universal formula is `1+p+p²+q+pq`, not a quintic polynomial.
For p³q it is `1+p+p²+p³+q+pq+p²q`. For pqr it is
`1+p+q+r+pq+pr+qr`: both latter types have seven terms but different degrees
(three and two). Equal term count is already a lossy shape summary.

## 5. Even term counts and sharp first mixed witnesses

Divisors pair as d and n/d, with a fixed point exactly when n is a square.
Equivalently `τ(n)` is odd iff all a_i are even. Therefore

**T(n) is even iff n is a square.**

This is a statement about the *number* of summands, not the parity of their
sum. For example the even source36 has eight terms whose sum55 is odd, while
the odd source225 has eight terms whose sum178 is even.

A mixed-prime square has at least two distinct squared primes, so its least
possible value is `2²3²=36`. The least **odd** mixed-prime square is
`3²5²=225`, with

`s(225)=1+3+5+9+15+25+45+75=178`.

Both have eight proper-divisor terms. Thus an interpretation in which all
mixed classes are required to have an *odd* number of summands is false;
“five or more” is valid, “five, seven, nine, ... only” is not. There is no
mixed-prime four-term or six-term case, because5 and7 are prime.

## 6. A complete abundance theorem for the five-term layer

Call n abundant, perfect or deficient as `s(n)>n`, `=n` or `<n`.

**Theorem.** If `T(n)=5`, the abundant sources are exactly12,18,20; the only
perfect source is28; all other sources in this layer are deficient.

**Proof.** Section3 gives n=p⁵ or n=p²q, p≠q. Every prime power is deficient:

`s(p^a)/p^a = 1/p + ... + 1/p^a < 1/(p−1) ≤ 1`.

For n=p²q, subtraction of p²q from its five-term sum gives

`s(n)−n = p²+p+1 − (p²−p−1)q`.

* If p=2, this is `7−q`. Since q is an odd prime, it is positive at q=3,5
  (sources12,20), zero at q=7 (source28), and negative for q≥11.
* If p=3, it is `13−5q`. It is positive at q=2 (source18), while every other
  allowed prime q is at least5 and makes it negative.
* If p≥5, q≥2 gives the upper bound `−p²+3p+3<0`. The expression is already
  −7 at p=5 and strictly decreases thereafter.

These cases exhaust all primes and prove the assertion. ∎

Consequences:20 is the unique five-term source with `s(n)=n+2`; 12 and18
have excesses4 and3, respectively. Similarly, among the three-term layer,
prime cubes are deficient, and

`s(pq)−pq = 2−(p−1)(q−1)`

shows that6 is the only perfect member and every other member is deficient.
Thus neither layer's abundance result is an extrapolation from a finite scan.

## 7. Evaluation loses the source class

The complete fiber `F6={6,25}` contains a distinct semiprime and a prime square.
The complete fiber **`F16={12,26}`** contains a five-term source and a three-term
source: `1+2+3+4+6=1+2+13=16`. In particular the scalar target16 does not
determine the number of canonical summands or the source class.

This is the earliest target where classes2 and3 coexist. To prove minimality,
a class3 source contains either p²q or pqr as a divisor. Proper-divisor sum
is increasing under divisibility, and is strictly increasing for a proper
multiple of a composite. The minimum possible `s(p²q)` is16, uniquely at
p=2,q=3: for p=2 the formula is `7+3q≥16`, and p≥3,q≥2 gives at least21.
The minimum `s(pqr)` for distinct primes is `s(30)=42`, by coordinatewise
monotonicity of `1+p+q+r+pq+pr+qr`. Hence no class3 source has target below16;
at16, source12 is present together with26. Completeness of the displayed
fibers is separately checked over the proved source cutoff `(t−1)²`.

Connection contract: source = integers with labeled prime factorizations;
map = evaluation of their proper-divisor polynomials at those primes;
preserved = their scalar aliquot sum; destroyed = exponent shape, term count
and prime labels; sidecar = prime factorization, or the exponent multiset if
only shape is wanted; decisive hostile =12 and26. No tournament or graph
claim follows merely from this shared loss-of-information pattern.

## 8. Literature boundary and terminology

The divisor-count and multiplicative-sigma formulas are standard; the above
proofs are supplied so this lane does not depend on a citation for them.
The fixed-excess condition `s(n)=n+2` is not quasiperfection, which instead
means `s(n)=n+1`.

For an **even** n with `s(n)=n+2`, deleting the proper divisor2 leaves a sum
equal to n, so it is near-perfect with redundant divisor2 under the definition
in [Pollack–Shevelev, *On perfect and near-perfect numbers* (2012), pp1–2](https://www.pollack-math.net/pnp-4.pdf).
The paper lists20 with that redundant divisor. This observation does not make
the expression `n+1+1` a proper-divisor expansion; the proper deletion is one
divisor2, not two copies of1. “Two-near-perfect” has a different deletion-count
meaning and is not implied by this repeated-unit expression.

The equation `σ(n)−2n=2` itself is studied in
[Li–Ma–Zhang, *Odd solutions of σ(n)−2n=2 have at least six distinct prime factors* (2013)](https://publi.math.unideb.hu/paper/1759/download/).
That paper uses “near-perfect” for this narrower equation, so the terminology
must not silently be mixed. The companion fixed-excess lane handles the
broader equation and its existing literature. No global novelty or complete
fixed-excess classification is asserted here. Sources retrieved2026-10-28.

## 9. Reproduction and hostile controls

Python3 standard library only:

```
python3 check_divisor_term_classes.py > classes_results.json
python3 -O check_divisor_term_classes.py > classes_results.optimized.json
cmp classes_results.json classes_results.optimized.json
```

The default source universe is all n≤100000. A proper-divisor sieve is
independent of factorization-based tau and sigma. Direct divisor-pair
enumeration independently verifies all sources through2000 and every named
hostile. Checks cover the complete three-way partition, both abundance layers,
term-count parity, all requested example sums, the complete small fibers,
and multiplicative partitions through term count31. The exponent-summary
collision minimum is checked by all integer partitions of Ω≤11. Explicit
exceptions keep checks active under optimization. Deterministic JSON output
contains no timestamps or wall-clock measurements. Generated JSON is an audit
receipt, not a mathematical input or required source artifact.

Finite checks certify exactly these universes; global classifications come
from the proofs above. No source-bound scan is used as an all-integer proof.
