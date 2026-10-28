# Complete aliquot fibers by addition shape, through target 2000

**Status:** elementary identities and completeness bounds **PROVED**; the
specified census and minimal witnesses **FINITE-EXACT**, independently checked.
No novelty, density, Goldbach, or infinite-family claim is made.

## Results first

Let `s(n)=sigma(n)-n` and `F_t={n>=1:s(n)=t}`. Partition sources, in priority
order, into **I:** prime powers; **II:** products of two distinct primes;
**III:** all other composites. Source 1 is separate. In the finite domain
`2<=t<=2000`, there are **28,770 source-target incidences**, and every source is
included—not merely sources below 2000.

| Class | Incidences | Singleton full fibers | Square-sigma incidences |
|---|---:|---:|---:|
| I | 332 | 148 | 2 |
| II | 26,647 | 87 | 283 |
| III | 1,791 | 170 | 74 |

There are 196 empty fibers, 405 singleton fibers, and 34 targets having sources
in all three classes. The maximum fiber size is 92, attained only at target
1891 in this window: 91 class-II sources and the class-III source 5885.
Class I has two sources only at target 31; all other nonempty I subfibers have
one source. The largest III subfiber has six sources, at targets 1483 and 1585.
These multiplicity statements are restricted to the stated target window.

The following minima are minima over **all integers t>=2**, not merely
unqualified search observations: each candidate exists, and complete fibers
for every smaller target were checked.

| Property | Least target | Complete fiber, partitioned |
|---|---:|---|
| II and III both occur | 16 | II: 26; III: 12 |
| More than one II source | 17 | II: 39, 55 |
| Nonempty fiber containing only III | 28 | III: 28 |
| More than one I source | 31 | I: 32, 125; II: 161, 209, 221 |
| I and III both occur | 40 | I: 81; II: 74; III: 44 |
| All three classes occur | 40 | I: 81; II: 74; III: 44 |
| More than one III source | 64 | II: 122; III: 56, 76 |

In particular, 28 is both the least only-III target and a singleton full fiber;
this is a concrete fact, not an implication from its being perfect.

Full presence-mask, exponent-shape, proper-divisor-term, omega, and parity
tables are in [CENSUS_OUTPUT.md](CENSUS_OUTPUT.md). Here omega means the number
of **distinct** prime divisors, not the prime factors counted with multiplicity.

## 1. What the addition rule records

For `n=product p_i^a_i`, the number of proper-divisor terms in the defining
addition is

`k(n)=tau(n)-1=product(a_i+1)-1`.

This counts the actual proper divisors once each; arbitrary regrouping of a
sum does not change the source class. Exactly three proper-divisor terms
means `tau(n)=4`, hence precisely `n=p^3` or `n=pq` with distinct primes.
Priority assigns the cubes to I and the distinct-semiprime case to II.
The third class has at least five terms: two terms force a prime square,
three terms have just been classified, and four terms force a prime fourth
power. Examples exposing coordinate loss are:

* Three terms: `8` lies in I, whereas `6` lies in II.
* Five terms: `32` lies in I, whereas `12` lies in III.
* Two distinct prime divisors: `6` lies in II, whereas `12` lies in III.

Thus the term count or omega alone does not recover the class. The full
exponent pattern does. The census has 14 cube sources in I and 26,647
distinct-semiprime sources in II, hence 26,661 incidences with exactly three
proper-divisor terms. Class-III shapes start with 486 incidences of `(2,1)`
(five terms); they are not an undifferentiated residual.

The minimal class-III target 16 also has a short structural proof. Every III
source has a divisor of shape `p^2 q` or `pqr`, with distinct primes. The
aliquot sum is nondecreasing under divisibility. For `p^2q`, the smallest
possibility when `p=2` is 12 with sum 16; when `p>=3`, the formula
`s(p^2q)=1+p+p^2+(1+p)q` is at least 21. A `pqr` divisor is at least the
three-prime pattern `(2,3,5)`, whose aliquot sum is 42. Thus III cannot occur
below target 16, and 12 realizes equality.

## 2. Why these are complete fibers

For `t>=2`, every source is composite. If p is its least prime divisor,
`n/p` is a proper divisor at least `sqrt(n)`, distinct from 1. Therefore

`t >= 1+n/p >= 1+sqrt(n)`, so **`n <= (t-1)^2`**.

Except for a prime square, p and `n/p` are distinct proper divisors, giving
the sharper **`4n < (t-1)^2`**. This is why excluding large sources without
retaining the prime-square branch would be wrong.

The first implementation is the Session82 increasing-prime-prefix inverse
algorithm, adapted here to a standard-library primality table. At a prefix
`m` with divisor sum S, appending a new largest prime q to exponent one gives

`s(mq)=s(m)q+S`.

If `s(m)>0`, the sole candidate is `q=(t-S)/s(m)`, subject to integrality,
primality, and strict prime ordering. Higher exponents are scanned until
their increasing values exceed t. The two lower bounds used to stop the
next-prime loop are `s(m)p^2+S(p+1)` for a terminal exponent at least two
and `s(m)p^2+2Sp+S` for a nonterminal prime. A recursive extension is retained
only if the next possible prime can still reach at most t. Monotonicity under
divisibility justifies the pruning. Every source has one increasing-prime
factorization path. Every prime divisor of a composite source is a proper
divisor, so it is at most `t-1`; the finite prime table is sufficient.

The second implementation is deliberately different: a direct divisor-addition
sieve computes `s(n)` for **every** `1<=n<=1999^2=3,996,001` and groups those
whose sums lie in 2..2000. Its loops use no inverse-factor pruning or shape
assumptions. Every one of the 1,999 full lists agrees exactly. Each source is
then re-encoded using factorized sigma, and every shape/parity/edge identity
is checked. The output stores only compact tables and the full-list digest,
not a large raw fiber file.

Boundary cases are not silently forced into this census: `F_0={1}`, while
`F_1` is the infinite set of primes. All class-I sources here therefore have
exponent at least two. Inherited hostile controls include `F_2=F_5=empty`,
`F_3={4}`, `F_7={8}`, `F_12={121}`, and `F_13={27,35}`.

## 3. Parity and the prime-gap sidecar

Factoring sigma proves that it is odd exactly when n is a square or twice a
square. Consequently:

* An odd source has even target exactly when it is a square.
* An even source has odd target exactly when it is a square or twice a square.
* A class-I source of an even target is `p^(2a)` for an odd prime p. The
  canonical square `(t-1)^2` occurs exactly when `t-1` is prime; higher even
  exponents, such as `81->40`, require separate handling.
* For an even target, class II is exactly the singleton candidate
  **`2(t-3)` when `t-3` is a prime greater than 2**, and empty otherwise.
  Thus even targets never have two distinct-semiprime sources.

For class II, with primes `p<q`,

`t=1+p+q`, and `(t-1)^2-4pq=(q-p)^2`.

For odd t the two primes are odd, and the positive half-gap r gives
`pq=h^2-r^2`, where `h=(t-1)/2` and `h-r,h+r` must both be prime.
These need not be consecutive primes. Repeated primes are forbidden:
`s(p^2)=p+1`, not `2p+1`. In particular allowing the diagonal `(2,2)` would
falsely manufacture a source for target 5.

There are no only-III odd targets in this finite window. That does not prove
the statement for all odd targets, a prime-untouchability assertion, or a
distinct-prime Goldbach statement. Similarly, the observed relative class
counts imply no universal existence or density result.

## 4. Square-sigma edges by class

An aliquot step `n->t` is a square-sum edge exactly when `n+t=sigma(n)` is a
square, provided the vertices are positive and distinct. Positivity is automatic
here. Distinctness is also automatic for square-sigma sources: if `n=t` and
`2n=z^2`, then n is twice a square, making sigma(n) odd by the parity criterion,
contrary to `sigma(n)=2n`.

In the complete target window, class I contributes only `81->40` and
`343->57`. This statement is finite-exact here; the earlier Session82 report's
global prime-power classification has a separately stated classical dependency
and is not inferred from this enumeration. The least class-II square-edge
target is 14, via `22->14`; the least class-III square-edge target is 74, via
`70->74`. Also `F_78={66}` is a class-III singleton square-edge fiber.

There is an elementary exact parameterization of the even-target II case:

**`t=12k^2+2`, `n=24k^2-2`, with `k>=1` and `12k^2-1` prime.**

Indeed `n=2q` with q an odd prime, and `sigma(n)=3(q+1)=z^2`. Then `z=3u`,
so `q=3u^2-1`. Since q is odd, u is even; write `u=2k`. Conversely these
formulas give `sigma(n)=36k^2`. Examples are `22->14`, `94->50`, and
`214->110`. No infinitude of primes represented by this quadratic is claimed.

For odd-target II sources, square sigma is equivalent to
`r^2+z^2=((t+1)/2)^2` in the centered prime-pair coordinates above. In particular
the hypotenuse `(t+1)/2` must have a prime divisor congruent to 1 modulo 4;
the primality of both centered endpoints remains a separate guard. The checker
tests the Pythagorean identity and this obstruction on every square-edge II
incidence. The obstruction is necessary, not sufficient.
For the obstruction, reduce a nondegenerate integer right triangle to a
primitive one. Its odd hypotenuse exceeds 1; any prime dividing it is coprime
to both legs, and their squared ratio is -1 modulo that prime, forcing the
prime to be 1 modulo 4. This prime also divides the original hypotenuse.

**A small hostile to count-based implication:** the complete fibers
`F_6={6,25}` and `F_14={22,169}` each contain one prime square and one distinct
semiprime. Their class-count, exponent-shape, term-count, omega and source-parity
profiles coincide. Nevertheless neither source of `F_6` has square sigma,
whereas `sigma(22)=36` is square. The missing coordinate is the actual
divisor-sum squarefree signature, not another inference from class counts.

## 5. Reproduction and stopping boundary

Python 3, standard library only; explicit failures remain active under `-O`.
Run from this directory:

```sh
python3 check_shape_census.py > CENSUS_OUTPUT.md
python3 -O check_shape_census.py > /tmp/census-optimized.md
cmp CENSUS_OUTPUT.md /tmp/census-optimized.md
```

Both full runs passed and produced identical bytes. The canonical digest is
SHA-256 of ASCII lines `t:n1,n2,...\n`, for increasing targets 2..2000 and
increasing sources within each fiber:

`01c598ee1831b6c1b255950803fd0271e66934fc408a7cec6c752663d83adfe9`.

Inheritance: Session82 [check_aliquot_fibers.py](../session82_decimal_aliquot_square_sum_20261028/check_aliquot_fibers.py),
[ALIQUOT_FIBERS.md](../session82_decimal_aliquot_square_sum_20261028/ALIQUOT_FIBERS.md), and
[SQUARE_SIGMA_EDGES.md](../session82_decimal_aliquot_square_sum_20261028/SQUARE_SIGMA_EDGES.md) supply
the audited inversion mechanism, bounds, and guarded gap/square connections.
This lane adds a fully source-sieved finite shape census and its independently
certified minima; it does not claim a new inversion algorithm. The frozen
universe ends at target 2000. No conclusion about larger fibers, global class
frequencies, orbit termination, or all prime targets follows.
