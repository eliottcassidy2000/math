# Exactly seven proper-divisor terms: a complete abundance classification

Date: 2026-10-28. Status: **PROVED** elementary classification and infinitude;
**FINITE-EXACT** independent verification. No literature-priority claim.

## Result

Among positive integers with exactly **seven distinct positive proper
divisors**, there are **no perfect numbers**, and the abundant ones are exactly

**`{24,40,54,56,70,88,104} ∪ {6q : q≥5 is prime}`.**

The two sets are disjoint. Every other seven-term source is deficient. Thus the
five-term abundance list `{12,18,20}` from
[DIVISOR_TERM_CLASSES.md](DIVISOR_TERM_CLASSES.md) is finite, but the seven-term
abundant class is infinite, with a simple exact reason rather than numerical
extrapolation.

**Correction lineage:** an exploratory message mistakenly multiplied `2*5*7`
as210. The correct exceptional source is **70**, whose proper divisors sum74.
The number210 has four distinct prime factors, **15** proper divisors, and
aliquot sum366; it is not in the seven-term layer at all. The product, layer
and aliquot value are all explicit hostile controls in the checker. The
erroneous numeral was corrected before this note was prepared or published.

## 1. Inherited mechanism and exact source shapes

Let `s(n)=σ(n)−n` and `T(n)=τ(n)−1`, where τ counts positive divisors including
n. The exponent-vector classification proved in
[DIVISOR_TERM_CLASSES.md](DIVISOR_TERM_CLASSES.md) gives

`T(n)=7 iff n=p⁷, p³q, or pqr`,

with distinct primes represented by distinct letters. This also follows
directly from the multiplicative partitions `8`, `4*2`, and `2*2*2` of τ=8.
Every case below is global in the primes, not bounded by the numerical check.

Live concepts: divisor-count layer / abundance / prime adjunction /
coprimality / fixed excess. The map from a perfect number to a prime multiple
explains the infinite part; the distinct-prime guard is load-bearing.

## 2. Proof for the three shapes

### Prime powers p⁷

Every prime power is deficient, since

`s(p^a)/p^a = 1/p+...+1/p^a < 1/(p−1) ≤ 1`.

### Two-prime shape p³q

Set `A=p³+p²+p+1` and `B=p³−p²−p−1`. Then

`s(p³q)−p³q = A−Bq`.

* For p=2, the excess is `15−q`. The allowed odd primes below15 are
  3,5,7,11,13, giving exactly24,40,56,88,104. Equality would require the
  nonprime q=15; larger primes give deficiency.
* For p=3, the excess is `40−14q`. Only q=2 is positive, giving54 with
  excess12. Every other allowed prime is at least5 and gives deficiency;
  equality has no integer q.
* For p≥5, use q≥2 and B>0 to bound the excess by
  `−p³+3p²+3p+3≤−2p²+3p+3<0`.

Thus this shape supplies exactly the six listed abundant values other than70,
and no perfect values.

### Squarefree three-prime shape pqr

If all three primes are odd, then

`σ(n)/n ≤ (1+1/3)(1+1/5)(1+1/7)=64/35<2`,

so n is deficient. Otherwise write `n=2pq` with odd primes p<q. Its excess is

`s(2pq)−2pq = 3+3p+3q−pq`.

* For p=3 this is **12**, for every prime q≥5. These are precisely the
  infinite family6q in the result.
* For p=5 this is `18−2q`. Since q>5 is prime, it is positive only at q=7,
  giving **70** with excess4. Equality would require the nonprime q=9.
* For p≥7, the coefficient of q is negative. Using q≥p gives the upper bound
  `3+6p−p²<0`, already −4 at p=7 and decreasing thereafter.

There are no perfect numbers in this shape. Combining the three cases proves
the classification. Infinitely many primes q≥5 exist, so the displayed
seven-term abundant family is genuinely infinite. ∎

## 3. Why the constant-excess family occurs

For any positive integer m and prime q not dividing m, multiplicativity gives

`s(mq)=q*s(m)+σ(m)=m+(q+1)*s(m)`,

and the divisor count gives

**`T(mq)=2T(m)+1`.**

The new proper divisors split disjointly into all divisors of m, and q times
the proper divisors of m. Coprimality makes this a genuine disjoint split.

If m is perfect, `s(m)=m`, so **`s(mq)−mq=2m`**. Taking m=6, with its three
proper divisors, yields seven-term sources6q with constant excess12. This is
an exact transport of a verified predicate, not a consequence of matching
term labels. The classical multiplicative-sigma identity is inherited, and
no novelty is claimed for this elementary consequence.

Guard failures are concrete: q=2 or3 divides6; the sources12 and18 then have
five terms, not seven, and excesses4 and3, not12. Also q must be prime for the
term-doubling rule: using q=9 gives54, which happens to have seven terms and
excess12 for a different factorization reason. The general theorem is not
extended merely because this particular composite example coincides.

The full classification is stronger than the family formula. For example54
also has excess12 but is not6 times a prime; the finite exceptional set must
not be silently discarded.

## 4. Exact square-sum-edge subfamily

An abundant aliquot step `n→s(n)` has positive distinct endpoints. Their sum
is square exactly when `σ(n)` is square. Within the seven-term abundant class,
the sources satisfying this additional predicate are exactly

**`{70} ∪ {6q : q=12k²−1 is prime, k≥1}`.**

Indeed, for the finite exceptions24,40,54,56,70,88,104 the sigma values are
60,90,120,120,144,180,210, respectively. Only70 gives a square. For the family
6q, multiplicativity gives `σ(6q)=12(q+1)`. If this is z², then6 divides z;
writing z=6u gives `q+1=3u²`. Since q is odd, u must be even, say u=2k.
Thus q=12k²−1. Conversely such a prime gives

`n=72k²−6`, `s(n)=72k²+6`, `σ(n)=144k²=(12k)²`.

This classification contains a **primality condition**; it does not prove
infinitely many such k. The underlying seven-term abundant family is infinite,
but its square-edge subfamily's infinitude is not established here.

The same prime condition appears in the distinct-semiprime source2q:
`σ(2q)=3(q+1)`, whereas `σ(6q)=4σ(2q)`. Since q≥5, multiplication by3 is
coprime to2q and multiplies sigma by the squareσ(3)=4. Thus it preserves and
reflects the square-sigma predicate in this guarded pair of families. Their
targets satisfy

`s(2q)=q+3=12k²+2`, `s(6q)=6s(2q)−6=72k²+6`.

This is a concrete instance of the inherited multiplicative-sigma mechanism,
not a new general closure theorem. The source shape changes from a three-term
distinct semiprime to a seven-term squarefree triple.

The smallest source in the classified square-edge subfamily is **66→78**
(k=1,q=11, total144). But the smallest *target* in this subfamily is **74**,
from70→74 (also total144). Source order and target order differ. This proves
the two minima only within **abundant seven-term sources**; a minimum over all
higher mixed shapes requires the separate complete-fiber census.

## 5. An exact counting corollary

Let `A_7(X)` count abundant seven-term sources at most real X≥0, and let π(y)
be the number of primes at most y. The disjoint classification gives

```
A_7(X) = max(π(X/6)−2,0)
         + #{n∈{24,40,54,56,70,88,104}: n≤X}.
```

For **X≥104**, this simplifies to **`A_7(X)=π(X/6)+5`**. No asymptotic theorem
or conjecture about primes is needed. This is a source count, not a count of
distinct aliquot targets: different sources can share targets.

## 6. Reproduction

Python3 standard library, standalone:

```
python3 check_seven_term_abundance.py > seven_term_results.json
python3 -O check_seven_term_abundance.py > seven_term_results.optimized.json
cmp seven_term_results.json seven_term_results.optimized.json
```

The default universe is all sources1..100000. One sieve independently sums and
counts proper divisors; a prime sieve generates the theorem's predicted
abundant set. Trial factorization classifies every actual seven-term source;
paired-divisor enumeration checks every source through2000 and all named
positive/hostile controls. The exact counting formula is verified at every
integer cutoff0..100000. Prime adjunction is checked on m=6,28,496 with
all eligible primes q≤100, and the repeated-prime boundary is checked on12,18.
The square-sigma corollary is checked over the entire abundant seven-term
source universe and against a separately generated quadratic-prime family;
source66/target78 and source70/target74 test the two distinct minima.
Checks use explicit exceptions and stay active under `python -O`.

This finite test is not the proof of infinitude or classification. Those are
the all-prime case splits above. Generated deterministic JSON is a local
verification receipt, not a source dependency.
