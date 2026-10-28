# Twenty, fixed aliquot excess, and inverse-fiber bridges

Date: 2026-10-28. Status: **PROVED elementary identities and restricted classifications**; **FINITE-EXACT census and complete small fibers**; **CITED classical classification**. No novelty claim for the two-prime excess-two classification.

## 1. What twenty does—and does not—classify

Write `s(n)=sigma(n)−n` and `Delta(n)=s(n)−n`. Twenty has

`proper divisors {1,2,4,5,10}, s(20)=22, Delta(20)=2`.

It has **five**, not three, proper-divisor terms. The user's clarified “exactly three proper divisor terms” is a different classification; fixed excess is a complementary coordinate, not a substitute for term count.

| n | s(n) | Delta(n) | number of proper divisors |
|---:|---:|---:|---:|
|12|16|4|5|
|18|21|3|5|
|20|22|2|5|
|24|36|12|7|
|28|28|0|5|

Thus20 and28 have the same divisor-count shape but different excess. The two coordinates must not be conflated.

Terminology matters. In [Pollack–Shevelev, *On perfect and near-perfect numbers*, 2012](https://www.pollack-math.net/pnp-4.pdf), near-perfect means `Delta(n)=d`, where d is one proper divisor of n. Thus12,18,20,24 have redundant divisors4,3,2,12 respectively;28 is perfect. Quasiperfect means excess+1; almost-perfect means excess−1, not+2. Powers of2 satisfy the latter equation. The label “2-near perfect” concerns omitting two divisors, **not** excess2; see the definitions in [Aryan et al., *On 2-Near Perfect Numbers*, 2024](https://math.colgate.edu/~integers/y61/y61.pdf).

The 2013 Li–Ma–Zhang paper cited below uses “near-perfect” for excess2 specifically. To avoid that nomenclature collision this note says **fixed excess two**. An even excess-two solution is near-perfect in Pollack–Shevelev's sense because2 is a proper divisor. A hypothetical odd excess-two solution would not be: its required redundant divisor2 does not divide n.

## 2. Exact two-prime factorization

For integers a,b>=1 and an odd prime p, put `R_b(p)=1+p+...+p^(b−1)`. Then

**`Delta(2^a p^b)+1 = (2^(a+1)−p) R_b(p)`.**

Proof: set A=2^(a+1). Since `sigma(2^a p^b)=(A−1)(R_b+p^b)` and `(p−1)R_b=p^b−1`,

`Delta+1=(A−1)R_b−p^b+1=(A−p)R_b`.

For b=1 this is particularly simple:

**`Delta(2^a p)=2^(a+1)−1−p`.**

Coprimality and primality are load-bearing. Putting p=9,a=1 into the prime formula gives−6, whereas n=18 has excess3.

### Complete classification for excess2 with at most two distinct prime factors

A positive integer n with at most two distinct prime factors satisfies `s(n)=n+2` **iff**

**`n=2^a(2^(a+1)−3)`, with a>=2 and `2^(a+1)−3` prime.**

Proof of completeness:

* n=1 has excess−1.
* For a prime power p^e, `sigma(p^e)/p^e < p/(p−1)<=2`, so excess is negative.
* For distinct odd primes p,q and positive exponents u,v,
  `sigma(p^u q^v)/(p^u q^v) < (3/2)(5/4)=15/8<2`.
* Thus an excess-two solution in this domain is `2^a p^b`. The identity gives `(2^(a+1)−p)R_b(p)=3`. Both factors are positive integers. If b>=2, R_b>=p+1>=4, impossible. Hence b=1 and p=2^(a+1)−3. Conversely direct substitution proves every guarded candidate.

This classification is already **Proposition2.3** of [Li–Ma–Zhang, *Odd solutions of sigma(n)−2n=2 have at least six distinct prime factors*, 2013](https://publi.math.unideb.hu/paper/1759/download/). The factorized proof above is independent and elementary. That paper's Remark2.4 also records the counterexample650 below. Its Theorem1.1 gives the stated six-prime necessary condition in the odd case; no deep proof from that theorem is used here.

Examples are20,104,464,1952,130304,522752. This is a **prime-conditioned family**, not a proved infinite family. It is not the global classification:

`650=2*5²*13`, `sigma(650)=3*31*14=1302=2*650+2`.

### Fixed general c>=0: unbounded prime sequence plus finite repeated-prime exceptions

Every nonnegative-excess n with at most two distinct prime factors has the shape `2^a p^b`, and it satisfies excess c iff

`h R_b(p)=c+1`, where `h=2^(a+1)−p`.

For b=1, `p=2^(a+1)−1−c` must be an odd prime. This can occur only for even c; a is unbounded and primality remains a separate guard.

For b>=2, `p+1<=R_b<=c+1`, hence p<=c. Also `3^(b−1)<=c+1`. Thus there are finitely many p,b. For each one, require `R_b|(c+1)`, put `h=(c+1)/R_b`, and check that `p+h` is a power of2 at least4. This is a complete finite algorithm for the **b>=2 subdomain**, with no source cutoff. It does not make the b=1 branch finite.

For c=0,1,2,4,8 there are no b>=2 exceptions. For c=3 the sole exception—and the sole two-prime solution—is18. For c=12 the sole b>=2 exception is54; the b=1 branch has `p=2^(a+1)−13` prime and starts with24. The checker derives these exception lists rather than merely searching n up to a chosen bound.

## 3. Prime dilation and its perfect-prefix exception

If q is prime and does not divide m, multiplicativity gives

**`Delta(mq)=Delta(m)q+sigma(m)`.**

For fixed m,c, a nonzero Delta(m) determines at most one q:

`q=(c−sigma(m))/Delta(m)`, guarded by integrality, primality and q not dividing m.

If m=M is perfect, the coefficient vanishes. Therefore

**`Delta(Mq)=2M` for every prime q not dividing M.**

This is an unconditional infinite family at a fixed excess: `Delta(6q)=12` and `Delta(28q)=56`. Euclid's infinitude of primes supplies infinitely many allowed q. These do **not** prove infinitely many excess-two solutions. The coprimality guard matters: q=3 with M=6 gives18, whose excess is3 rather than12; q=2 gives12, whose excess is4.

## 4. Two prime tails: fixed excess versus inverse fibers

Fix m>=1 and distinct primes p,q with `gcd(m,pq)=1`. Let S=sigma(m), A=S−m=s(m), and D=2m−S. Direct expansion gives two exact reductions:

* If A!=0, `s(mpq)=t` iff **`(Ap+S)(Aq+S)=At+mS`**.
* If D!=0, `Delta(mpq)=c` iff **`(Dp−S)(Dq−S)=2mS−Dc`**.

When t and c are the actual sum and excess, both displayed identities hold even if the coefficient is zero; the reverse implication requires the stated nonzero coefficient. These are factor-pair identities with separate congruence, positivity, distinct-primality and coprimality guards. They are not universal parametrizations without those guards.

For inverse fibers, m=1 has A=0 and must be treated separately: `p+q=t−1`. For m>1, A>0 and the right side is positive, so positive factor pairs give a finite reconstruction for fixed m,t. This is compatible with Session82's factor-prefix inverse method.

For fixed excess, if D=0 the right approach is the unfactored relation `c=S(p+q+1)`. If D!=0 and the right side is nonzero, integer factor pairs give a finite reconstruction for fixed m,c; negative factors must not be silently discarded. If the right side is zero, one factor can vanish: `Dp=S`. Under the guards this says **mp is perfect**, leaving q free and recovering the preceding infinite dilation family. A blanket assertion of finitely many sources from “factor pairs” would miss this zero-product branch.

In the special m=1 case the excess equation is `(p−1)(q−1)=2−c`. For distinct primes, no positive excess is possible; at c=0 the only product is6. This explains why positive fixed excess is not a feature of the distinct-semiprime three-term family.

## 5. Complete inverse fibers and the coordinate lost by observing s alone

Session82's bound applies: for t>=2 every source is composite; its least prime factor gives a proper divisor at least sqrt(n), so `s(n)>=1+sqrt(n)` and **n<=(t−1)²**. A sieve through this bound and independent divisor-pair enumeration produce:

| t | complete fiber `{n:s(n)=t}` |
|---:|---|
|16|12,26|
|20|34,361|
|21|18,51,91|
|22|20,38|
|28|28|
|36|24|

For instance the observation22 does not retain “excess2”: it comes from20 or38, with excesses2 and−16. At a fixed target t and prescribed excess c, the sole possible source is t−c, but it must still satisfy the divisor-sum equation. This retains the source coordinate rather than treating the fiber label as a structural type.

The adjacent [square-excess obstruction](SQUARE_EXCESS_OBSTRUCTION.md) proves that no even excess-two source makes a square-sum aliquot edge, and that a hypothetical odd intersection must have at least356 distinct prime factors. That additional predicate must not be inferred from excess alone.

## 6. Verification and boundaries

Standard-library Python3:

```
python3 check_fixed_excess.py
python3 -O check_fixed_excess.py --output /tmp/fixed_excess_optimized.json
cmp fixed_excess_results.json /tmp/fixed_excess_optimized.json
```

The checker independently checks the first10,000 proper-divisor sums by divisor pairs, tests the guarded identities, derives the complete repeated-prime exceptions, and computes all small fibers through their proved source cutoffs. A separate finite census through1,000,000 gives excess-two sources exactly

`20,104,464,650,1952,130304,522752`.

That census is not a global theorem. The global theorem is restricted to at most two distinct prime factors; the odd existence question and unrestricted excess-two classification are not settled. Output is deterministic and explicit checks survive Python optimization. The repository and Session82 source files were not modified by this lane.

Inheritance: factored sigma and Session82 prime-prefix inversion. Live board: proper-divisor count / fixed excess / prime-support count / inverse fiber / square-edge predicate. Canonical hostile:650 breaks the unguarded two-prime formula;18 breaks noncoprime perfect-prefix dilation;20 and38 break descent of excess through s. Least-used sidecar: whether a prefix is perfect, causing the linear coefficient to vanish and creating a genuine infinite branch.
