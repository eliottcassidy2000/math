# Fixed excess two and square-sum aliquot edges: a 356-prime obstruction

Date: 2026-10-28. Status: **PROVED elementary necessary condition**, with a **FINITE-EXACT rational certificate**. This is a derived connection, not a claim of literature priority. No existence or nonexistence theorem is asserted.

## Result

Let `s(n)=sigma(n)−n` for positive integers. If

` s(n)=n+2 ` and ` n+s(n) is a square, `

then **n is odd, has at least 356 distinct prime factors, and has at least 1156 decimal digits**.

In particular, **no even fixed-excess-two source is a square-sum aliquot edge**. The familiar `20→22` has total42, not a square; the entire even family `2^a(2^(a+1)−3)` and the outside-family source650 obey the same obstruction.

The 356 is the exact crossing of the specified prime-product bound. It is **not** an existence result, an upper bound, or a proof that 356 is attainable. The possible odd intersection remains unresolved here.

## Proof

The unit n=1 fails the premise, since s(1)=0 is not 1+2. Thus any solution has n>1.

First n=1 fails the equation because sigma(1)=1. Suppose now `sigma(n)=2n+2=z²`. The even square z² has even root, say `z=2k`, so

`n=2k²−1`, `sigma(n)=4k²`.

Thus n is odd. If a prime p divides n, then p is odd and

`(2k)² ≡ 2 (mod p)`.

Hence 2 is a quadratic residue modulo p. The supplementary law

`(2/p)=(-1)^((p²−1)/8)`

therefore forces **every prime factor of n to be 1 or 7 modulo8**. This standard law also follows from Gauss's lemma: the relevant count is `(p−1)/2−floor(p/4)`, whose parity is even exactly in these two residue classes.

Let `q1<q2<...` be all primes congruent to1 or7 modulo8, and write `r=omega(n)`. Ordering n's distinct prime factors and using the fact that `x/(x−1)` decreases for x>1 gives

`2 < sigma(n)/n = product_(p^e||n) (1+p^(-1)+...+p^(-e))`

`  < product_(p|n) p/(p−1) <= product_(j=1..r) qj/(qj−1).`

The strict upper inequality holds because all prime exponents are finite. Exact integer arithmetic certifies

| r | q_r | rigorous enclosure for the product |
|---:|---:|---|
|355|5431|`1.999836 < product < 1.999837 < 2`|
|356|5441|`2 < 2.000203 < product < 2.000204`|

Thus r<=355 is impossible. Finally n contains at least356 distinct allowed primes, so

`n >= product_(j=1..356) qj >= 10^1155`.

The same exact computation shows this prime product has1156 decimal digits. Repeated prime powers can only increase n.

## Guards and hostile controls

* `n=2k²−1` is necessary, not sufficient: k=2 gives n=7, but sigma(7)=8 rather than16. The divisor-sum equation cannot be discarded after the norm substitution.
* The all-primes residue restriction uses both predicates. It does not apply to a generic fixed-excess-two source:650 has the prime factor5. It also does not apply to a generic square-sigma source:115 has sigma144 but has excess−86.
* The 2013 Li–Ma–Zhang theorem that an odd solution of sigma(n)=2n+2 has at least6 distinct primes concerns the larger class without the square-sigma restriction. This result is not a replacement for that theorem: its much larger bound uses the extra predicate.
* The elementary congruence excludes3 and5 automatically. It does not exclude all sufficiently large prime supports; the prime-product envelope already exceeds2 at356 terms.
* Distinct source and successor are automatic under excess2, and n>1 with positive successor is automatic for an actual solution. Thus there is no hidden loop/zero-vertex issue in the positive simple square-sum graph.

## Reproduction

Standard-library Python3, no input artifacts:

```
python3 check_square_excess.py
python3 -O check_square_excess.py --output /tmp/square_excess_optimized.json
cmp square_excess_results.json /tmp/square_excess_optimized.json
```

The checker generates every prime through5441 both by a sieve and by independent trial division; compares the full lists; restricts to the required residue classes; identifies the first exact rational crossing; records the reduced numerators, denominators and signed gaps from2; verifies the decimal enclosures with rational inequalities; and records the exact lower-bound integer. It also directly checks the quadratic-residue criterion for every odd prime through5441. Explicit exceptions survive optimized Python. Generated JSON is a reproducibility receipt, not an input dependency. A separately implemented independent audit in `../audit-square-excess/` reproduces the full allowed-prime list, both rational boundary inequalities, and the lower-bound integer; the audited result is PASS.

Connection contract: source predicates are fixed aliquot excess and square successor total; map is `n=2k²−1`; preserved necessary data are the residue classes of all prime factors; lost data are factor exponents and the exact divisor-sum equation; repair is to retain both original equations. Cheapest hostile is n=7; cheapest decisive numerical gate is the exact355-prime product, not a floating-point approximation.
