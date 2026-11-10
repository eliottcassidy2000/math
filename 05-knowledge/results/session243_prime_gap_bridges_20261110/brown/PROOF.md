# Brown numbers: a finite exact inverse with the prime-power sidecar

**PROVED elementary equivalence; FINITE-EXACT census n=4,...,60.**
No assertion that these are all Brown numbers, that every Brown root is prime,
or that a scalar square root can be realized by a tournament operation.

## 1. Inherited object, not the other Brocard conjecture

The repository's `collatz_doubling_tower_brocard_20260927.md` Section3 already
uses n! = (m−1)(m+1) and the divisor involution. Its root examples are
(n,m)=(4,5),(5,11),(7,71). This is Brocard's **problem**, distinct from four
primes between consecutive prime squares. The present contribution is the
explicit coprime prime-power partition inverse and its interface comparison.
No literature-firstness is claimed for this elementary reformulation.

## 2. The exact inverse

Fix n>=4. Write n! = 2^e O, with O odd and e>=3, and put C=n!/4.
Let the indivisible blocks be 2^(e−2), and r^v_r(n!) for every odd prime r<=n.

There is a bijection between positive solutions m to n!+1=m² and selections
of these blocks into two factors A,B satisfying

    AB=C,  gcd(A,B)=1,  |A−B|=1,

with the even block assigned to A. The inverse is m=A+B. Thus only the
2^(pi(n)−1) selections of the odd blocks need be checked, not all divisors
of C and certainly not all integers below sqrt(n!).

Proof: m is odd. The numbers (m−1)/2 and (m+1)/2 are adjacent positive
integers, hence coprime, and their product is C. Exactly one is even, giving
a unique naming A,B. Coprimality forces every whole prime power of C onto
one side, proving completeness of the selections. Conversely adjacent A,B
have (A+B)²−1 = (A−B)²−1+4AB = n!, proving the inverse. There is no ordering
or multiplicity ambiguity after placing the even block on A. The boundary
n>=4 is retained because e>=3 ensures a nontrivial even block.

This is a finite, terminating **membership algorithm at a specified n**.
The number of sign selections grows with n; the existence of this exact
inverse does not decide whether any solution occurs after an arbitrary bound.

## 3. What the primorial wheel forgets

A Brown root m is coprime to n!, so every prime factor of m exceeds n. This
is the same small-prime roughness predicate as for n!+1 or P_k+1, but it is
far weaker than the Brown equation. The stronger sign conditions are

    m = +1 or −1 (mod r^v_r(n!)) for every odd prime r<=n;
    m = +1 or −1 (mod 2^(e−1)).

The odd-prime statement follows because gcd(m−1,m+1)=2, so an odd prime
power cannot split across the two factors. Of the two even neighbors,
exactly one has 2-adic valuation1 and the other has valuation e−1. These
conditions are also sufficient for n! | m²−1. They establish divisibility,
not equality. For example n=4 and m=23 satisfy all sign congruences, but
23²−1=22·4!, whereas at n=5 the value m=7 survives gcd(7,5!)=1 yet fails the sign
condition at5. Both losses are explicit: the wheel drops signs and exponents, and
sign congruences drop the quotient (m²−1)/n!.

The inverse in Section2 restores exactly the lost quotient through the
adjacency condition |A−B|=1. Its arithmetic factors are not paths and do
not make a tournament merely by relabeling the blocks.

## 4. The Chan–Pak consequence and its limit

An actual tournament with H=m can be joined ordinally to a second copy,
constructing H=m². The converse is false:49 is physically attained, while7
is excluded by the inherited all-order omission theorem. Likewise ordinary
factor extraction from an attained count is not a graph constructor.

For the known Brown roots5,11,71, independent physical suppliers can be
combined to witness25,121,5041. The factorial identity then becomes a true
relation between these **already supplied** path counts; it is not how the
root suppliers are obtained. The exceptional Euclid prime2·3+1=7 demonstrates
that primality plus small-prime roughness alone cannot imply attainability.

Connection contract: Brown solution -> coprime adjacent prime-power partition
preserves the exact equation and has an explicit inverse. Brown solution ->
rough integer preserves only exclusion of small prime divisors and loses
factor exponents, signs and equality quotient. A supplied tournament -> its
double ordinal join preserves squaring, but reverse passage loses physical
realizability. No map from a general Brown solution to a tournament is proved.

## 5. Validation

`verify.py` uses only integer arithmetic and explicit runtime gates. It checks
all block selections at every n=4,...,60, independently compares the result
with direct integer-square-root testing, and checks positive/hostile controls.
The 57-value finite census finds exactly the three inherited examples. Normal
and optimized modes must give identical JSON. `--fail` deliberately fails.
