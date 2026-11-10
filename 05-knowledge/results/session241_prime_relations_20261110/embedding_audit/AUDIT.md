# Independent hostile audit: arithmetic embeddings and repunit chart

**Status: complete construction and repunit proofs PASS; finite controls PASS.**

The auditor wrote a separate verifier without importing producer code or results.
The proofs reviewed are `../constructions/PROOF.md` and `../repunit_chart/PROOF.md`. The inherited fan source is
`01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, statement4.
No new literature-priority claim is made.

## Inheritance and audit board

Closest mechanism: near-transitive path subset bijection, ordinal-product
factorization, and THM-4748's strong fans. Hostiles: injective weights15,51,85
with no private primes; normalized repunit R3=21; induced E3 inside E4 despite
3 not dividing5. Corrected near miss: numeric multiplicativity does not mean
canonical ordered tournament representatives multiply without sorting.
Least-used sidecar: the prime-index dictionary needed to interpret a Fermat
factor as an ordinary source prime.

Board: valuation supports; strong-component words; gcd normalization;
induced containment; family-relative vertex cost. The construction preserves
specified scalar predicates and loses addition, numerical primality, and the
raw induced order. It does not solve Chan–Pak or force a divisibility order
onto tournament vertices.

## 1. Construction proof review

For a>=1, reversal of the endpoint arc in the transitive tournament of order
a+2 leaves one increasing Hamilton path and 2^a paths using the backward arc.
Each latter path is specified by the subset of interior vertices placed
before the backward arc. The a=0 boundary is excluded: a two-vertex tournament
has one Hamilton path, not2.

For positive a,b, their common divisor g leaves coprime quotients A,B.
If A,B have opposite parity, any common divisor of 2^a+1 and 2^b+1 divides2
and is odd. If A,B are odd, Bezout gives a linear combination with odd sum of
coefficients; invertibility modulo any common divisor gives 2^g=-1 there.
This proves the exact gcd formula and hence the distinct-v2 criterion.
For k>=1 distinct valuations sorted t_j obey t_j>=j, giving exponent sum
at least2^k-1. Equality forces exponents1,2,...,2^(k-1). The cost statement
is correct within this specific bank, not among all tournaments.

Fermat numbers are pairwise coprime whether or not individually prime.
Assigning one to each ordinary prime gives a multiplicative injection and
preserves gcd, lcm, coprimality, and divisibility. Actual counts arise by
ordinal sums. The producer correctly distinguishes count multiplicativity
from ordered component-word multiplication, which needs canonical sorting.

## 2. Infinite families and private-prime quantifiers

Let I be any index set and w_i positive integers greater than1. The domain is
N^(I), the **finite-support** vectors. Each coordinate has a private prime
if some prime dividing w_i divides none of the weights w_j with j!=i.
This is a condition relative to the entire family, not just a current prefix.

Private primes suffice to read every coordinate inequality from valuations.
Conversely, if coordinate i has no private prime, the finite set of prime
factors of w_i can each be covered by a different-from-i weight. Taking
sufficient powers of this finite selection produces a finite-support vector
with zero i-coordinate whose image is divisible by w_i. Thus reflection
fails. No compactness assumption, uniform exponent bound, or finitely many
weights assumption is needed. The producer's countable scope is sound and
can be strengthened to arbitrary I (a successful family is necessarily
countable because distinct coordinates require distinct private primes).

Gcd or lcm preservation forces pairwise coprimality by testing two distinct
unit vectors. Conversely disjoint prime supports make min/max commute with
the encoding. The >1 hypothesis matters to injectivity and coprimality
reflection. Our verifier deliberately includes unit weights as hostile
controls: their empty prime support yields a divisibility-reflection failure.

Injectivity alone merely requires finite rational independence of valuation
columns. Weights15,51,85 are independent (determinant -2), yet15 divides51*85.
Thus multiplicative independence cannot replace private primes.

## 3. Repunit chart: independent all-order derivation

For n>=2, THM-4748's F_(2n) has chain vertices1,...,2n+1 and apex0
losing exactly to n,n+1,2n+1. It is strong, has order2n+2, and
H_n=4^n-1. Set the n=1 carrier to C3, with H_1=3.
The explicit cycle0->1->...->2n+1->0 verifies strength for n>=2.

For positive a,b and integer q>1, Euclidean subtraction on exponents gives

    gcd(q^a-1,q^b-1)=q^gcd(a,b)-1.

Indeed when b>=a, subtract q^(b-a)(q^a-1) from q^b-1 to obtain
q^(b-a)-1, and iterate the Euclidean algorithm. The base case is gcd(x,0)=x.
It follows that H_a divides H_b exactly when a divides b: one direction is
the geometric sum; the reverse uses the gcd identity and strict monotonicity.

Set R_n=H_n/3, a positive integer. Since gcd(3x,3y)=3gcd(x,y),

    gcd(R_a,R_b)=R_gcd(a,b).

Again strict monotonicity proves R_a divides R_b iff a divides b, including
a=1 where R_1=1. This normalization is arithmetic only: R_3=21, an inherited
forbidden tournament H-value. It must not be represented as a new H-spectrum
construction. H_1=3 is not a multiplicative unit; R_1=1 is.

Finally gcd(H_a,H_n)=3 iff gcd(a,n)=1. Counting a in1,...,n therefore gives
Euler phi(n), including phi(1)=1. The range endpoint matters at n=1.
The chart preserves gcd/divisibility but is not multiplicative: R_2=5,
R_4=85 !=25. Ordinary lcm is also not preserved: lcm(R_2,R_3)=105 !=R_6=1365.
One may transport operations inside the chart, but these are not ordinary
integer lcm/multiplication or native tournament operations.

## 4. Exact universe and reproduction

Run from repository root:

    python3 05-knowledge/results/session241_prime_relations_20261110/embedding_audit/verify.py

This writes `verification.json`. Its independent recursion counts paths by
remaining vertices and last vertex, without producer imports. Universes:

- Actual near-transitive adjacency: a=1,...,9 (order<=11).
- Actual fan adjacency and spanning cycles: n=1,...,6 (order<=14).
- Near gcd and repunit gcd/divisibility: every pair1<=a,b<=100.
- Totient interpretation: n=1,...,100.
- Bank optima: every k-subset of exponents1,...,16, k=1,...,5.
- Private-prime obstruction: all454 multisets of weights1,...,12 of size1,...,3;
 675 coordinates with no private prime each get an explicit finite witness.

Finite tests are controls, not proofs of all-order or infinite-family claims.
No theorem IDs are assigned by this audit. Git integration belongs to parent.

## 5. Complete repunit text follow-up review

The complete `../repunit_chart/PROOF.md` is accepted. The all-index image-totient
classification is sound: if an odd prime r divides n, every prime q dividing
2^r-1 has order exactly r for base2, since the only other divisor of prime r
is1 and no prime divides2-1. Thus r divides q-1, q divides4^n-1, and phi(4^n-1)
has an odd prime divisor. Otherwise n=2^j, and 4^n-1=F0*...*Fj. For j<=4
these are the five known Fermat primes, so their totients multiply to a power
of2. For j>=5 the divisor641 contributes640 to the totient, excluding a pure
power of2. No assumption about future Fermat primes is used.

For an odd prime p and d=ord_p(4), visibility is exactly d|n. When n=dt,
odd-prime LTE applies to (4^d)^t-1 because p divides4^d-1 and p does not
divide4. This yields the displayed valuation formula without exceptional
odd primes; p=3 gives v3(h_n)=1+v3(n).

The divisor-count and incidence-Mobius statements follow from the proven
poset isomorphism on each finite divisor interval. They are transport
statements rather than numerical equalities at h_n. The shifted addition
law is an exact expansion of (h_a+1)(h_b+1).

Independent exhaustive injection search finds the induced Q2-to-Q3 map
[0,1,3,4,5,7], and verifies all36 ordered adjacency entries. The source index2
does not divide3, and15 does not divide63, establishing the claimed hostile.

The sparse chart count is floor(log_4(X+1)) for **real X>=0** (in particular
for positive integer thresholds). This domain is necessary; negative X was
not intended and does not satisfy the unqualified expression. The parent
was asked to state the domain explicitly. Monotonicity gives the formula
by solving4^n-1<=X with n>=1; there are no rounding or duplicate-value issues.

The independent verifier now additionally checks phi(h_n) for n=1,...,16;
the641 obstruction at n=32; prime visibility and LTE for all25 odd primes
through101 and n=1,...,200 (5,000 pairs); the induced witness via its own
exhaustive search; and the sparse-count formula at every integer threshold
0,...,4096 using exact integer powers instead of floating-point logarithms.
These checks do not import the producer's algorithm or data.
