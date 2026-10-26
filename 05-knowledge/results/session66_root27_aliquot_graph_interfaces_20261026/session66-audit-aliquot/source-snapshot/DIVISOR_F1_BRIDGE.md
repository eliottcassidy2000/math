# A weighted monoid carrier, the 14/15 split, and complement defect

2026-10-26. **PROVED elementary synthesis**, with **FINITE-EXACT** controls.
No claim of a new classical divisor theorem, an F1 field structure on inverse
aliquot fibers, or a result on the existence/value of untouchable density.

## Inheritance and target

The user's examples are fibers of `s(n)=sum_{d|n,d<n} d` for positive integers
n, with `s(1)=0`. We treat inverse fibers as sets, not a new arithmetic
operation on their output labels. The closest inherited mechanism is the
separation of a multiplicative carrier from an additive evaluation in
`normseven_f1_20261025/KUMMER_AND_F1.md` and `TWIST_AVERAGE.md`.
The canonical hostile here is `4 -> 3` versus `9 -> 4`, despite isomorphic
unlabelled divisor chains. The least-used coordinate is the prime-labelled
weight of each divisor, not merely the exponent vector or its support size.

Concept board: prime monoid / additive evaluation / inverse fibers /
complement involution / minimum orientation defect / native rank guards.
This note covers the first five. It does not supply a root-27 rank.

## 1. The precise F1-compatible construction

Positive integers under multiplication form the free commutative monoid on
the primes, by unique factorization. For `n=product p_i^a_i`, its divisor
interval is the box `0<=e_i<=a_i`. In the monoid ring define

    A_a(X_1,...,X_r) = product_i(1+X_i+...+X_i^a_i) - product_i X_i^a_i.

The subtraction simply deletes one monomial; all remaining coefficients
are one. For n=1 this polynomial is zero. Thus

    A_a(p_1,...,p_r)=s(n),
    A_a(1,...,1)=tau(n)-1.

Proof: the monomials biject with proper divisors. Evaluating at prime
weights sums their integer values; augmentation at one counts them.
This is a monoid carrier plus an additive evaluation, not an assertion
that the two evaluations coincide. It uses the elementary monoid-ring
construction in Deitmar, *Schemes over F1*, §1/Theorem1.1,
https://arxiv.org/html/math/0404185v7 . The paper defines its F1 as the
trivial monoid and base extension by the monoid ring; we do not identify
that object with the zero ring, or identify the arithmetic numbers 0 and1.

**No-go for forgetting prime weights.** Permuting the generators 2 and3
of the free monoid maps4 to9, while s(4)=3 and s(9)=4. Consequently s
cannot factor through the quotient by arbitrary prime permutations, nor
through the isomorphism class of the unweighted divisor interval.
An arithmetic evaluation retaining each prime weight restores it.

**Fiber antichain.** If m=kn, k>=2, then

    s(m)>=s(n)+n,        s(m)>=k*s(n)+1.

For the first inequality all proper divisors of n and also n are distinct
proper divisors of m. For the second, multiply each proper divisor of n
by k, then adjoin1, which is distinct from these multiples. Hence each
fiber s^{-1}(a) is a divisibility antichain. This does not imply small or
finite fibers: s^{-1}(1) is the infinite set of primes. It gives no density
for missing fibers. The equal-output relation is not a monoid congruence:
s(2)=s(3)=1 but s(2*2)=3 differs from s(3*2)=6. Thus multiplying
representatives does not define a well-defined multiplication on fibers.

## 2. Why 14 and15 fit together inside the four-prime cube

Let N=2*3*5*7=210. Its divisor box is the Boolean four-cube of size16.
Its proper-divisor set has15 elements. Its **interior** divisor set

    D°(N)={d:d|N,1<d<N}

has14 elements and is stable under `iota(d)=N/d`. Its seven pairs are

    (2,105), (3,70), (5,42), (6,35), (7,30), (10,21), (14,15).

In particular the proposed split `{3,5}|{2,7}` gives the pair `(15,14)`.
Among the three unordered two-prime/two-prime splits their product pairs
are `(6,35),(10,21),(14,15)`, with gaps29,11,1. The last is the unique
closest complementary pair among *all* seven, by direct ordered divisor
comparison. This is an exact property of210, not a universal phenomenon
for any four primes.

The complement pairs retain the additive weights too:

    s(210)=1+107+73+47+41+37+31+29=366.

**Boundary hostile:** complementation is NOT a permutation of the
15-element proper-divisor set, because it sends1 to210. One cannot
apply a permutation-defect theorem on those15 vertices without altering
the set or operation. Full divisors16 or interior divisors14 are valid.
This 15 is neither the15 Ogg primes nor the15 Richelot kernels; no
bijection with their arithmetic predicates is claimed.

## 3. Complementation and the minimum defect theorem

For N>1 let V=D°(N), v=tau(N)-2 and f=1 if N is a square, else0. The
complement involution has f fixed points and b=(v-f)/2 swapped pairs.
For an orientation T of all unordered pairs of V, define

    D_iota(T)=#{(u->v) in T: iota(u)->iota(v) is not in T}.

Then

    min_T D_iota(T)=b.

This is a specialization of Session65's involution-defect lemma, also
proved here: every swapped vertex pair forces its own edge to reverse,
contributing one. Every other unordered-pair orbit has length two (or
has both endpoints fixed); choose its orientations equivariantly. Thus
no other reversals are necessary and the bound is attained. For two fixed
endpoints there is no forced reversal. The polynomial version from that
lemma is, with M=b(f+b-1),

    sum_T z^D_iota(T) = 2^[binom(f,2)+b+M] z^b (1+z²)^M.

Here f<=1, but the inherited formula holds for arbitrary involutions.
For N210 on14 interior vertices: f0,b7,M42, polynomial
`2^49 z^7 (1+z²)^42`; the minimum is7, the maximum91. The finite checker
constructs an explicit minimizing tournament, without searching2^91 cases.

**Crucial scope:** iota is divisor complementation, not arithmetic
Frobenius. The transferable theorem concerns a specified permutation;
the supersingular geometric interpretation does not transfer just because
the same defect formula applies. Nor is the theorem a canonical choice
of arithmetic tournament. Indeed the canonical *numerically increasing*
tournament on interior divisors has defect binom(v,2), since N/d reverses
all numerical inequalities. At210 it has91, not7. The minimum in the
statement must never be omitted.

The complement-invariant trace `d+N/d` is equal on paired vertices. It
cannot orient those seven pairs without extra information. This is an
intrinsic tie obstruction, not license to manufacture a tournament.

## 4. Interfaces and limits

- Divisor box -> proper-divisor polynomial: preserves every labelled
  divisor and its weight. Augmentation loses prime weights; retain the
  evaluation homomorphism to recover aliquot sums.
- Full divisor box -> interior: preserves complementation; deleting only
  the top does not. The boundary choice is essential.
- Complement permutation -> minimum tournament defect: preserves its
  pair-orbit obstruction, not actual arithmetic orientation or dynamics.
- Fiber antichains -> candidate pruning: divisibility excludes duplicate
  output targets along a divisor chain; neither no odd exceptions nor a
  density formula follows.

The general research pattern is already in META-PATTERNS: controlled
forgetting requires a sidecar, and the object must be corrected before
sharpening a bound. No new meta-pattern is promoted from this single lane.

## Reproduction and source scope

`python3 -I -B divisor_carrier.py` and `python3 -I -B -O divisor_carrier.py`.
The stdlib-only checker verifies polynomial/divisor evaluation, complement
fixed points, a minimizing orientation and actual increasing-orientation
defect for N2..2000; 9750 strict-divisibility controls; exhaustive all
orientations for interiors of N2,4,6,12,24,30; exact210 pairs and hostile
prime relabel/boundary/tie controls. Explicit checks survive optimization.
Finite tests support the proofs, not their universal quantifiers.
Source pins are in SOURCES.json; external source is cited, not bundled.
