# Adaptive prime guards and the lattice cost of a first Collatz rise

**PROVED scoped:** the capacity inequalities and conditional implications below.
**CITED:** Korselt's criterion and the Mavrotheris/van der Blij star-lattice
congruence. **FINITE-EXACT:** the prime census and all reported constants.
No Giuga--Carmichael integer is constructed. Universal Collatz and the actual
giant Mersenne terminal remain **OPEN**. These are filters on a proposed joint
arithmetic certificate, not proofs of convergence for sources failing it.

Artifacts: [checker](../../04-computation/experiments/joint_prime_capacity_20261014g.py),
[receipts](joint_prime_capacity_20261014g.json),
[adaptive completion](adaptive_giuga_completion_20261014g.md),
[Mersenne children](mersenne_lattice_child_guards_20261014g.md), and
[eightfold capacity](lattice_eight_capacity_20261014g.md).

## 1. Inheritance and portfolio

The closest proved mechanisms are the [signed two-prime compiler and star
lattice](giuga_block_lattice_20261014f.md), the [full old-prime CRT
interface](giuga_korselt_defects_20261014f.md), and the [native 8911 descent
ray](collatz_successor_guards_20261014f.md). The canonical hostile is 8911:
it is Carmichael and initially rises, but its Giuga index is nonintegral.
The corrected near miss is applying the odd-source lattice congruence to the
even source 30, or treating an exact eighth-step descent as a signature theorem.
The least-used sidecar is the pattern of *absent* small prime factors.

Anchor: adaptive completion with every prime obligation. Niche: source-specific
reciprocal capacity. Wildcard: the intrinsic forbidden-pair graph on the first
few admissible primes. The live board is source / prime support / index /
normality / characteristic norm / actual terminal. Relevant META-PATTERNS are
"Type every analogy and every implication" and "Audit saturation and basis
covariance before naming a lattice bridge". The graph here records a genuine
binary incompatibility, with no invented tournament orientation.

Incoming commits included audited Fermat/circulant corrections (MISTAKE-608)
and [pre-audit D19/Wilson work](d19_parity_shift_synthesis_20261014.md).
They reinforce retaining the exact
level and parity of a congruence. No pre-audit claim is a dependency here.

## 2. The complete joint source has a strong first-rise restriction

Let a composite squarefree integer N satisfy both

\[
 p\mid N/p-1,\qquad p-1\mid N-1\qquad(p\mid N).
\]

Write r for its number of prime factors and

\[
 \kappa=\sum_{p\mid N}\frac1p-\frac1N\in\mathbb Z_{>0}.
\]

Korselt makes N odd. For N>1, the odd accelerated Collatz map U initially
increases exactly when N=3 mod4. In that case every prime factor is 3 mod4:
a factor 1 mod4 would force 4|(N-1). Consequently r is odd. Reducing
`kappa N=sum(N/p)-1` modulo two shows that kappa is even, hence at least two.
In particular

\[
 \sum_{p\mid N}\frac1p>2.                                      \tag{1}
\]

The support is also **normal**: two of its primes cannot satisfy p|(q-1),
since that would make p divide both N and N-1. Normality is necessary, not
sufficient for Korselt.

The star matrix with centre kappa, prime diagonal leaves, and centre-leaf
entries one has determinant -1. For odd N its central vector is characteristic.
The cited theorem gives

\[
 r\equiv\kappa+1\pmod8.                                       \tag{2}
\]

This application is credited to [Mavrotheris, *Characteristic-Vector
Congruences for Weak Primary Pseudoperfect and Giuga numbers*](https://zenodo.org/records/22085830),
using [van der Blij's congruence](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/vdb.pdf).
The complete integral-index hypothesis is essential; arbitrary Collatz sources
do not acquire it by choosing a nearby integer for the centre.

## 3. A finite support graph gives a bound on all larger primes

The first four primes 3 mod4 are 3,7,11,19. Join two when one divides the other
minus one. The edges are precisely {3,7} and {3,19}: **a three-vertex path plus
an isolated vertex**. Its ten independent sets are the possible exact marked
support patterns. Here that familiar four-vertex object has a concrete role:
it partitions every normal support, including the absent-prime information.

For a pattern C, admit its marked primes and every later prime p=3 mod4 such
that no q in C divides p-1; exclude all other marked primes. This deliberately
relaxes incompatibilities among later primes, so every actual support in the
branch remains admissible.

Let L=10^7, h be the number of admissible primes at most L, and q the first
admissible prime above L. Compute the directed rational upper bound

\[
 H=2^{-64}\sum_{p\le L,\ p\text{ admissible}}
                 \left\lceil\frac{2^{64}}p\right\rceil.
\]

For any r-element admissible support, sorting the primes gives

\[
 \sum\frac1p\le H+\frac{\max(0,r-h)}q.                         \tag{3}
\]

For r>=h, the largest possible reciprocal sum uses all h smaller admissible
primes before any larger ones. Replacing its larger primes by q only increases
that sum. For r<h, H alone is an upper bound. Thus (3) applies to **every**
larger prime, not just the finite sieve.

For a required index b, put

\[
 R_C=h+\left\lfloor(b-H)q\right\rfloor.                        \tag{4}
\]

If H<b, every r<=R_C is excluded because the reciprocal sum is at most b,
whereas a Giuga integer of index at least b needs a sum strictly greater than b.
All sums, ceilings, floors and comparisons are exact integers in the checker.

For the ten first-rise branches, b=2. The weakest branch has C={7,19},

\[
 h=261735,\quad q=10000079,\quad
 H=\frac{15547960639220426166}{18446744073709551616},
 \quad R_C=11833260.
\]

The next six ranks still have upper sum less than four. Since the index is
positive and even it would have to be two, so (2) permits only r=3 mod8.
The first such rank above R_C is 11833267. Every other branch gives a stronger
bound (all ten receipts are saved). Therefore:

**Conditional descent theorem.** Any full Giuga--Korselt integer with fewer
than **11,833,267 distinct prime factors** satisfies U(N)<N. Equivalently, a
joint source that initially rises needs at least that many factors.

This is not a ROOT certificate for an arbitrary source and not an optimality
or literature-priority claim. The bound combines a proved tail inequality with
a finite exact census. Its cutoff can be increased without changing the proof.

## 4. Odd Mersenne exponents force the index to a multiple of four

Let N=2^E-1 with odd E>=3, and suppose it is a full joint source. The previous
section gives p=3 mod4 for every p|N. Also ord_p(2)|E is odd, hence divides
(p-1)/2. Euler's criterion and the supplementary quadratic reciprocity formula
for 2 force p=7 mod8. Therefore N=7 mod8, r is odd, and N/p=1 mod8. Reducing
the cofactor identity modulo eight gives

\[
 -\kappa\equiv r-1\pmod8.
\]

Together with (2), this gives 2kappa=0 mod8, hence

\[
 \boxed{\kappa\equiv0\pmod4,\qquad\sum1/p>4.}                 \tag{5}
\]

The [Mersenne companion](mersenne_lattice_child_guards_20261014g.md) independently
derives this same source-specific restriction. This is the genuine interaction
between the lattice and Mersenne arithmetic; it is not an inference from the
number of steps in the 8911 certificate.

Repeat the exact capacity argument with primes 7 mod8, b=4 and marked primes
7,23,31,47. The only forbidden marked pair is {23,47}, giving twelve branches.
The weakest branch is C={7,23,31}, with

\[
 h=127864,\quad q=10000223,\quad
 H=\frac{8760151626555056513}{18446744073709551616},
 \quad R_C=35379762.
\]

Until the next permitted rank, the sum is less than eight; (5) pins the index
at four and (2) pins the rank at 5 mod8. Thus any such joint odd-exponent
Mersenne number has at least **35,379,765 distinct prime factors**.
All its factors are at least seven. Since 7^5>2^14,
`2^E>N>=7^r>2^(14r/5)`. Oddness of E gives the additional necessary condition
**E>=99,063,343**. No factorization of an enormous Mersenne integer is used.

The actual uncovered minimum E=99708991741 is beyond this bound. More sharply,
the companion rejects its Korselt eligibility by an exact modular test on E
itself, and its signed-index eligibility by an order-group bound. Grouping
prime divisors by d=ord_p(2)|E gives a reciprocal contribution less than1/5
per nonempty d>1 group: their product divides 2^d-1 and each prime is at least
2d+1. Hence S(M_E)<(tau(E)-1)/5. A composite M_E with tau(E)<=6 cannot be Giuga;
a full joint odd-exponent source must have tau(E)>=22. The actual E has four
divisors, so neither signed integral-index lattice exists at that source.
In particular, a squarefree exponent in this hypothetical joint class needs
at least five prime factors (2^4<22), and a prime-power exponent ell^a needs
a>=21. These are conditional exponent restrictions, not existence statements.
These exclusions do not ground its Collatz child. A successful grounding still
needs a legal source-to-child word and an actual terminal proof.

## 5. What eight and 56 preserve

The [even-core companion](lattice_eight_capacity_20261014g.md) retains the right
characteristic vector for source30 and now improves the earlier56-new-prime
threshold to **57**. A proposed56-block must contain the36 primes7..167 and
all remaining primes are at most793. At the enlarged mandatory core M, every
old-prime condition forces the remaining product B into one class modulo M.
Its least positive representative exceeds the product of the20 largest
available primes, which proves impossibility. This proof uses retained guards
and capacity; it does not require van der Blij. The separate even-lattice
condition forces Q=31 mod240 for the original56-prime block and survives as a
necessary filter. At57,24 new primes7..103 are mandatory,33 remain, and the next
prime is between107 and185. The exact pending residue and index are saved as an
OPEN completion obligation, not a successful extension.

The earlier56 bound concerns a *proper multiple* of an even Giuga source;
it is not the odd-source rank congruence. The 8911 eight-step receipt
counts actual trajectory steps. Equations (2) and (5) instead concern a lattice
signature and the number of prime leaves. These coordinates can coexist in a
marked certificate, but identifying them loses information.

The [adaptive compiler](adaptive_giuga_completion_20261014g.md) now retains
the partial product, integer deficit, target index, remaining prime count,
normal-support constraints and the pending suffix-product guard. Partial blocks
may fail a terminal condition that a later prime repairs. The count must still
be fixed for its proved termination guarantee; allowing an unbounded count
remains a separate coverage obligation.

## 6. Verification

Run `python -B -X utf8 04-computation/experiments/joint_prime_capacity_20261014g.py`
and the same command with `-O`. Both give **21,043 checks**, identical saved JSON,
LF SHA256 `652a101df6341edc5d3918b219966989e37f76fa5145431c86efe895eabc9469`.
The prime sieve is checked against independent trial division through20000;
successor primes and every eligible intervening integer are checked separately.
Small directed sums are compared with exact Fraction arithmetic. The complete
ten- and twelve-branch universes, graph types, floors and rank exclusions are
retained in the receipt. No floating-point inequality enters a proof.

The hostile 8911 has index1874/8911, so its first rise does not contradict either
conditional theorem. No result here assumes that a Giuga--Korselt source exists.

Independent audit: the lattice lane checked the full proof and API, replayed
optimized execution, and independently used an odd-only sieve through10001000
to reproduce all22 branch counts, successor primes and exact ceiling sums.
The adaptive and even-core packets received independent proof/API/output audits
as well. The Mersenne companion records the audit of its separate source reader.
