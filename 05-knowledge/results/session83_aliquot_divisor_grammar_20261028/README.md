# Session83: aliquot fibers, divisor-term classes, and guarded composition

Date:2026-10-28. **Research completed as proofs, independently checked finite computations, and explicitly scoped literature consequences.** This packet follows the clarification that the second class means **exactly three proper-divisor terms**, rather than arbitrary three-summand rewritings. All eleven checkers use the Python standard library.

## Start here

For a positive source n, write s(n)=σ(n)-n and T(n)=τ(n)-1. A target row is F_t={n:s(n)=t}. Source1 is separate; its row is0. The row1 contains all primes and is infinite. For n>=2 the disjoint partition is:

1. **I:** prime powers p^a, including primes.
2. **II:** products pq of distinct primes.
3. **III:** every other composite.

Exactly three proper divisors means p³ or pq. The prime-power priority puts p³ in I. All III sources have at least five proper divisors; the count need not be odd. For example36 has eight. The equality s(20)=22=20+1+1 does not move20 into II: its five proper divisors are1,2,4,5,10. Fixed excess s(n)-n=2 is a separate useful property.

**Main synthesis:** retain the prime-labelled divisor box, source, divisor sum, divisor count and coprimality guard while composing. Evaluate into aliquot rows only after deciding which later predicates must remain observable. This is an organizing principle, not a claim of a new universal classification theorem or an integer–tournament isomorphism.

## 1. Complete structural classifications

[DIVISOR_TERM_CLASSES.md](DIVISOR_TERM_CLASSES.md) proves that the exponent shapes of sources with k proper divisors are exactly the multiplicative partitions of k+1, with1 subtracted from each factor. In particular three-term shapes are p³,pq; five-term shapes are p⁵,p²q; seven-term shapes are p⁷,p³q,pqr.

The [low-term atlas](LOW_TERM_ATLAS.md) gives a complete all-size abundance classification for T(n)<=9:

- Five terms: abundant sources exactly12,18,20; perfect source28.
- Seven terms: abundant sources exactly24,40,54,56,70,88,104, together with6q for every prime q>=5; no perfect source.
- Eight terms: abundant sources exactly36,100,196; no perfect source.
- Nine terms: ten explicit abundant sources, and perfect source496.
- Across all layers through nine terms, the only perfect sources are6,28,496.

The [seven-term proof](SEVEN_TERM_ABUNDANCE.md) explains the infinite branch: adjoining a new prime to a perfect source M gives excess2M, and transforms its term count T into2T+1. Hence6, which has three proper divisors, generates the seven-term family6q with constant excess12. This is genuine infinitude, not an extrapolation from a search.

## 2. The entire fiber census through target2000

[FIBER_SHAPE_CENSUS.md](FIBER_SHAPE_CENSUS.md) and [CENSUS_OUTPUT.md](CENSUS_OUTPUT.md) contain the full finite results. Every fiber was compared with an independent proper-divisor sieve through3,996,001, the proved complete source cutoff. All1,999 lists agree exactly.

- 28,770 source-target incidences: I332, II26,647, III1,791.
- 196 empty fibers;405 singleton fibers;34 rows containing all three classes.
- First II/III overlap: F16={12,26}.
- First nonempty only-III row: F28={28}.
- First all-three row: **F40={44,74,81}**, with classes III,II,I respectively.
- First repeated-I row: F31 contains32 and125.

These minima are certified globally by checking every smaller target with complete fibers. The aggregate counts remain restricted to the stated window.

[Source density and fiber weights](SOURCE_DENSITY_AND_FIBER_WEIGHTS.md) proves that III has natural density one among source integers, although II dominates this bounded-target census. It supplies exact source Dirichlet-series identities and an elementary CRT/Chebyshev proof. **Source density and the density of empty target rows are different questions.** No asymptotic untouchable-density claim follows.

## 3. Twenty and a strong square-edge obstruction

[FIXED_EXCESS_AND_FIBERS.md](FIXED_EXCESS_AND_FIBERS.md) gives an independent proof of the classical complete classification of excess-two sources with at most two distinct prime factors:

\[
n=2^a(2^{a+1}-3),\qquad a\ge2,\quad 2^{a+1}-3\text{ prime}.
\]

The theorem is credited to Li–Ma–Zhang(2013), not presented as new. The source650=2*5²*13 shows why the two-prime restriction cannot be removed. Fixed-excess factor-pair formulas retain their zero-coefficient and zero-product exceptions; a perfect prefix produces genuine infinite families rather than a finite search.

The separately audited [square-excess obstruction](SQUARE_EXCESS_OBSTRUCTION.md) proves:

> If s(n)=n+2 and n+s(n) is square, then n is odd, has at least356 distinct prime factors, and has at least1156 decimal digits.

The mechanism is n=2k²-1, forcing every prime divisor to be1 or7 modulo8. The abundancy envelope from the first355 allowed primes is strictly below2; the first356 cross2. All comparisons use exact integers/rationals. This is a necessary condition, not an existence theorem or an attainable bound. In particular20→22 is not a square-sum edge.

## 4. The concrete2–3–11 and2–5–7 connection

[CARRIER_AND_144_BRIDGE.md](CARRIER_AND_144_BRIDGE.md) ties the factor grammar directly to the earlier triple themes:

\[
\sigma(22)=36,\quad\sigma(66)=\sigma(70)=144,
\]

where22=2*11,66=2*3*11,70=2*5*7. Multiplying22 by the coprime factor3 multiplies sigma byσ(3)=4, lifting the three-term edge22→14 to the seven-term edge66→78. The different triple70 yields70→74. Prime-exponent parity vectors of sigma explain the cancellations; a single “square/non-square” bit does not compose.

The complete square-total fiber is

\[
\sigma^{-1}(144)=\{66,70,94,115,119\}.
\]

The abundant seven-term square-sigma sources are exactly

\[
\{70\}\cup\{6(12k^2-1):k\ge1,\ 12k^2-1\text{ prime}\}.
\]

Infinitude of this quadratic-prime subfamily remains unproved here. The least III source making a square-sum aliquot edge is66, while the least corresponding target is74, from70. Those minima differ because s is not monotone.

## 5. Repunits,11,25, and the carefully typed7/21 connection

[PRIME_ROWS_AND_REPUNITS.md](PRIME_ROWS_AND_REPUNITS.md) proves that a prime target from p^a requires a prime exponent a dividing the target minus1. The collision

\[
s(32)=s(125)=31
\]

is a classical Goormaghtigh repunit collision, now read as two prime-power predecessors with different term counts. Grantham's published finite result, together with exact prime-base checks at8191, implies that31 is the only prime target below10^700 with two such predecessors. This enormous bound is **cited**, not locally recomputed. The unrestricted conjecture remains open. The paper's abc-conditional finiteness provides a genuine conditional connection, not a proof of abc.

[PERFECT_TERM_LAYERS.md](PERFECT_TERM_LAYERS.md) proves: for prime r, a perfect source with exactly2r-1 proper divisors exists iff2^r-1 is prime, and is then unique. Thus the familiar sources6,28,496,8128 occupy layers3,5,9,13. At r=11,2047=23*89 blocks layer21; at r=13,33550336 occupies layer25.

There is also no perfect seven-term source, by the complete shape proof. **Seven and21 are excluded proper-divisor term counts here, not excluded tournament statistics.** Their mechanisms differ, and no equivalence with the tournament Hamiltonian-path exclusions is supplied. Other term counts are excluded too.

## 6. What remains from the1001 anchor

The earlier [Session82 characterization](../session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md) remains inherited: among squarefree products of three primes,1001 uniquely feeds a prime cube whose next aliquot edge has square sum. The chain is1001→343→57, with343+57=400. The unrestricted fiber of343 has19 sources; squarefree-three-prime uniqueness must not be promoted to unrestricted uniqueness.

[LATTICE_COMPANIONS.md](LATTICE_COMPANIONS.md) gives a base7 carry certificate for the seven proper divisors of1001 summing to7³. It also shows why the source Boolean cube and the successor divisor chain are not isomorphic.

The other square edge1001+1001000=1001² belongs to the universal companion C(n)=n(n-1). With t=s(n),u=s(n-1),

\[
s(C(n))-C(t)=t(n-t)+(n+t)u.
\]

This gives exact noncommutation obstructions. Square adjacency sometimes survives under s, as264 and69432 map to456 and120648 with total348², but it does not survive universally. The earlier orbit of1000001 through21→11→1 remains verified background; this session does not extend the stopped1001000 factorization frontier.

## 7. Reproduction, audits, and limits

Run the eleven standalone checkers in this directory, normally and with `python3 -O`:

- `check_divisor_term_classes.py`
- `check_square_excess.py`
- `check_lattices.py`
- `check_shape_census.py`
- `check_prime_rows.py`
- `check_source_density.py`
- `check_perfect_layers.py`
- `check_fixed_excess.py`
- `check_seven_term_abundance.py`
- `check_low_term_atlas.py`
- `check_carrier_144.py`

Some emit local deterministic JSON receipts; no receipt is needed as input. `check_shape_census.py` emits the readable census table. Normal and optimized outputs were compared; explicit checks remain active under optimization. Independent audits used different methods, including a2,000,000-source multiplicative sieve for the low-layer atlas and a separate trial-division prime list for the356-prime obstruction. The source and target universes are specified in the individual notes.

This packet makes no claim to solve prime untouchability, odd perfect-number existence, general fixed-excess classification, the Goormaghtigh conjecture, quadratic-prime infinitude, or aliquot-orbit termination. The assertion that2 and5 are the only prime untouchables remains a conjectural target, not a conclusion of these class decompositions. Elementary results are presented with proofs without asserting literature priority; inherited and published results are credited separately.
