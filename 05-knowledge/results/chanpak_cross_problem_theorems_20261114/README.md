# What the Hamiltonian spectrum can and cannot borrow from other problems

These are all-order results for specified tournament constructions.
No claim of Chan–Pak completeness, nor a wider finite coverage census.

1. [Prime divisors, Fermat primes, and the square-class obstruction](PRIME_WRAP.md).
   A cyclic wrap has N=h+2Q. Quadratic reciprocity proves that **even
   infinitely many sources of finite square-class rank miss infinitely
   many prime and 7p targets**. For a finite palette, multiplicative
   order gives an exact divisibility test and connects to the known
   Kummer density theory. Transitive-source prime outputs are exactly
   the Fermat primes other than 3. Prime-power solvability need not
   glue across primes: 3 and 5 divide separate transitive-wrap
   outputs but 15 never divides one.
2. [An explicit unbounded-rank source bank](FERMAT_RANK.md).
   One-return sources at dyadic widths supply independent square
   classes, regardless of whether there are infinitely many Fermat
   primes. The resulting outer-wrap bank sees density one of primes
   as divisors, not as exact outputs. Two refinements: [almost every
   fixed-length prime tuple has a simultaneous witness](SIMULTANEOUS.md),
   yet [the exact count set has a sharp subpolynomial size](BINARY_BUDGET.md).
   With arbitrary rather than Fermat-indexed widths a [five-ratio
   identity](FIVE_RATIO_IDENTITY.md) gives simultaneous unit coverage
   at any odd modulus coprime to 3·5·7·11·17 and at all individual
   odd prime powers. A genuine obstruction survives at 255 even though
   each of its two-prime factors admits a witness.
3. [A modular compiler at fixed Hamiltonian cycle count](FIXED_CYCLE_CONGRUENCES.md).
   The complete scalar residue obstruction for a fixed painted source
   is gcd(h,S). The planted cycle count remains h even when the
   congruence orbit is full. A seed with cycle count 5 and path count
   69 matches every permissible residue but hits at most 1/16 of all
   odds in upper density. This separates height control from the
   sort of finite-quotient mixing used in thin-orbit problems.
4. [From path counts to Yang–Lee zeros](ZERO_GEOMETRY.md).
   The Hamiltonian count is a hard-core partition function at
   activity two. Repeated cores inside a cyclic wrap give a
   computable zero-activity transition. Triangle cores give a
   vertical Cauchy zero distribution plus a single exponentially
   small zero.
5. [Encodings versus full conflict graphs](CONFLICT_GRAPHS.md).
   Every finite graph occurs as a chosen odd-cycle conflict graph
   of a strong tournament, but the *full* conflict graph has a
   universal vertex in each nonempty component. This is why the
   multivariate encoding cannot be read as scalar spectrum coverage.

The most useful operational distinction: **a character can reject
a bank, and an unbounded bank can remove almost every divisibility
obstruction, but the absolute target is a separate inverse problem.**
The next target is a selective or multiple-visit physical construction
whose square-class sign and positive height can be chosen together.

Run `python probe.py` and `python replay_five_ratio.py` for independent
finite controls on the formulas.
