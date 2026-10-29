# Independent structural audit — Session 107

Date: 2026-10-29. Scope: character/fiber argument, squareclass constraints, Monster valuation and Paley quotient, and arithmetic orientation guards. No publication by this auditor.

## Verdict

**PASS** for the reviewed character and group proofs, subject to their explicit fixed-root / fixed-modulus / finite-tournament scopes. Independent exact tests pass. No new classification, infinite Sylvester-prime claim, or common source of Wieferich lifting and Monster primes is justified or asserted.

## Character proof details checked

1. The quadratic large-sieve restriction gives an operator norm, because arbitrary complex coefficient vectors are allowed. Taking a square root gives `B_X = O_epsilon(X^(1/2+epsilon))`, after renaming epsilon. This is `o(N)` for fixed epsilon less than 1/2.
2. Telescoping each edge of a simple labelled tournament is valid even with intersecting fibers. After fixing the other vertices, all surviving factors separate as `c f(x_i)g(x_j)`, with absolute values at most one. A multiple-edge or genuinely higher-order observable would need a different argument. The stated sum `B_X/2 sum (n_i n_j)^(-1/2)` follows.
3. The conditional collision bound is valid, including its constant. To spell out an otherwise terse step, write `c=P(collision)`, `t=2^(-m)`, `Z=product W`. Then

   `|P(F | distinct)-t| <= (|E Z-t| + |E[Z 1_collision]-tc|)/(1-c)`.

   Since `0<=Z<=1`, the second numerator term is at most `max(t,1-t)c<=c`. Using `c<=D<1` proves the displayed `(A+D)/(1-D)` bound. One should not instead separately charge both an unconditioned collision error and another full conditioning error.
4. The squareclass image is affine, not generally linear: `y=s+Vz`. On q=3 mod4, -1 is fixed negative; 2 remains a free bit controlled by q mod8. Odd prime bits are independent CRT coordinates after the reciprocity signs are retained. PNT in each **fixed** reduced progression yields equal asymptotic color densities. Dirichlet alone would give existence, not density.
5. Square 4 is identically positive. Repeated 2 gives a duplicate bit. The named prime bases 1093 and13 are not carrier vertices of the 3 mod4 prime tournament. Primes dividing bases are excluded to avoid zero symbols. No growing-base or growing-modulus conclusion follows.
6. The sparse Sylvester prime-term sequence is a valid hostile to an unsupported arbitrary-subsequence inference: a later term is 1 modulo every earlier term; among the prime terms congruent 3 mod4, reciprocity makes all edges coherently ordered. The proof does not require infinitely many such terms.

## Group structure: independently checked primary sources

The official [ATLAS Monster page](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) and [Baby Monster page](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/) verify the order exponents used in the valuation chain. Direct subtraction produces `13→2→5→3→7→4`, not a literal prime cycle. Identifying13 and4 is extra quotient data; the centralizer normalization by2 changes that map.

The [Breuer–Magaard–Wilson computational supplement](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf), §4.13 (printed p20) and §5 (p22), explicitly identifies the 7A normalizer as a subdirect product of `7:6` and `He.2`. Thus quotienting by its normal Held subgroup yields `7:6`; the orientation-preserving degree-seven action is `7:3`. Section2 also verifies the stated restriction decomposition of the196883-dimensional character and distinguishes the central-faithful96256 constituent. These imported character-table results were checked in the primary source, not recomputed here.

The elementary obstruction is sound: an involution in a tournament automorphism group would swap two vertices and reverse their unique arc. Hence its automorphism group has odd order. A nontrivial action of an even-order simple group would be faithful, contradiction. The local quotient does not violate this: only7:3 preserves orientation, while the extra coset reverses all arcs. Both local affine groups have nontrivial cyclic quotients and therefore minimum complex representation degree1. This does not contradict spectral cancellation for a specially chosen Paley connection set.

## Independent exact verification

`python3 check_independent.py` checks CRT signatures using a separately implemented Jacobi algorithm on **all reduced classes** of each tested modulus, reciprocity and supplemental laws for94 odd primes, all5040 Paley7 permutations, the order-ratio arithmetic, character-degree sums, Sylvester product recurrence, paired13/139 orientation signs, and the opposite quadratic signs of the two named Wieferich primes. It uses runtime checks rather than removable assertions.

These finite tests guard arithmetic and conventions; the unbounded statements rest on the proofs and cited analytic/group inputs. No entire Monster character-table reconstruction was attempted. No denied tool action or external mutation occurred.

## Euclid-offset conditioning audit

Also reviewed `conditioning/EUCLID_OFFSET_SYNDROMES.md`. The multiplicative character syndrome is exact for arbitrary exponents; the claimed rank loss correctly requires a nontrivial homogeneous relation, not merely a nonzero affine offset. The Sylvester transitive-subtournament deduction is valid for **odd prime terms** (indices at least1). An initial wording ambiguity that could include s0=2 was reported to the parent and repaired before final hashing. The double-exponential term lower bound and O(log log X) sparsity argument are valid and do not assert sparsity of the union of all prime divisors of arbitrarily late terms.

For Fermat factors, exact order `2^(i+1)` implies p=1 mod8 when i>=2; the negative2-character rows at roots3 and5 are correctly retained. The three displayed later Fermat primes give a symmetric residue clique, not a tournament. Independent tests additionally verify primality of3263443, the displayed s6 factorization, all six earlier-root zero syndromes, F5's displayed factorization and five root syndromes, and the residue clique.

The Gowers normalized mixing inequality was checked against the Cauchy–Schwarz bound in the proof of Theorem3.3 (printed p12); the source supports the absolute-error formulation. Aricheta Table1 explicitly identifies the odd Baby primes with the rationality property at level X_0^+(2). The Monster minimum196883 remains imported ordinary-character-table data, not an independently recomputed minimum; the cited2024 verification paper verifies that full table. The Baby minimum4371 additionally follows directly from the low-degree list and central-involution values in the computational supplement. No modular minimum is substituted for either complex degree.
