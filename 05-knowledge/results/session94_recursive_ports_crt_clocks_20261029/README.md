# Recursive port inventories, sparse CRT certificates, and source-preserving clocks

2026-10-29. Initial source snapshot `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.
The session develops full recurrences for specified observations, a sharper
re-marking protocol, and an exact arithmetic coupling to critical ternary
branching. Classical species, branching and classification inputs are credited.
No external-priority, new A000568 record, or Collatz-convergence claim is made.

## 1. Shared ports require cycle data beyond every fixed-point moment

[Weighted cycle carriers](ports/WEIGHTED_CYCLE_CARRIERS.md) give full recursive
cycle indices for odd cyclic tree actions. In the ternary tournament family,

    Z_(A_(k+1)) = (Z_(A_k)^3 + 2 Z_(A_k)[p_i→p_(3i)]) / 3,
    Z_(C_k) = p_1 product_(j<k) Z_(A_j)[p_i→p_i²].

The second substitution doubles **cycle multiplicities**, not cycle lengths.
It retains the same automorphism acting on both mirrored copies. General odd
cyclic alphabets retain every divisor of the root rotation order; for nonternary
alphabets the note states the permutation-action scope explicitly.

Weighted palettes are substituted into these indices, with the paired colors
forming a Cartesian-product palette. At nine vertices, coupled binary palettes
give192 classes; independently quotienting the two sides incorrectly gives128.

Even the entire fixed-point-moment hierarchy does not determine arbitrary
attachments. The regular actions C9 and C3×C3 have identical fixed-point laws
and identical all-order named-port moments, yet their three-subset orbit counts
are10 and12. The full cycle inventory is a necessary distinction.

For unordered r-port subsets in the ternary family, the exact count U(k,r) is
a degree-r polynomial in depth k, leading coefficient2^r. A saturated LCA-depth
construction gives nonnegative integer coefficients

    U(k,r)=sum_(ell<=r) beta(r,ell) binom(k,ell).

This yields O(R³) arithmetic preprocessing and O(R²) evaluation of all counts
through r=R at any depth, rather than stepping through every level. It was
tested at depth one million; integer bit costs remain charged. Examples:

    U(k,1)=2k+1,
    U(k,2)=k(4k−1),
    U(k,3)=k(2k−1)(4k−3).

The last requires a three-cycle correction; it cannot be obtained by treating
fixed vertices as the only invariant subsets. This is a recurrence for this
specified family, not unrestricted A000568.

## 2. Recursive A000568-related observations without a prime-class bank

[Observation jets](enumeration/RECURSIVE_OBSERVATION_JETS.md) recover three
aggregate inventories from the existing colored Burnside channels:

1. Ordered named ports, with repetitions, and how many strong components they
   occupy. The component decomposition is canonically ordered. For two ports,
   the same-component series is R2−2R1²/A.
2. Distinct named children of a prime-root decomposition. A triangular
   derivative recurrence extracts these from a small bank of scalar colored
   channels, rather than constructing the full aggregate prime cycle index.
3. Selected root-cycle derivative channels, using a divisor-triangular
   recurrence and gcd-weighted cycle data. Ordinary fixed-point moments alone
   do not provide these channels.

At fixed derivative order J, the child-slot compiler uses J+1 inherited scalar
channels and polynomial formal-series algebra. Its scope is the requested
observations: it does not recover prime representatives, individual groups,
transporters or a complete metagraph bank.

The order11/12 tables retain that distinction. The eleventh-cycle derivative
has constant32/11, matching the inherited order11 symmetry ledger; the
ninth-cycle derivative has degree-three coefficient16/3, matching the relevant
order12 C9 contribution. These are normalized cycle-index coefficients, not
fractional tournament counts. Literal low-order module/port controls and a
separate full-index comparison verify the extraction.

A000568 and strong-count prefixes already reached350 in the repository. This
session supplies recursive query compilers, not a new coefficient record or
a claimed measured speedup over the best scalar implementation.

## 3. A much shorter legal loop for every retained correction

[The portrait protocol](loop_protocol/LINEAR_LOOP_PROTOCOL.md) improves the
inherited3[1+(k−1)3^k] re-marking bound to

    3(3^k−2k−1)/4.

At27 vertices this is15 steps instead of165. Every paired internal-node
rotation is one conjugated based triangle. A descending-depth portrait
normal form uses each nonzero rotation once, including its inverse, and
retains the exact action on all named ports.

This counts actual square-root re-marking edges. The conjugator constructs
the triangle's endpoint marks; it is not silently performed as a free move.
Computing and storing those full marks still costs resources. The bound is
not claimed optimal and is not a bound on black native-complement moves.

For a uniform centralizer source, c=(3^k−2k−1)/4 portrait entries are independent
uniform trits. The prescribed protocol cost is3 times Binomial(c,2/3), with
mean2c. A different source law must retain its own weights.

## 4. The critical ternary law has an exact Collatz-source realization

[Valuation trees and clocks](dynamics/VALUATION_TREE_CLOCKS.md) establish a
specific coupling, not a dynamical analogy. On the correctly punctured odd
2-adic domain with Haar source measure, successive accelerated valuations are
independent with probability P(a=j)=2^(-j). Their exact factorization is

    a = 1+b+2G,
    b ~ Bernoulli(1/3), G ~ geometric-on-{0,1,...}(3/4), independently.

Give a tree vertex three children when b=1 and zero otherwise. This is exactly
the critical ternary law governing the fixed-point genealogy of uniform
recursive automorphisms. It does not specify a full uniform automorphism or
its holonomy: that requires additional phase/source information.

Retain the valuation at every breadth-first vertex. If F(u,z) records total
population and division clock, its exact recursive law is

    F = u(2z + z² F³)/(4−z²).

Conditional on total population N, the division clock has mean2N−1/3 and
variance16N/9. The standard Collatz clock also includes its odd steps. The
ordered residual labels and the unused source tail are retained when an exact
source is required; aggregate tree size and total clock do not determine
chronological affine carry.

The auxiliary tree stops for a natural-density-one set of positive odd
sources. **That is not Collatz descent or convergence.** On the explicit
cylinder x≡41 mod64, the tree stops after four accelerated steps but the
endpoint is(81x+103)/32>x. The seed trajectory is41→31→47→71→107. Conversely,
seed1 has an infinite auxiliary tree although its native orbit is the familiar
root cycle.

The11-cell history supplies a further rooted hostile: the rational period-five
cycle starting1/13 has valuation word(4) repeated; the integer cycle can start
at−5 or−7, producing(1,2) or(2,1) repeated. The auxiliary tree can die or survive
after changing the root of the same arithmetic cycle. The map is not an orbit
conjugacy that forgets the named source.

Negation relates3n+1 and3n−1 on the stated domains, but not positive-domain
convergence. Antihydra retains a different counter barrier. A stale whole-space
coding sentence was repaired by propagating the **already-proved Session64**
punctured-domain theorem; see [correction provenance](dynamics/DOMAIN_REPAIR_PROVENANCE.md).
No shared ledger was changed.

## 5. Compatible CRT charts beyond76, with a nine-candidate inverse

[The CRT carrier theorem](crt/COMPATIBLE_CRT_CARRIERS.md) gives, for odd n,

    Z[phi]/(phi^n−1) ≅ Z/L_n             if3 does not divide n,
                        F4 × Z/(L_n/4) if3 divides n.

Odd factors occur on a selected split branch. Prime powers retain their
precision: at period55 the Z/121 chart uses phi=37 of order55, while reduction
to11 has order5. Swapping in the conjugate would change the clock.

For a ring projection pi and only rotation/reflection, the **smallest closed
observer** is the reflected double(pi(x),pi(Rx)). Rotation updates its two
coordinates by phi and phi^(-1), and reflection swaps them. It need not recover
the whole source. Adding arithmetic translation+1 is a different menu and can
force a much larger refinement.

[The inverse reader and sparse certificates](crt/NINE_CARRIES_AND_ALTERNATING_GROUPS.md)
make that boundary concrete. Every odd-period ring state can be decoded by
at most nine integer carry candidates and exact greedy golden extraction.
The proof uses both real embeddings and charges O(n) algebraic operations,
not constant-cost integer arithmetic.

This produced complete sparse primitivity certificates at periods25 and55.
The latter has312,119,004,989 states, but its independent verifier checks only
44 legal word occurrences, plus factorization and all22 possible nontrivial
divisor block systems. It does not import the decoder or enumerate the state
space. Period25 uses six word pairs.

Jones's cited CFSG-dependent full-cycle classification, after exact exclusion
of every projective degree and checking generator signs, then gives

    G25 = A167761,
    G55 = A312119004989.

These are two finite-specific composite-degree results, not universal
composite primitivity. They illustrate the difference between a huge exact
operation-closed carrier and a short checkable certificate about it.

## 6. What the recurrences unify, and what they do not

The common mechanism is to choose the required observation first, then keep
the data its composition law actually uses:

- cycle multiplicities for weighted attachments;
- derivative channels for aggregate decomposition queries;
- paired portraits for legal marking corrections;
- source valuations and chronological residuals for arithmetic trees;
- reflected charts and exact selected roots for CRT clocks.

Equal cardinalities, matching moments or a short recurrence do not by themselves
identify these carriers. The regular C9/C3² hostile and the growing-but-stopped
Collatz cylinder make those limits explicit.

[Independent audit](audit/AUDIT.md) records proof reviews and separate finite
implementations. All released checkers use standard-library exact arithmetic,
run normally and with `python3 -O`. Two small JSON files are mathematical
primitivity certificates required by their independent verifiers; bulky run
receipts are not required inputs. The known all-odd reader and Jones theorem
remain explicit dependencies.

The session does not solve positive-integer Collatz convergence, Lucas-prime
infinitude, the unrestricted prime tournament atlas, or a general metagraph
connectivity problem. It supplies full recursive observation laws, a sharper
legal protocol, and independently verifiable large composite certificates.
