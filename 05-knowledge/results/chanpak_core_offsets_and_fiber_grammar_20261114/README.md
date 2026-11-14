# Changing offset cores and exact Hamiltonian-fiber grammar

Two different recursions are needed: one chooses the actual tournament for a prescribed absolute count; the other counts all isomorphism classes without discarding their continuation data.

## Physical construction

**Proved.** At every integer height t >= 1, strong tournaments realize the eleven consecutive odd counts **13790·2^t + 885,887,…,905**, and the nine consecutive odd counts **288·4^t − 11, −9,…,+5**. The core and contacts are selected by the residual and genuinely change; in the affine bank there are nine values of the core Hamiltonian count and ten values of its ordered pair-path count. At height 1 the affine interval is 28465…28485. The arithmetic inverse tests for an exact power of 2 or 4 and returns a literal tournament. A finite union of such rays is sparse; universal placement remains open.

A different mechanism gives an exact two-count return at every pair of independent wing heights. The full-interface banks of weights 65,109,161 have unit moves whose large effects cancel because **65·109+65·161−109·161=1**. Arbitrarily marked transitive spacers preserve the difference, not the absolute starting count. Repetition of the older 9/45 banks instead has an interval-dependent 3-adic obstruction.

## Recursive spectrum and the deletion metagraph

**Proved.** Treat each class as its unique ordered word of strongly connected components. Each pointed strong class C supplies one local rule per vertex orbit: replace the contiguous component word of C−v by C. Singleton insertion is the only equal-H move, so neutral growth for a fixed nontrivial backbone is a Pascal lattice in the lengths of its gaps. Fixed positive count jumps are much more rigid: a jump J can involve a marked strong block of order at most J/m+2, where m is the exterior Hamiltonian-count product and an odd divisor of J/2. Consequently every fixed positive one-vertex jump has a finite set of absolute endpoint pairs, independent of height. Jump two occurs only at 1→3 and 3→5. This does not constrain two unrelated tournaments to have distant counts.

The exact ordinary class series is the inherited ordered-component series V=1/(1−S), with multiplicative count labels. Its new pointed lift is E=(Delta V)^2 D; its fixed part under converse is Ec=(Sigma V) Dc, with Sigma doubling the order and squaring the exterior count. Ordinary isoclass multiplicities, pointed edge multiplicities, and the simple metagraph are distinct quotients. The positive inventory correction stratifies the class series by number of unavailable large components; all odd fibers through 73 at all orders can be counted from the strong order-eight catalog. Forgetting count labels gives A000568, but the recurrence does not remove the need for new strong atoms at higher counts.

**Open:** the universal absolute selector, a complete simpler description of the unrestricted strong alphabet, and Chan–Pak. The exact windows are not a cutoff argument.

## Proofs and replay

- [Changed-core windows and exact inverse](WINDOWS.md), [literal decoder](window.py): forty tournament recounts over the two heights.
- [Three-bank cancellation and the 3-adic obstruction](ORDINAL_SELECTOR.md), [independent endpoint transfer](ordinal_controls.py).
- [Pointed SCC rewrite and converse recurrence](METAGRAPH.md), [complete small-layer comparison](replay.py).
- [Positive catalog corrections](POSITIVE_CORRECTIONS.md).
- [Bounded-gap rigidity for arbitrary induced inclusion](INDUCED_GAPS.md). At unbounded absolute height even many vertex deletions cannot connect consecutive odd counts.

The strong affine witnesses have at most t+14 vertices; the strong quadratic ones have at most 2t+10. Adding neutral singleton SCCs gives a witness at every larger desired order without altering its count (but not preserving strength).

## How the spectrum recurses inside strong classes

The [profile-refined modular census](MODULE_PROFILE_THEOREM.md) lifts the canonical prime-quotient Burnside recurrence to actual H fibers and independently reproduces the complete profiles through order six ([replay](module_profile_verify.py), [receipt](module_profile_verification.json), [audit](MODULE_PROFILE_AUDIT.md)). One isotype is chosen per automorphism cycle, not one per vertex of the cycle. Full path-cover polynomials are exactly the observational states for whole-module substitution; H alone already fails at order five. At fixed H=h their support has the sharp bound 5^(2j)<=h^3 on coefficient j. This groups the 56 order-six classes into 24 profiles but does not remove the prime-quotient input. A different marked sidecar, the [triangle dimer formula](TRIANGLE_DIMERS.md), gives exact weights in one infinite profile family and shows why many physical profiles there still give no usable offset.
