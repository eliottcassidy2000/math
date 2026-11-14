# Session 325: absolute target overlap and the full 99-offset bank

**Chan–Pak remains open.** This session proves exact common targets and a full offset alphabet, not just congruence alignment. It also identifies why fixed banks cannot be turned into a universal absolute inverse by padding.

## New exact coverage

[The full physical offset alphabet](FULL_99_DIGIT.md) strongly realizes every number 99·2^m+c for m≥1 and c=1,3,…,197. These are **99 consecutive odd targets at every scale**. Every digit has a verified literal rail with a phase shift of at most four; the finite initial bridge consists of 115 literal strong witnesses. At the specific target 99·2^m+23, two distinct stars and a unit-end Busch seam meet at every scale. The two stars agree for arbitrary input profiles, not just transitive reservoirs.

For a general potential source with H coprime to 99, the full offset bank determines its unique odd digit modulo 198 and forces the exact two-slot count Q=(N−cH)/99. Singleton phase cost is explicit. Fixed H can supply only a finite union of dyadic Q rays. That is the remaining constructive obligation: choose a growing **physical** source profile (H,Q), not merely a congruent one.

## Independent absolute intersections and state control

- [Primitive selective-chart overlap](ABSOLUTE_OVERLAP.md): two different strong contexts realize exactly 396561 Q_V(2)−571 H(V) for every nonempty V. The common family contains the prime 792551 and 7·3625619, with a prime cofactor. [Independent review](SELECTIVE_REVIEW.md) and direct order-19 replays are supplied.
- [Ready overlaps](READY_OVERLAPS.md), [the full binary-cut corridor](FULL_CUT_CORRIDOR.md), and [the exact branching inverse](BINARY_SEAM_INVERSE.md) give full-state coalescence of nonconverse cores, a same-size asynchronous overlap, and a common active interface under which three consecutive odd labels stay consecutive. The broader corridor still has a strict absolute density cap. A target can now be tested against a ready corridor rather than assumed to belong to it.
- [The seven rail and hidden order](RAIL_OVERLAP_AND_HIDDEN_ORDER.md): strong tournaments at the same absolute H can have an entire exact corridor of marked start counts. The exponent decoder reconstructs the singleton placement, while the arithmetic 3-chains for different dyadic exponents remain disjoint.

## Structural limits

[Exact exponential-chart overlap](EXPONENTIAL_OVERLAP.md) gives affine-lattice descriptions (hence ultimately periodic exponent sets) for finite exponential charts. Pure dyadic overlap reduces effectively to a bounded gap test and linear equations. Two fixed rays with different offsets meet at most twice. The two all-transitive selective charts are exactly disjoint.

[The asynchronous growth obstruction](FINITE_GROWTH_OBSTRUCTION.md) shows that after the primitive selective collision, even independently growing both source lengths and padding different reservoirs drawn from a finite bank can work at only finitely many pairs of source lengths. The safe continuation is arbitrary common reservoir growth at fixed source state; unbounded source growth needs new reservoirs or a new interface.

The lateral source pass used reproducible random sampling rather than only earlier Chan–Pak papers; the [historical bridge](RAIL_OVERLAP_AND_HIDDEN_ORDER.md#1-a-reproducible-lateral-source-pass) names the sampled Fourier/Zeckendorf, Smith, Mazur, and planar Jacobian mechanisms and what each actually transfers. See also the orientation correction to S324 in [the digit proof](FULL_99_DIGIT.md#3-three-genuinely-different-cuts-at-992m23).

Run the `replay_*.py` files next to these notes with Python 3. For the large duplicate prime-cofactor check use `python3 replay_overlap.py --full`. The bank is an exact certificate, not an exhaustive search claim.
