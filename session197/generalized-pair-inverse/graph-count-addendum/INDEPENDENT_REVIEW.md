# Independent audit: generalized pair graph-count limitation

**PASS, no correction.** The presentation count is a valid upper bound on the unrooted graph image. It is smaller than the complete C2 class count at order14, and asymptotically negligible relative to that class. This has no negative implication for scalar Hamilton-path coverage.

## Frozen binding

- `session197/generalized-pair-count-bound/PROOF.md`: `45b6dd975497a9feca61e28677d3477c7e8310efe612ca13c0720e7064c79440`.
- verify.py: `0646191f213a9077ecbe3a40040441b2f824301f21e15e2d2e19465ba69f1c84`.
- VALIDATION.json: `4ee36a82a71c67f7ae63fed496b2f1022064fb16a52753fd47f580c2664411da`.

The source was copied before replay because its verifier writes the receipt. Normal and optimized runs reproduce the frozen bytes. No original file was changed.

The exact named dependencies were read and their hashes checked:

- generalized-pair-inverse/PROOF.md: `9df4e06c88d28e159953601611bbcc6fb90f00f6766b5c09b03b3dfc7c76cbd3`.
- FIBONACCI_CLASS.md: `d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e`.
- C2_COMPLETE.md: `4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d`.

Only the actual family definitions, source domains and vertex increments are needed from the generalized-family packet. Its scalar identities and inverse algorithm undergo a separate review and are not used to derive this counting obstruction.

## Structural checks

The rooted strong-source count is B=z+z³/Delta, with the singleton included and no strong order-two source. The earlier Fibonacci packet's caution about an unrooted bridge subclass does not weaken its rooted start-one enumeration. For a weak start-one graph, the marked root lies in the initial SCC and every later SCC is a singleton. That intrinsic decomposition identifies both the rooted core and the tail length, so W=B/(1-z) counts rooted isomorphism classes exactly. Algebraically W=z(1-2z)/Delta.

A suffix of positive length m has an arbitrary first m-1 binary digits and a final1. Including the empty suffix gives S=(1-z)/(1-2z). The ordinary gadget adds its root as well, giving zS. The literal generalized-family definitions give these increments:

- N: two gate vertices, then a nonempty suffix; source order at least2.
- M: three gate vertices, then any suffix; source order at least2.
- Q: two gate vertices and one ordinary gadget; source order at least2.
- P: four gate vertices, then any suffix; singleton source permitted.

Thus the presentation generating function is exactly(W-z)[z²(S-1)+2z³S]+Wz⁴S, simplifying to4z⁵(1-z)²/[Delta(1-2z)]. This counts presentations, not distinct outputs. A rooted source isomorphism extends over the deterministic named gate and suffix, so no further vertex-label or automorphism multiplicity is required. Conversely multiple presentations may give one graph; this only strengthens the upper bound. Forgetting the output root likewise cannot increase the graph count.

The complete C2 generating function is precisely the independently proved z⁵+z⁵(1-z)(1-z²)/Delta², including the regular-five exception. At order14 its count is44,394. The family-presentation contributions are N16,687, M10,434, P4,181, Q10,434, totaling41,736. Therefore at least2,658 C2 graph classes cannot admit one of these presentations, regardless of any finite search.

At alpha=(3-sqrt5)/2, the presentation function has a simple pole with nonzero numerator; its other poles have larger modulus. The complete C2 function has a double pole with positive nonzero leading coefficient. Partial fractions give O(rho^n) versus Theta(n rho^n), with rho=1/alpha. Their ratio tends zero, proving missing graph classes at every sufficiently large order. No cancellation removes either relevant pole.

This invalidates a graph-isomorphism-preserving reduction of all C2 graphs to one generalized gate. A graph outside the family may share its H with a graph inside, at the same or another order. Scalar-only absorption, regular target ranges and cofinite scalar coverage remain untouched.

## Independent arithmetic controls

The standalone `check.py` does not import the source. It derives B from the canonical rooted decomposition recurrence rather than the Fibonacci closed form, adds the weak tails, and separately accumulates each of the four gate contributions. Multiplying the resulting coefficient array by its claimed denominator recovers exactly the numerator4z⁵-8z⁶+4z⁷ through order300. A separate convolution of two reciprocal-Delta sequences gives the C2 count. The order14 values and the continued finite gaps through300 agree.

- check.py: `ceaa2bf9b5fdd50a578303787a4a1197a29a4f70d4bba3c3363da5923f7efa9c`.
- independent-normal.json and independent-optimized.json: `3a662eeb6a7bf05fd19c7d6e464cd93eebbb14f15319d66999f241e76194678d`.

Run `python3 check.py` and `python3 -O check.py`. Both are standalone, do not write files, and agree exactly. The finite coefficient controls support the implementation; the all-order conclusion follows from the exact presentation count and the pole comparison.
