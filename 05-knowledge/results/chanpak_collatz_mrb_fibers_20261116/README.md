# One arrow, exact fibers, and the absolute remainder

Session 2026-11-16. This packet combines the decimal Collatz identity, finite tournament count fibers, and the MRB idea of bounding an entire weighted tail. The result is a sharper distinction between a programmable residue and a constructive exact target. **Neither Collatz nor Chan–Pak is settled.**

## The real content of 1024 = 1000 + 24

In exponent coordinates the odd Collatz arrow v→3v+1 is exactly the map y→2y³ on y=2^v; the even arrow is a square root. This embeds the full functional graph, with exactly one outgoing arc, into the positive powers of two. The decimal rewrite holds for just **1→4 and 3→10**: even the unconstrained equation 2^u=10^v+3·2^w has only the two positive solutions (4,1,1),(10,3,3). The elementary proof uses a fixed-gap (Pillai-type) equation and the Fermat factor 641. This is not an application of [Mihăilescu’s gap-one Catalan theorem](https://www.math.kth.se/4ecm/program/mihailescu.pdf).

There is also an exact Collatz inverse lesson: fix a parity word and a positive endpoint. The affine inverse is a valid positive source whenever its integrality condition holds, because nonpositive starts stay nonpositive under the parity rule. In a tournament decoder the analogue would be a proven, physically realizable source region whose complement cannot enter it under any legal digit. [Tower proof and inverse certificate](COLLATZ_TOWERS.md).

## Carmichael and actual fibers: a scalar is not the whole vertex

The [totient image question and inverse-multiplicity question](https://arxiv.org/abs/math/9907204) are distinct. For tournaments, ordered strong components give an exact polynomial-in-order formula for each fixed H once the order is sufficiently large, with coefficients separating true strong multiplicity from singleton padding. The same method handles the converse quotient. It also isolates Chan–Pak's exact remaining atomic tasks after explicit finite seeds 49,63,343: strong witnesses for the odd primes other than 7, and for 7p with p an odd prime other than 3,7. The proof targets those atoms; it does not infer that they exist. [Multiplicity formula and strong finite seeds](SCC_MULTIPLICITY.md).

There is a literal counterexample to using H alone as the return-graph state: two equal-height states with H=33,423 share the same rank-four image. Their 136 distinct exact one-block successors on either side are **disjoint**, though both sets cover every odd residue mod 32. The exact lattice of future readouts is diag(4,2,2,1), so the entire modular continuation language can nevertheless be tracked with just **32,768** reachable states instead of 524,288. [Equal-count fiber and minimal observable quotient](FUNCTIONAL_FIBERS.md).

## The MRB idea improves the absolute obstruction

The older MRB contour audit bounds a whole weighted tail, not a single attractive term. The tournament version weights **all** isolated return triples by H^(-s). Both entire triple sums are certified below 1 by directed integer rounding. Retaining the common *live suffix* sharpens this further: the absolute spectra are **o(X^(25/36))** for the 111-action mixed-prime bank and **o(X^(23/30))** for all 168 words. Even ordinal-sum closure stays this sparse. [Live-suffix pressure and exact checks](SUFFIX_PRESSURE.md); [isolated-tail proof, MRB source and changing-core diagnostic](MRB_PRESSURE.md).

The practical synthesis: use a sufficient state to make each move predictable, and require a physical certificate for the inverse (not just that its final residue agrees). Then check the weighted growth of the whole proposed language. The present bank has exact decimal steering and recovers its image, but its weight is still too small to represent every absolute target.
