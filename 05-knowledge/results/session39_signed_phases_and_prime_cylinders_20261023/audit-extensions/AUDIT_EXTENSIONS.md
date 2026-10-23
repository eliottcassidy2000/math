# Session39 extension audit

Status: PROVED scoped mathematical checks; FINITE-EXACT independent controls.
This extends, and does not modify, the first audit checkpoint.

## Closing a prefix to a target

The closed-target theorem passes. For a prefix of length r, carry B and
valuation sum A, append b and set C=2^A+3B. For y odd and coprime to3,
the congruence 2^(A+b)y=epsilon C modulo3^(r+1) has a unique b class
modulo2*3^r, since C is a3-adic unit and2 generates that unit group.
The resulting starting numerator divided by3^(r+1) is odd integral and
positive for sufficiently large b. The already audited prefix-integrality
and backwards-oddness argument authenticates the exact extended valuation
word, rather than merely a rational branch composition. Images divisible
by3 are impossible, explaining the target restriction. Empty prefixes
are covered by the same argument.

Choosing minus targets1,5,17 puts starts from each known basin in every
finite valuation cylinder. Therefore no always-terminating valuation-only
observer can classify these basins, even on the promised union of those
basins. This does not establish convergence, complete cycle classification,
a density lower bound, or undecidability for algorithms allowed the original
integer. Combining prime ubiquity with basin ubiquity does not establish
prime starts in their intersections; the target construction is exponential,
not a Dirichlet arithmetic progression.

Independent check_closed_targets.py constructs backwards and tests
integer divisibility at every step; it does not import the producer or
use its total numerator congruence as the construction. All510 cases pass:
lengths0..3, alphabet1..4, both signs, targets1,5,17. The producer's1210
phases/3630trajectory receipt also replays exactly. Its cited inverse-orbit
precedent is credited, with no priority claim.

## Marked affine cells

AFFINE_CONNECTIONS.md and GIRTH_PHASE_HOSTILES.md were read completely.
The following claims pass:

- Every equitable singleton/triple tournament partition is in the exact
  uniform/permutation/complement carrier; singleton-induced arcs and uniform
  orientations must remain included.
- Rotation gauge retains edge signs. Joint based holonomies are classified
  modulo one common root translation, not independent cycle conjugacies.
- Point sections are common fixed points. Full unrestricted frames require
  identity holonomies; orientation-preserving frames also require positive
  signs on all individual edges.
- Girth at least five gives phase-independent vertexwise exact E in both
  orientations. The proof correctly uses no triangles and no double common
  neighbor, not the false assertion that radius-two balls are trees.
- The cycle family is strong, intrinsically cell-colored by distinct scores,
  self-converse, and has exactly C3 automorphisms. All three phases are
  nonisomorphic. For girth at least five they have the same exact E.
- Forest gauges give local phase blindness. Zero/nonzero first separate
  through some cell-marked motif at g vertices; the two nonzero phases first
  separate through some cell-marked motif at g+1. These are not asserted as
  first unmarked-motif thresholds.
- One actual mark, with intrinsic cell colors supplied, propagates to a
  discrete directed neighbor refinement in at most1+eccentricity rounds.
  This does not identify lack of a section with failure of reconstruction.

The first audit's independent graph and affine checks support the central
examples. The final expanded producer is replayed against its frozen
receipt, including all512 binary3-by-3 block recognition controls.

## Pins and reproduction

PRODUCER_PINS.json specifies the precise approved bytes. Run
python3 check_closed_targets.py and compare with closed_target_independent.json.
The producer receipt replay output is also saved as closed_target_producer.json.
For the marked producer, run its check_marked_connections.py without --write;
it checks its own frozen receipt without changing it. Scope, sources and
non-implications remain those declared in the approved prose.

## Native odd-depth coverage boundary

COVERAGE_CLOCK_BOUNDARY.md passes mathematical audit. For both signs,
U(x)<=(3x+1)/2<=2x, so the i-th state from n<=X is at most2^iX.
The valuation letter bound floor(log2(3*2^iX+1)) follows directly.
A fixed word and target determine at most one starting integer via its
nonzero affine slope. Counting words, including the zero-length word,
gives the asserted product-sum upper bound. Every letter is bounded by
floor(log2X)+r+3, validating the coarser expression uniformly in y.
For r=o(logX/loglogX), its logarithm is o(logX), hence coverage is X^o(1).
This is a native accelerated-iteration horizon obstruction, not algorithmic
complexity or a statement about unbounded-depth annulus budgets. A fixed
finite target set contributes its cardinality factor; moving target sets
need that factor explicitly. The producer's180coverage cases and40920
step-height/valuation controls replay exactly. The asymptotic conclusion
is proved, not inferred from that small finite test.
