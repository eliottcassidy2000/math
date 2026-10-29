# Packing, golden observation, and the operators a representation must retain

2026-10-29. Source snapshot `90e704a8f35614d2085778d1b518f03ff0d750ca`.
This session connects the earlier eleven-square, Zeckendorf, even-graph and
trijection work through explicit arithmetic maps, guarded operators and metric
lifts. Classical inputs are credited; new applications and finite controls
are scoped. No external-priority or major-conjecture closure is claimed.

## The synthesis

The useful unifier is not that several objects contain11,7,5 or the golden ratio.
It is the computation of **what survives a change of representation under the
actual next operation**. Three concrete answers now meet:

- The kernel of a modular observation matrix determines exact arithmetic memory.
  Its determinant is the golden norm.
- A hard feasibility restriction leaves a positive boundary Gram in the flip
  operator. Its rank is controlled by signed incidence and bipartiteness.
- The toroidal seven-point model supplies exact pairwise geodesic distances,
  but its shortest displacements need not integrate into one planar placement.
  Cycle sums detect the missing information.

These use familiar observability, incidence and quotient mechanisms. The
contribution is their exact realization on the repository's objects, including
constructive repairs and sharp failures, rather than a claim of one new theory
that resolves every problem in the corpus.

## 1. The golden norm is an exact memory invariant

The [Zeckendorf observer theorem](zeckendorf/NORM_CONTROLLED_OBSERVERS.md)
starts from the native two-register digit update. For the modular observation
`aX+bY`, the exact retained register is the image of

    [[a,b],[b,a+b]],       determinant = a²+ab−b².

Smith invariants give the minimum complete automaton size for any reachable
single residue target, including the previous-digit legality state and the
rejecting sink. The ordinary Fibonacci divisibility minimum is classical;
the note credits Charlier–Rampersad–Rigo–Waxweiler and supplies its constructive
reachability argument within the more general observer formulation.

Concrete consequences:

- Ordinary integer value modulo11 needs243 complete states.
- The golden eigen-observer X+4Y modulo11 needs23.
- Keeping4 when moving to modulus121 costs2663 states.
- Correctly lifting the eigenvalue to37 restores243 states.

The prime-power behavior is exact. Inert primes such as2,3,7 require two
registers for primitive observers; split primes such as11 permit a scalar
Hensel-lifted channel. At the ramified prime5 the minimum register size at
precision5^k is5^(2k−1), a different precision defect.

The same indefinite form appears in the golden packing clock: with
N=[[2,1],[1,−2]] and Q=[[0,1],[1,1]], QᵗNQ=−N. Thus appending a zero preserves
the modular memory index while reversing the norm's sign. This is an actual
shared form, not a claim that it is a positive Euclidean metric.

Joint moduli also give a positive result: compatibility of BOTH registers and
ONE shared legality bit recovers exactly2 lcm(m,n)² live states. Separate
marginal state counts alone would miss this compatibility condition.

## 2. The cyclic golden reader needs an operation-specific repair

The inherited all-period theorem identifies cyclic no-11 words with
R_n=Z[phi]/(phi^n−1) at odd n. At even n, exactly three words read as zero:
the empty support and the two alternating words.

[The dynamic refinement](golden_fibres/GOLDEN_FIBRE_DYNAMICS.md) proves that,
for every even n>=4, the ring value alone cannot carry named legal bit toggles.
More strongly, the current ring value together with its value after **one
fixed named toggle** reconstructs the entire word. A three-way tag only at
zero is the exact state repair; it adds two states, not a duplicated ring.

Randomizing the attempted toggle does not remove the loss: transitions back
to the zero fibre have probability0 from the empty word but1/2 from an
alternating word. The partition is not strongly lumpable. Period2 is explicitly
exceptional because its ring is a singleton and the observer reveals nothing.

Odd period here is an arithmetic fibre statement, not odd-order group descent:
R_3 already has4 elements. The Session90 tournament-automorphism hypothesis
cannot simply be transplanted to this quotient.

## 3. Hard guards turn Fibonacci counting into a much larger operator problem

Let H be a constraint graph, I(H) its number of independent sets, B the sum
of all free bit flips, and P the independent-set projector. The legal operator
is A=PBP. Its missing two-step return is

    D = PB²P − A² = PB(1−P)BP.

[The all-graph rank theorem](operators/GUARDED_FLIP_LEAKAGE.md) proves over Q/R

    rank D = I(H) − 2^b,

where b counts bipartite connected components, including isolates. The proof
identifies D with a signless exchange-graph Laplacian plus pinned states.
It also permits any nonzero real named flip weights through an exact product
gauge. This preserves rank, not necessarily an integer lattice or modular rank.

A two-state Fibonacci transfer counts independent sets of a path, but its
exact boundary-return Gram has rank F_(n+2)−2. A short counting formula therefore
does not imply a short linear interface for later operations. This is a precise
continuation of the even-graph Fourier work: the free Walsh transform remains
valid, but it does not erase the feasibility projector.

The legal operator is also the standard independent-set/PXP Hamiltonian in
Rydberg-blockade models, with primary sources cited. That is an exact operator
identification; it is not a QCD model or a physical eleven-square packing.

### The eleven-state pentagonal example

For H=C5 there are11 feasible words. Under phi=4 modulo11, legal insertions
form a15-arc directed subgraph of Paley11. The undirected state graph is a
wheel whose five rim edges are subdivided. The ambient unguarded ±phi^i
translations would instead give K11.

[Its exact leakage blocks](operators/LUCAS_FIVE_GOLDEN_OPERATOR.md) are

    0,       2I+Adj(C5),       3I+Adj(C5).

The smallest positive eigenvalue is phi^(−2). Their Smith forms give the full
cokernel Z ⊕ Z/5 ⊕ Z/100. Modulo5 the second block is nonzero nilpotent, with
successive ranks3,1,0. Thus an odd-prime defect can belong to a guarded Gram
even though Session90's dyadic Fourier coordinate change remains invertible
at every odd prime. The two operators must not be confused.

[An integral golden frame](zeckendorf/GOLDEN_FIVE_CYCLE_FRAME.md) identifies
the zero-sum lattice of the singleton block with two copies of Q², using a
unimodular basis. The index5 occurs only when splitting off the constant mode
from the full integer lattice. This makes the connection to the arithmetic
register exact and locates the characteristic-five failure.

At even cycle lengths, leakage rank happens to equal |R_n|, but the kernels
are different. The operator note gives the hostile and a repair: the golden
pushforward plus the two alternating-word Gram observations reconstructs
all amplitudes. Equal ranks alone do not supply that map.

## 4. Eleven-square arithmetic and seven-point geometry have different symmetries

The [packing and metric note](packing/GOLDEN_METRIC_AND_PALEY_TORI.md)
recovers the exact problem: eleven congruent unit squares, arbitrary rotations,
and disjoint interiors in a square container. It cites the current reported
optimality proof and its native-certificate trust model; this session does not
rebuild that formalization. The older supplied golden paper concerns a marked
Trump contact family, not all packings, and was not independently reread here.

The golden quotient has11 cells and a multiplier of order5. Its lifted matrix
Q preserves no nonzero symmetric bilinear form and hence no Euclidean metric;
Q⁵ may fix the quotient cells while its lift still expands. A marked arithmetic
torus is not automatically a metric symmetry of a physical packing.

There is a sharp positive contrast. Among odd Paley affine groups
p:((p−1)/2), precisely p=3 and7 admit faithful flat two-torus isometry actions.
At7 an explicit Eisenstein metric realizes the group7:3 of order21 and seven
pairwise equidistant torus points. Their triangular quotient has7 vertices,
21 edges and14 faces. The faces are exactly the cyclic triples of Paley7 and
split into two Fano planes, whose binary spans are two embedded [7,4,3] Hamming
codes with intersection {zero, all-one}.

This is the classical toroidal K7 construction in an explicit dictionary,
not a planar unit-distance K7. Shortest edge lifts already have a triangle
whose displacement sum is a nonzero lattice vector. The missing cycle
coordinate prevents a common planar lift. This is a literal geometric version
of the spanning-tree integration versus cycle-consistency distinction in our
trijections, rather than a shared numerical seven.

## 5. Geometric Euler repair and the LRC boundary

[Geometric padding](geometry/GEOMETRIC_EULER_REPAIR.md) extends Session90's
sparse graph construction. Distinct rigid trees with integer-coordinate
unit-distance realizations can repair a connected core to e=v−1 while
preserving its abstract automorphism group, Royle sign, faithful realizability
and optional minimum point separation. The added size is O(beta log beta),
matching the earlier lower bound within the disjoint odd-symmetry contract.

The containing square and density are not preserved. Nor are ambient Euclidean
symmetries: a core may retain abstract D8 while the displayed enlarged point
set has trivial Euclidean symmetry. Side contacts versus corner contacts of
unit squares give another explicit predicate mismatch.

[The LRC refinement](geometry/LRC_REFLECTION_AND_GOLDEN_GAP.md) proves that
mixed-parity positive-threshold good-time persistence is entirely doubled by
reflection, so even all cross-threshold signed traces remain blind. The
classical opposite-slope contact argument instead gives rational maximizing
clocks of denominator at most2M, where M is the largest integer velocity.
Its elementary quadratic-norm consequence is

    |maximum loneliness − (2−phi)| > 1/(12M²).

This excludes a literal finite golden optimum, while allowing limits and
properly normalized comparisons. It is not a new LRC lower bound or a closure
of LRC(14).

## Verification and open directions

[Independent audit](audit/INDEPENDENT_AUDIT.md) records separate automaton
minimization, exact modular and incidence calculations, geometric distances,
word fibres, contact clocks and counterexamples. All released checkers use
standard-library exact arithmetic and active checks under normal Python and
`python3 -O`. Generated receipts are local outputs, not required fixtures.

The corpus gains explicit observation quotients, a general guard-rank theorem,
an integral golden operator bridge, and geometry-preserving constructions.
Natural unlabelled correspondences, unrestricted closed tournament carriers,
LRC(14), planar unit-distance/chromatic problems and positive-integer Collatz
convergence remain open in the scopes discussed here. No physical packing
optimality is inferred from a finite-ring count or a shared golden eigenvalue.
