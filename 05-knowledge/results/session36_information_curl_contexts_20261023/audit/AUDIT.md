# Session36 independent mathematical audit

Verdict: PASS within the scopes below, after the explicitly recorded repairs. No Git operation was performed. This report does not assert literature novelty.

## 1. Proof audit

### Integral metric curl and finite theta recovery

Reviewed `session36-lattice/METRIC_CURL_RECOVERY.md`.

- Equality of restricted forms forces a difference gradient whose potential has at most two levels. The tournament signs then force precisely a uniform cut reversal. The exact fixed-label fibre consists of the cyclic rotations of the ordered SCCs, with cardinality equal to their number.
- The root-lattice isometry argument is valid, including the rank-zero/rank-one exceptions. A global minus sign cancels twice in an alternating pairing; it does not automatically identify the converse tournament.
- The pencil identity `det(tG-K)=1^T adj(tI-M)1` has the correct covolume factor. The real orthogonal spectrum and the integral Smith type are genuinely weaker than the simultaneous integral metric pairing in the stated examples.
- At genus n−1 the fixed Gram matrix I+J describes integral bases, since its determinant equals the lattice determinant. Consequently equality of the joint Gram/alternating coefficient family at that Gram matrix recovers the simultaneous isometry class. This is not a theorem about an undecorated scalar theta series.
- The universal two-root distribution for skew conference matrices is correct. The counts at pairing ±4 follow from directed four-cycles; the tight-root-frame second moment then gives pairing ±2 and zero. The displayed conference adjacency polynomial gives the claimed four-cycle count.
- The marked-vertex outdegree resolves the remaining cyclic cut. A lattice-only root marking cannot do so.

Required repair found and confirmed: “odd rank” for even tournament order was replaced by “odd dimension”; alternating ranks are always even. The score in the final marking statement is now explicitly outdegree.

### One-hole module and quotient contexts

Reviewed `session36-context/RESULT.md` through its current main module/context results.

- The complete path-cover profile is exactly the state preserved by every one-hole uniform module context. The cyclic contexts have the claimed triangular response, with top coefficient (m!)².
- At known order n, r−1 probes, r=ceil(n/2), determine the profile using parity and the universal leading coefficient. The exact rational connection rank is r, witnessed by the near-transitive reversed-arc family.
- The linear/affine minimality wording is correctly restricted. It is not a lower bound for distinguishing the finite set of realized profiles with arbitrary integer-valued probes.
- The quotient-context Stirling inverse has the correct direction and factorial factors. Nonempty-factor observations directly recover the full-support walk coefficients.
- The signed-forest expansion follows by consecutive-arc product expansion; “directed forests” must be read as signed orientations in the complete vertex set, not merely paths following the tournament's positive arcs.

The later C3-structure/walk-completeness extension is separately audited in `C3_WALK_ADDENDUM.md`.

### Rooted curl decoder and enumeration

Reviewed `session36-root/CURL_REALIZABILITY.md`, `CURL_CONNECTED_GRAPH_ENUMERATION.md`, and `UNLABELLED_CURL_NECKLACES.md`.

- The four sign constraints are necessary and sufficient. Pin propagation followed by the linearly ordered free SCCs gives exactly b+1≤n realizations. No hidden higher tournament compatibility condition remains.
- The all-triangles-local-alphabet hostile K=(-3,-3,-3) is valid, cocyclic, and globally inconsistent. The latest script explicitly checks this documented hostile.
- The labelled logarithm follows from free rotations of disjoint labelled SCC blocks. Its equality with connected labelled graph counts is formal and all-order; the EGF has radius zero.
- Hamiltonian-weighted logarithms are correct over the stated commutative Q-algebra.
- The unlabelled Burnside formula, inverse Dirichlet coefficients, primitive-necklace Möbius formula, and permutation stabilizer `r·product Aut(S_i)` are correct.
- The even-graph bridge uses the published automorphism-edge-reversal definition, not even vertex degree. Closure under disjoint union and the Euler-product argument correctly identify primitive strong necklaces with connected even graph counts. This primitive-necklace/connected-even equality and an explicit rank codec were already developed in the October11 notes `tournament_even_bijection_20261011c.md` §§2–4 and `tournament_atoms_abc_and_plane_20261011c.md`; Session36 adds the metric-curl interpretation, not that count identity.
- The weighted unlabelled formula requires the Adams operation s→ds. The inverse retains that operation. Omitting it first fails at n=6, with the repeated C3 witness.
- There is no S_n-equivariant bijection between labelled curl classes and connected labelled graphs for n≥3: the latter has the fixed complete graph, while an invariant alternating triangle tensor must vanish and no tournament tensor does. The n=1,2 boundary is sharp.

The decoder hostile and Q-algebra wording repairs have been confirmed in the current files.

## 2. Independent finite certificates

All scripts below are dependency-free except imports among the audit's own files; none imports producer or repository mathematics code. All substantive checks use explicit exceptions, not removable assertions.

1. `independent_audit.py` / `independent_receipt.json`:
   - Every labelled tournament of orders1 through6, total33,867.
   - Every alphabet-valid rooted curl input through order6, total1,052,742; this extends the producer's invalid-input census by an entire order.
   - Direct bitset sign-assignment feasibility, Tarjan SCCs, explicit component rotations, constant H, and Hamiltonian-weighted logarithms.
   - Direct connected labelled graph enumeration independently gives curl-class counts1,1,4,38,728,26704.
   - 5,625 arbitrary module contexts, 2,120 cyclic contexts, and exact rank witnesses through order8, using independent permutation/cut and maximal-run word enumeration.
2. `independent_theta.py` / `theta_receipt.json`:
   - Direct coordinate-vector star-root basis counts for n2..6:2,12,48,240,1440.
   - All128 specified Paley8 switches have the same genus-two root coefficient distribution.
   - The selected genus-three coefficient is36 on112 switches with H333 and84 on16 switches with H189.
3. `independent_necklaces.py` / `necklace_receipt.json`:
   - Complete isomorphism-orbit removal through n6, independent SCC-rotation quotient, and all76 tournament-isotype stabilizers.
   - Curl metric-isotype counts1,1,2,3,9,45.
   - Adams-weighted checks at powers0,1,2 through n6. The naive omission produces1130 instead of1133 for the first moment, and34561 instead of34597 for the second.
4. `independent_even_graphs.py` / `even_graph_receipt.json`:
   - Complete simple-graph isomorphism-orbit and automorphism-reversal-parity enumeration through n6.
   - Even graph type counts1,1,2,4,12,56; connected even type counts1,0,1,2,8,43, agreeing with primitive necklaces.
5. The C3 extension has its own exact certificate detailed in the addendum.

Normal and `python -O` output files are byte-identical for all audit scripts. The finite tests support conventions and implementations; the general conclusions depend on the proofs above, not extrapolation.

## 3. Source check and exclusions

Independently opened the primary article “Tournaments and Even Graphs are Equinumerous,” arXiv:2204.01947v2, and confirmed its Theorem1.1 and definition:
https://arxiv.org/html/2204.01947v2

The context addendum separately verifies the interval-inversion literature theorem. The underlying external rank88 lattice construction and incidental exploratory score-response files were not independently reverified in this audit. No claim about arbitrary metagraph wiring, polynomial-time Hamiltonian counting, universal genus-two theta equality, or literature priority is certified.

## 4. Cross-lane consequence

The metric curl class determines the entire path-cover profile, not only H: the path-colour polynomial is multiplicative across ordered SCCs, so cyclic rotation preserves it. Hence metric-curl equivalence implies equivalence under every one-hole uniform module context. This does not make the two equivalence relations identical.
