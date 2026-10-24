# Six-slot moves generate every integer charge relation, but not every native route

## 1. The exact carrier

Use the authenticated 28-source supply with targets U=1965 and V=2234. Each
source contributes one occupied path slot. Write b_s=1 for its V path and b_s=0
for its U path. Every common owner clock is retained. The independent pinned
supply certificate proves these are the complete domains once the original
U/V margins are fixed: source 1923's extra path goes to 1559 and cannot occur
when that target's required nonnegative count is zero.

All 56 admissible options have distinct normalized black resources. Source
matrices, owner frames, complement, converse and actual lifts remain in the
supply certificate. No raw hidden lift is extra capacity. These are selected
packet assignments, not a change to the deterministic full-atlas complement.

Let g_s be the original paired-charge label in

    G = (Z/2)^4 plus Z/128.

The exact guarded assignment fiber is specified by

    sum_s b_s = k,       sum_s b_s g_s = q.

The common target labels agree in G. Thus the source-row current changes vanish;
the displayed total group sum tests all target-row dyadic defects. This uses the
separately audited original marked harmonic-bank theorem, not an assumed new
low-modulus pattern. Individual ordinary per-port characters are a different
predicate; the selected X-mod-8 phase remains unchanged because X_U=X_V.

## 2. Four-slot obstruction and six-slot lattice saturation

Enumerate every pair of disjoint two-subsets with equal G sum. There are exactly
27 such four-slot trades. Their signed incidence vectors have rational rank 23
in Z^28. In particular every one preserves the exact integer quantity

    b_1638 + b_1974.

This is a linear obstruction to routing endpoints that change that sum. It
persists under any legal context for these four-slot trades, regardless of the
number of temporary states or which other slots are occupied.

There are exactly 1917 disjoint three-versus-three equal-charge trades. Adjoin
these six-slot vectors to the four-slot vectors. Their integer span is EXACTLY

    K = { z in Z^28 : sum z_s=0 and sum z_s g_s=0 in G }.

Here is a short independently checkable lattice certificate. Delete the last
coordinate of the zero-sum lattice, identifying it with Z^27. The differences
g_s-g_last generate all 2048 elements of G, verified by complete finite group
closure. Therefore K has index 2048 in Z^27. Every enumerated move lies in K.
Among the first 48 actual move vectors, saved integer linear combinations give
27 basis columns with determinant of absolute value 2048. Direct multiplication
checks every combination, and exact Bareiss elimination checks the determinant.
The generated sublattice therefore has the same index as K, proving equality.
No integer-normal-form algorithm is trusted by the verifier.

Consequently the four-slot rational obstruction disappears at six slots, and
there is no additional additive integer or modular obstruction beyond the
original total degree and G charge. This is NOT yet binary-capacity connectivity.
The saved integer coefficients may be negative or require repeated uses that
cannot be realized inside a current unit-capacity assignment.

## 3. A complete native binary fiber still needs eight simultaneous changes

In the full 28-source supply, fix k=4 and group charge q=(7,15), where 7 encodes
the four parity bits and 15 is the cyclic coordinate modulo 128. Exactly two
four-subsets have this charge:

    A = {176,1425,1549,2087},
    B = {1394,1923,2837,3369}.

These are disjoint. Put the sources of A, respectively B, at target V and all
other sources at U. All 28 slots, not just these eight, are allowed to vary
subject to their original pins and guards. Nevertheless these are the only
two states in the complete fiber. Therefore no route using at most six changed
slots per event connects them, even though all six-slot moves together generate
the entire integer charge kernel. Their direct eight-slot switch is legal, so
the exact bottleneck is eight.

The 20 other occupied source slots are already included in this claim. This
is not an assertion about adding new sources, changing their initial margins,
releasing pins, or enlarging the original metagraph. The actual slot paths and
all normalized resources are supplied by the independent complete-domain audit.

## 4. Full low-cardinality census

The verifier independently enumerates all subsets at cardinalities 2, 3 and 4,
groups them by actual G charge, and connects two states precisely when their
Hamming distance is at most six. This tests every admissible two-, four- or
six-slot replacement in the declared supply, not merely a selected move list.

- Cardinality 2: 378 assignments, 351 nonempty fibers; all connected.
- Cardinality 3: 3276 assignments, 1646 nonempty fibers; all connected.
- Cardinality 4: 20475 assignments, all 2048 charge fibers nonempty.
  Of these, 1851 are connected, 175 have two components, 21 have three and one
  has four components under the six-slot event bound.

Thus 197 different cardinality-four fibers exhibit this capacity-sensitive
failure. The two-state example above is a particularly transparent certificate.
The census is finite-exact; it does not classify higher cardinalities.

## 5. A general contextual decomposition rule

The correct object for legal routing is not merely the integer move lattice.
Let A be a finite set of actual admissible assignments, retaining all shared
owner coordinates, normalized resource capacities and guards. For each legal
native event x -> y, form e_x-e_y in the free vector space with basis A. The
span of these differences has codimension equal to the number of connected
components of the undirected legal-state graph.

Proof: every difference has coefficient sum zero separately on each component.
Conversely, a path from x to a chosen component root expresses e_x-e_root as a
sum of event differences. These differences span all vectors with zero sum in
each component. Thus the component sums give exactly the quotient. QED.

For binary source choices, a formal trade A <-> B is usable only in a context C
for which both A union C and B union C are actual legal states and every resource
is distinct. When the supply is jointly compatible and the trade sets are
disjoint, C must avoid A union B. If shared variables or target-side attachments
remain coupled, they are additional context predicates, not discarded metadata.

This is why ordinary polynomial or integer-lattice manipulation can overstate
native reachability. Multiplication of an abstract binomial by a repeated port
may be an algebraically valid operation but an inadmissible context. Simply
quotienting by square-zero variables is not a substitute: one side of a binomial
can vanish while the other remains, producing a relation to zero rather than a
legal two-state move. The actual context condition must remain explicit.

The state-graph component invariants can be nonlinear functions of the original
occupation bits. The native two-state hostile has such a separation even though
the global six-slot move lattice has no extra additive invariant. This connects
the electrical viewpoint to decomposition without conflating two graphs: the
original metagraph supplies the harmonic guards; the legal assignment graph
supplies the remaining reachability components.

A physical analogy is hard-core occupancy versus unrestricted multiplicity.
It is only an analogy: these commutative resource choices do not automatically
carry fermionic exchange signs or a quantum field theory.

## 6. Reproduction and boundaries

    python3 -B check.py --repo REPO --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out normal.json
    python3 -O -B check.py --repo REPO --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out optimized.json

Both modes pass 10184 checks and produce identical receipts. The witness
contains only exact integers. The optional generator is not needed to verify
it. No universal six-slot Markov basis, complete atlas rewiring, or all-order
native supply theorem is asserted.
