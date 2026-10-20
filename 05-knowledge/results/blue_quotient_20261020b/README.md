# Blue-quotient connectivity and Petersen reductions beyond flowers

20 October 2026, continuation b. Baseline: d9674b6a00.

**Current continuation:** [all-order descent and cyclic-five reductions](../black_descent_20261020c/README.md)
proves the global black-descent obligation left open here, supplies an
SCC interchange and gives double-star, Loupekine and infinite
cyclic-five nonflower Petersen reductions. Universal connectivity and
unavoidable snark reductions remain open.

**FINITE-EXACT:** the black-active core is connected at tournament order
nine; the full metagraph has five bridges, all leaf edges.
**PROVED:** every strong black-active tournament has a black closed
Hamiltonian path, and the rooted Kempe reductions below preserve full
labelled colouring support.
**OPEN:** all-order black-core connectivity and an unavoidable admissible
reduction for arbitrary snarks.

## 1. The remaining blue quotient, resolved at order nine

The [complete census and certificates](census.md) cover all 191,536
ordinary tournament classes, or 97,144 converse-merged nodes. An explicit
black spanning tree connects all 97,138 black-active nodes. Each of its
97,137 edges carries a Hamiltonian word and a partner-isomorphism witness.
A second canonicalization algorithm independently checks the complete
class universe, every tree edge, and all 65,536 blue tilings.

After black contraction, the blue quotient is exactly:

- one black-active hub with 32,760 blue loops;
- five pure-blue leaves, each attached by one edge;
- one pure-blue trivalent node, attached to the hub by three parallel edges.

There are therefore exactly **five full-metagraph bridges**, all incident
to leaves. The blue subgraph alone has **19 additional non-leaf bridges**;
all disappear when their endpoints are connected through black paths.
This is a concrete warning against declaring a bridge from a partial
colour subgraph.

The complete black edge set was unnecessary: complete vertex coverage
and an authenticated spanning tree suffice for connectivity. This is a
certificate of the finite assertion, not a sample-based extrapolation.

## 2. Structural progress toward an all-order proof

The [closed-path theorem](structure.md) proves that a strong tournament
has every closed Hamiltonian path blue exactly when it is C3 or the
regular five-vertex tournament R5. The odd-order proof reuses the
[audited prefix rerouting argument](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md)
and checks that its witnesses also close. The even-order proof constructs
a new cycle through the two orbits of a rotation and cuts it at
equal-score endpoints.

Two further proved mechanisms retain actual Hamiltonian words:

- Along a fixed Hamiltonian cycle, the set of blue cuts is an exact
  coset of a rotation subgroup. If any cut is black, at least two are.
- Starting from an open black path, two black path reversals realize a
  controlled two-arc change: its endpoint shortcut and one path edge.

For a path from a to b, the cyclic-triangle change is

    c3(T')-c3(T) = s(a)-s(b)-1.

Strict descent is proved for strong non-self-converse tournaments, regular
black-active tournaments, and tournaments with an all-blue-cut cycle.
The general black closed-path descent claim passed all classes through
order seven in this packet; it is now proved at every order in
[continuation c](../black_descent_20261020c/cycles.md).

The surviving difficulty is precise. A fixed cycle can have no descending
black cut: an explicit order-nine example has six black cuts, all with
s(b)=s(a)-1. Changing the cycle restores descent in that example.
Moreover, even a universal strong descent lemma would leave multiple
nonstrong SCC-word basins to connect. The local two-arc move can return
to the same merged node. These failures are retained as tests for a
future normalization proof.

## 3. Actual Petersen reductions for a nonflower infinite family

The [snark note](snarks.md) supplies two elementary rooted minor lemmas.
An isochromatic four-pole has two disjoint Kempe paths giving its
equal-pair wire minor. A heterochromatic four-pole has two disjoint
paths and a connector giving its two-adjacent-junction minor.
Both minors preserve the full labelled boundary support.

Apply the first lemma to the eight-vertex four-pole obtained by deleting
two adjacent Petersen vertices. Its extension vector is exactly twice
the two-wire indicator. Inserting it at two independent edges therefore
multiplies every retained boundary extension count by two. Reversing
the insertion gives an actual minor back to the preceding graph.

This yields the classical iterated Petersen dot-product family

    |V(G_r)|=10+8r,       G_r -> G_(r-1) -> ... -> Petersen.

Every step preserves the colouring obstruction. The exact tests include
all eight specified members from 18 through 74 vertices, with independent
perfect-matching colour counts and exhaustive cyclic-cut checks. These
are a deterministic family, not a census of all snarks at those orders.
The relation to the established four-pole factorization literature is
explicitly credited in the note.

The completed graph class is an essential side condition. Applying the
same wire replacement inside Petersen itself produces the looped
dumbbell, preserving uncolourability but leaving the class of simple
bridgeless snarks. The successful family avoids this failure because
each reduction restores its already specified predecessor exactly.

## 4. The cubic-vertex repair is realizable, and now fully measurable

The repaired seven-port kernel J consists of one all-distinct cubic
junction and two equal-colour wires. Two Petersen-pole insertions turn
it into a **connected 17-vertex seven-pole Q** with

    f_Q=4J,     support size 54,     total extension count 216.

Two authenticated rooted minor reductions return Q to J. This makes
the parity repair a concrete connected object.

The [seven-port measurement theorem](junctions.md) also finds exactly
what the 105 labelled J closures miss. On the 91 colour-permutation
orbits of parity-compatible words, their rank is 84. The seven missing
vectors w_i distinguish singleton port i in patterns (5,1,1) and
(3,3,1). For an actual pole P they have the Kempe interpretation

    6 <f_P,w_i> = -N_i,

where N_i counts colourings with a unique boundary colour at i and three
mixed-endpoint paths in the other two colours.

Seven additional closures K_i, each a five-internal-vertex path tree
with middle port i, satisfy

    <K_i,w_j> = -4 delta_ij.

Thus the **105+7=112 physical closure counts determine all 91 extension
coordinates**. The exact rank uses integer null vectors and a nonzero
minor modulo 101; independent tree-colouring and cut-XOR calculations
agree on all 2,187 words.

The five-vertex tree cannot replace J: its support has size 96.
More strongly, deleting two disjoint internal tree edges gives a
54-word J support **disjoint** from its original support. The connected
cyclic gadget succeeds where this naive tree simplification fails.

## 5. The Pythagorean connection retained by these reductions

For nonzero boundary vectors f,g, the vectors

    U=f tensor g,       V=g tensor f

have equal norms. Hence U+V and U-V are perpendicular, with squared lengths

    2(||f||^2||g||^2 + <f,g>^2),
    2(||f||^2||g||^2 - <f,g>^2).

Their nondegenerate isosceles case is exactly disjoint support:
<f,g>=0. General support-preserving reductions retain that endpoint.
The special Petersen insertion does more: f changes to 2f, so every
length doubles and the normalized triangle is unchanged. The connected
Q construction has scale factor four.

This is a precise correspondence at a marked interface. It supplies
neither a general graph-isomorphism dictionary nor a new proof of
Kuratowski's theorem or the general Petersen-minor statement.
The earlier primitive sequence through (3,4,5) remains valid but approaches
a flat triangle; it must not be confused with this isosceles limit.

## 6. Connection to the newly inherited tournament work

The [matrix-valued distance constraints](../merged_metagraph_distance_moment_constraints_20261020/PROOF.md)
add Walsh Gram positivity and exact fibre-compression loss. Its explicit
impossible rewiring preserves the entire complement graph, so these
constraints cannot by themselves settle black connectivity.

There is nevertheless a useful common mechanism. Both the tournament
quotient and the pole closure map lose coordinates. For seven-poles we
have now identified the loss exactly and supplied seven physical probes
that recover it. For tournament connectivity, a cut indicator h on
merged fibres must satisfy

    h(x)=h(complement x) whenever x is not reflection-fixed.

Equivalently its complement-antisymmetric part is supported on the blue
subspace. Classifying these compatible **Boolean fibre indicators**
remains the missing all-order step; positive Gram matrices alone do not
classify them.

The subsequently received
[endpoint-transport theorem](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md)
prescribes the entire coloured triangle-current row and folded-score
capacities from each tournament's Hamiltonian-path endpoint inventory.
It strengthens the necessary conditions on a proposed black-component
cut. Its inventory still omits closedness and destination identity
within a score shell, so it does not settle our descent or connectivity
obligation. Its elementary Hamiltonicity proof also provides a direct
audited input to our closed-path theorem.

The incoming
[weighted marked-branch correction](../weighted_branch_carrier_audit_20261020/WEIGHTED_BRANCH_AUDIT.md)
was checked for relevance. It concerns decorated tournament switching
classes, not the ordinary/converse-merged classes used here; none of
our counts or proofs depends on its corrected fixed-point formulas.

## 7. Reproduction and next precise obligations

- [Census](census.md): full rebuild and independent replay harness,
  compressed witnesses and source hashes.
- [Structure](structure.md): all small-class and literal hostile controls,
  including the even two-orbit construction.
- [Snarks](snarks.md): exact supports, rooted minor models, colourable
  controls, 18–74-vertex snarks and tree failures.
- [Junctions](junctions.md): all seven-port kernels, complete measurements
  and independent tree checks.
- [Audit ledger](AUDIT.md): independent checks, inherited corrections and
  exact limits of the conclusions.

At this checkpoint the next targets were cycle-changing descent,
boundary-compatible SCC normalization, and rooted replacements inside
cyclically five-connected cores. See
[continuation c](../black_descent_20261020c/README.md) for the proved
descent, a uniform SCC move, and cyclic-five reductions; general
normalization and unavoidable snark reductions remain open.
