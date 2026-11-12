# From one five-state pentagon to an all-scale coherence grammar

**Status:** PROVED in the marked transitive-tail family specified below;
FINITE-EXACT controls in `RESULTS.json`; known polytope identification CITED.
No universal tournament classification or Chan–Pak coverage is claimed.

## 0. Inheritance and what changes

The [ready-return pentagon](../../ready_return_pentagon_20261103/PENTAGON.md)
already constructs the five source tails

    E — A — AB — BA — B — E

by insertions, one arc reversal and deletions, commuting with its 11-prefix.
It expressly leaves higher operational coherence open. The earlier
[Collatz pentagon note](../../collatz_pentagon_operator_coherence_20260930.md)
already distinguishes a cocycle/commutator identity from a physical
associator. We retain those distinctions.

Here the five-cycle is embedded in a complete all-size operation grammar.
The ambient graphs are **stellohedra**, a classical object, not a newly
invented polytope. The new application is their marked-tournament realization,
a self-contained quantitative coherence normal form, and precise failures
of pentagon-only or scalar-only continuation. This is a family-relative answer
to the old operational-coherence question, not an associator for every
possible tournament attachment.

## 1. Partial-permutation growth states

Fix a nonempty tournament D on a named vertex set B. Fix k distinct potential
new names 1,...,k. For each name a specify exactly which core vertices it
beats, its **signature** S_a subset B. Repeated signatures are allowed, but
the names remain distinguished.

For any injective word w=(a_1,...,a_m), including the empty word, let T(w)
consist of D and those m named vertices. Keep every prescribed core incidence;
orient the tail transitively, a_i→a_j iff i<j. This is a genuine tournament.

Allowed elementary edges are:

1. append an unused name to the tail;
2. delete the last name;
3. exchange adjacent tail names, reversing exactly their mutual arc.

No other arc changes. In particular this is a **graded marked** operation
metagraph. It is not the fixed-order unmarked metagraph, and its elementary
arc reversal is not asserted to be a single move in every restricted
Hamiltonian-spine gadget calculus.

**Theorem H1.** Its marked states and edges are exactly the graph G_k of
partial permutations with those three moves. Every vertex has degree k and

    v_k = Σ_{m=0}^k k!/(k−m)!,       v_0=1,       v_k=k v_{k−1}+1.

For k=0,...,7 the counts are 1,2,5,16,65,326,1957,13700.

**Proof.** The named vertex set recovers the used names and their transitive
tail order uniquely. The core incidences are fixed, so this map is injective.
For m>0 there are m−1 adjacent swaps, one deletion and k−m appends. For m=0
there are k appends. Counting words by length gives the sum; counting by the
first name gives the recurrence. ∎

This is a complete recursive description of this family at every size.
Repeated signatures may collapse different states after forgetting names;
no unlabelled class-count formula follows from v_k.

### The classical polytope identification

Let the star have hub 0 and leaves 1,...,k. A tube is a nonempty proper
connected induced vertex set. Compatible tubes are nested, or disjoint with
disconnected union. For w=(a_1,...,a_m), form the k tubes:

- `{0} ∪ ([k] minus {a_1,...,a_j})` for j=1,...,m;
- `{a}` for each unused leaf a.

They are a maximal tubing. Conversely the central tubes of any maximal
tubing form a chain; maximality forces its successive discarded leaves one
at a time, while the remaining leaves supply singleton tubes. This recovers
w. Exchanging one tube exchanges adjacent word positions or adds/removes
the last position. Thus G_k is the one-skeleton of the star
**graph-associahedron**, the k-stellohedron.

The underlying classical construction is due to Carr–Devadoss,
[Coxeter complexes and graph-associahedra](https://arxiv.org/abs/math/0407229);
the partial-permutation description and terminology are in Postnikov,
Reiner and Williams, [Faces of generalized permutohedra, §10.4](https://math.mit.edu/~apost/papers/genpermfaces.pdf).
The explicit tournament realization and the normal-form proof below do not
require a geometric realization of that polytope.

## 2. Two finite relation schemas control every scale

For any prefix p and distinct unused names a,b, take the **tail pentagon**

    p → pa → pab → pba → pb → p.

For any adjacent swap w↔s_i w and an unused name c, take the
**append/swap transport square**

    w → s_i w → (s_i w)c → wc → w.

These diagrams are actual named tournament diagrams. Every loop restores
all vertices present at its base state and every arc; no claim is made that
H is constant along an edge or that readiness defects stay fixed.

**Theorem H2 (all-scale coherence).** Backtrack cancellation, these pentagons
and these transport squares generate every closed walk in every G_k. In
particular the permutation braid hexagons and commuting-swap squares are
consequences of these two schemas. The same is true in the increasing union
with countably many available tail names, because every finite loop uses
only finitely many names.

**Proof.** Give each word w its canonical path P(w), which appends its names
in order from E. For an append or deletion edge w↔v, the two canonical paths
agree after a backtrack is removed. Suppose v is obtained from w by swapping
positions i,i+1, indexed from 0. Write w=pab s. Use a transport square for
each letter of suffix s to move the swap back to the state pab. There the
pentagon identifies the append path through pa with the append path through
pb followed by the swap. Hence P(w), the edge w→v and reverse P(v) fill using

    |s| transport squares + 1 pentagon = |w|−i−1 cells.

For an arbitrary loop concatenate these fillings for its edges. Adjacent
canonical paths cancel. This fills the loop. ∎

The observation is stronger than “there is a pentagon at every size”: its
**transport law** is part of the complete relation system. New states keep
appearing, but no new type of coherence relation is required in this marked
family. The first pentagon occurs at k=2; transport squares already become
necessary at k=3.

### A minimal relation spine, not merely a spanning list

**Theorem H2b.** Attach one cell for each distinct tail pentagon and each
transport square. This two-dimensional complex collapses to the canonical
append/delete tree. Its cell boundaries form an integral basis of the graph
cycle lattice, with no redundant relation. If P_k and S_k count its pentagons
and squares, then

    P_k=(v_k−k−1)/2,
    S_k=((k−3)v_k+k+3)/2,
    P_k+S_k=|E(G_k)|−v_k+1.

**Proof.** The append/delete edges form a spanning tree: every nonempty word
has the unique parent obtained by removing its last letter. Every remaining
edge is an adjacent swap. If the swap exchanges the final two positions, it
is the unique top swap of its tail pentagon. Otherwise remove the common
last letter; the swap is the unique top-length edge of a transport square,
whose other swap has smaller word length. These assignments are bijections
between non-tree edges and cells. Order cells by descending top word length.
Each top swap then belongs to only its own remaining cell, so it is a free
edge and the pair can be collapsed. All non-tree edges and all cells vanish,
leaving the tree. In cycle coordinates relative to that tree the boundary
matrix is triangular with diagonal ±1, proving the integral-basis claim.
There is one terminal swap for each word of length at least two, counted
with its two orientations, giving P_k. Subtract this from the number of
non-tree edges to obtain S_k. ∎

Thus the six pentagons plus three transport squares at k=3 are exactly nine
independent compatibility equations. The braid hexagon is not an additional
independent equation. The number of equations grows, but the two rule types
and their elimination order do not.

## 3. Exact and approximate integration

Give every directed edge an antisymmetric vector increment δ(u,v), in a
normed vector space. Define a potential from an anchor F(E) by

    F(w) = F(E) + sum of δ along P(w).

**Theorem H3 (local-to-global error invoice).** If the circulation of δ on
each tail pentagon and transport square has norm at most ε, then an adjacent
swap at positions i,i+1 of a word w satisfies

    ||δ(w,s_iw) − (F(s_iw)−F(w))|| ≤ (|w|−i−1) ε.

Append/delete edges satisfy exact equality. Every closed walk γ therefore
has

    ||sum_γ δ|| ≤ ε Σ_{swap edges of γ} (word_length−swap_position−1).

In particular zero local defects imply that δ is globally the exact
difference of a potential, unique after the anchor is set.

**Proof.** At an unsuffixed swap pab→pba its discrepancy is exactly the
pentagon circulation. Adding one common final letter subtracts the transport
square circulation from that discrepancy. Induction over the suffix proves
the displayed bound; sum discrepancies around γ, whose potential differences
telescope. ∎

For exact group-valued transports, replace additive circulation by ordered
products and reversal by inverse. H2 likewise says that identities on all
displayed faces imply path independence. The approximate vector bound does
**not** transfer unchanged to matrices: operator norms and conjugation losses
would have to be retained. Nor are tournament insertion maps automatically
invertible linear maps on a chosen response vector; that is a separate
hypothesis whenever the transport language is applied.

**A useful negative dimension statement.** Flat scalar increments on G_k
have v_k−1 degrees of freedom; adding the potential anchor gives v_k. Every
assignment of v_k vertex values supplies an exact form, and the connected
graph makes the kernel of the difference map precisely the constants.
Thus the pentagon case has five potential parameters, but the next case has
sixteen. Coherence alone cannot turn the whole recursive family into a
five-parameter counting model. A further analytic closure law would be needed.

## 4. Two small hostiles, with different mechanisms

### H4a. All pentagons can pass while higher coherence fails

On G_3 put δ(123,213)=1, put −1 on the reverse edge and zero elsewhere.
There are six tail pentagons, and none uses this first-position swap at
length three; every pentagon circulation is zero. But

    12 → 21 → 213 → 123 → 12

has circulation −1, while the full-permutation braid hexagon

    123 → 213 → 231 → 321 → 312 → 132 → 123

has circulation 1. The missing constraint is a transport square, not another
copy of the same isolated five-state model. This is the smallest k permitting
the failure; G_2 is just its pentagon.

### H4b. A physical H quotient is not a continuation state

Use core vertices x=0,y=1 with y→x. Tail A beats x and loses to y;
B and C beat y and lose to x. Then

    H(T(AB)) = H(T(BA)) = 5,
    H(T(ABC)) = 11,      H(T(BAC)) = 13.

All four are literal tournaments; the append operation uses the same named
C and the same core signature. The old graphs have four vertices, the new
ones five. Direct permutation enumeration in `verify.py` checks these counts.
Consequently identifying states merely because H agrees does not give a
well-defined append action. The failure exists within one fixed-signature
stellohedral family, not only between unrelated constructions.

This is a control illustrating an already known general quotient hazard,
not a claim to have first discovered that H loses continuation data. It also
explains why a full-state returning loop is compatible with rigidity of a
fixed-deletion binary fiber: the loop changes tail size or arcs, leaves that
fixed fiber, and returns. It supplies no nontrivial motion inside a rigid
fiber.

## 5. The bridge to stress trajectories and other histories

| Source → target | Preserved predicate | Lost information / required sidecar |
|---|---|---|
| Ready-return five-state cycle → G_2 | actual insert/delete/flip diagram | larger G_k need not retain the ready-family joint defect |
| Marked fixed-signature tournaments → stellohedron | all three operation adjacencies | unmarked isomorphism quotient may merge states; preserve names/signatures |
| Local pentagon and square defects → global transport | exact path independence; quantitative additive error bound | a base anchor and the square transport law are indispensable |
| Stress-potential integration → discrete edge-potential integration | closed data plus global topology control existence of a potential | stress PDE and analytic amplitude law are not graph flatness; stress domains may have genuine periods |
| Scalar path counts → a putative machine state | current output only | append fails H4b; retain enough continuation/interface data |

The mechanism recurring across the user's five themes is not that all sets
of five are secretly the same. It is that a small local compatibility diagram
may govern arbitrarily many configurations **when a transport rule proves
that it does**. Here that rule is exact and uniform. In the stress system the
finite-dimensional Hermitian amplitude law provides additional analytic
closure; in this tournament operation family, flatness leaves v_k free
potential values instead. This contrast is a productive obstruction rather
than a reason to abandon the analogy.

## 6. Next constructive questions, explicitly OPEN

1. Find a response tensor for a useful Chan–Pak source family on which the
   stellohedral edges have closed, physically realizable transport, retaining
   the absolute count anchor. H2 would then audit its path independence using
   only pentagons and squares, but existence of such a compact tensor is not
   proved here.
2. Classify signature assignments for which the H quotient is a right
   congruence under every allowed append. H4b excludes arbitrary signatures;
   all-equal source/sink signatures provide elementary degenerate controls.
3. Compare response-tensor dimension with the stress-style finite-jet closure:
   identify extra recurrences beyond flatness which reduce v_k without losing
   physical target realization. The local-to-global error invoice can test
   proposed approximate recurrences against a precise absolute defect.

No target interval, universal odd-count coverage, or unrestricted prime-core
classification follows from these coherence results.
