# Binary endpoint intervals: implication ideals and neutral packet contraction

Status: **PROVED, abstract all-size finite-record theorems** with **FINITE-EXACT** controls. No native electrical-current realization, tournament supply theorem, novelty claim, or prize closure is asserted. The graph and Boolean arguments are elementary; inherited sources are archived in `inputs/` and pinned in `SOURCEPINS.json`.

## 1. Inheritance, carrier, and loss ledger

Session49 `nonreturn-repairs/THEOREM.md`, Sections 1–4, already proves the general endpoint-interval DAG, additive neutral partitions, the conditional Davenport bound, and the synchronized-block antimatroid criterion. Its Section 5 shows why an undirected connected fibre does not give a strict route. Session49 `hole-theory/THEOREM.md`, Sections 6–7, requires full witnesses, simultaneous slot/clock costs, and retained nonreturn history. These are dependencies, not new claims here.

The new specialization identifies the exact endpoint interval of arbitrary Boolean unary/binary relations and locates the missing obstruction in a proposed combination of small synchronization blocks with finite-group zero-sum theory: **neutral packet contraction may create directed cycles even when the original implication graph is acyclic**.

Each of finitely many labelled slots has two distinct full records. Keep source, target, actual path, lift, normalized resource, pins, population columns, charges, and actual changed-clock cost. For no resource reacquisition we assume all alternative per-slot reservations are globally disjoint, as in Session49's binary disjoint-reservoir contract. For per-slot record nonreturn, binary distinct full records suffice. A fixed-key lift model does not authenticate resource replacement. Here all models and charges are ABSTRACT unless separately authenticated by a native carrier.

Fix accepted endpoints x,y. Coordinates where they agree remain fixed in a strict binary nonreturn route. On each differing slot i use a switch bit z_i, equal to 0 at x_i and 1 at y_i. Its additive defect is

    d_i = guard_i(y_i) - guard_i(x_i) in A,

where A is an abelian group retaining every required charge and exact integer margin, for example Z^r × Q. Let c_i>0 be its literal changed-clock count. Acceptance is **exactly** the declared unary/binary relations and zero total additive defect; extra unencoded predicates are outside the theorem. In particular a finite quotient cannot replace an integer margin without a proved equivalence. No additional microscopic transition restriction is inferred from accepted arrivals.

The map to a switch set forgets no record when its two endpoint records and all fixed data remain attached. The subsequent SCC quotient forgets internal slot counts, clocks, and separate defects unless their exact sums and the member lists are retained. We retain all of them. Historical retired resources must be checked against the continuing contract before accepting any acquisition; choosing new endpoints never resets history.

## 2. Exact 2-CNF endpoint-interval classification

**Theorem 1 (implication ideals).** Restrict arbitrary Boolean unary and binary relation tables to the endpoint interval and apply the switch gauge above. The relation-only accepted states are precisely the order ideals of a finite poset of strongly connected implication blocks. Each block switches all its physical member slots together. The full accepted states are exactly the ideals whose aggregate defect is zero.

**Proof.** Every local relation accepts the two endpoint restrictions. A unary relation on a differing slot therefore accepts both switch values and imposes no restriction. A relation with an endpoint-fixed coordinate also accepts both values of its remaining differing coordinate. On two distinct differing coordinates its table contains (0,0) and (1,1). There are exactly four possibilities: all four pairs; forbid (1,0); forbid (0,1); or forbid both mixed pairs. They impose respectively no condition, the implication z_i⇒z_j, the implication z_j⇒z_i, or both implications. Repeated-coordinate binary scopes reduce to unary scopes.

Equivalently, express every forbidden Boolean pair by its two-literal clause. After substitution every surviving non-tautological clause has one positive and one negative literal: two positives fail at all-zero and two negatives fail at all-one. Thus the classification covers all 2-CNF formulas as well as all Boolean binary relation tables.

Put an edge i→j whenever selecting i requires selecting j. Along a directed path, selecting its first vertex forces all later vertices. Vertices in one SCC consequently have equal switch bits. Conversely any union of SCCs closed under outgoing reachability satisfies all implications. Order the SCCs by C≤D when selecting D forces selecting C. This is a partial order, and those unions are exactly its downward-closed sets, or ideals. If a block C has defect d(C)=∑_{i∈C}d_i, then the assembly at an ideal I has defect ∑_{C∈I}d(C). It is fully accepted exactly when this vanishes. QED.

Both the empty ideal and the full ideal P are accepted; in particular ∑_P d(C)=0. The relation-only family is a distributive lattice under union/intersection. Its zero-defect subfamily need not be closed under union, intersection, or deletion of a single block. Do not apply the relation-only lattice theorem to the fully guarded family without checking this distinction.

## 3. Indivisibility, exact costs, and the correct relative interval

Give each SCC block C its exact weights

    s(C) = number of physical differing slots in C,
    w(C) = sum of their actual changed-clock counts.

Let Z(P,d) be the zero-defect ideals. A strict route is exactly a chain from ∅ to P in Z(P,d). An event I⊊J has exact cost

    (s(J\I), w(J\I)).

All such routes have total physical work (s(P),w(P)); their event bottlenecks may differ.

**Theorem 2 (relative neutral ideals).** For accepted ideals I⊊J, the intermediate relation-valid states correspond bijectively to the ideals K of the induced poset P[J\I], via K↦I∪K. Intermediate accepted states correspond to those K with ∑_K d=0. Hence I→J is indivisible into smaller accepted monotone events if and only if J\I has no nonempty proper zero-defect ideal in its induced poset. These are exactly the cover relations of Z(P,d).

**Proof.** If I∪K is ideal, every predecessor in J\I of an element of K belongs to K. Conversely, if K is ideal in J\I and v∈K, every predecessor of v lies in J because J is ideal; one outside I therefore lies in K. Predecessors of I are already in I. Since I has defect zero, I∪K is accepted precisely when K does. QED.

Thus “minimal nonempty zero-defect ideals” describes first events from ∅; general events require the **relative** induced interval. A global neutral subset not closed under predecessors is not a legal event from ∅.

Every route can be refined until all its events are covers. Positive additive weights mean this never increases either bottleneck coordinate. Therefore the exact attainable joint event-budget region is the upward closure of the bottleneck pairs of cover chains. A finite DAG dynamic program computes the Pareto frontier; separately minimizing slots and clocks need not give an attainable pair. This is an exact finite algorithm, not a polynomial-size claim.

## 4. Exact neutral acyclic packetization and a constructive repair

A packetization is a partition Π of the SCC blocks into nonempty parts. Contract each part, retaining an edge A→B for every implication from a block in A to a block in B, with A≠B. This **packet dependency digraph** can have cycles even though the original SCC graph is a DAG.

**Theorem 3 (packet criterion).** A strict accepted route is equivalent to a packetization in which every packet has zero defect and the packet dependency digraph is acyclic, together with an ordering that executes required packets before packets depending on them. Its exact event costs are the sums of s,w on the packets. Conversely every strict accepted route gives such an acyclic neutral packetization.

**Proof.** In a route, take the blocks first switched in each event as its packet. Subtracting two accepted guards shows that packet has zero defect. If A depends on B, B is either the same packet or must occur earlier. A cycle among distinct packets would demand a strict cyclic ordering, impossible. Conversely an ordering of an acyclic dependency graph that executes requirements first makes every prefix a relation-valid ideal. Each prefix has zero defect because all packet defects vanish. No physical slot switches twice. QED.

**Corollary 3a (cycle-merging compiler).** Start from any partition into zero-defect packets, without requiring its dependency digraph to be acyclic. Compute that digraph's SCCs and merge the packets inside each SCC. Every merged packet remains neutral and the new dependency graph is acyclic. Executing the merged packets in requirement-first order gives an accepted strict route. Its exact event costs are the sums over each packet-SCC, not the old per-packet costs.

Every acyclic packetization is already a fixed point of this compiler. Hence minimizing the largest merged costs over all neutral packetizations is another exact formulation of the route bottleneck problem. The compiler does not promise that an arbitrary zero-sum extraction has small merged SCCs.

**Scoped sufficient bound.** Suppose every initial SCC block is exact-margin neutral, each has at most s physical slots and clock cost at most w, and its remaining charge is in a finite abelian Q. Partition its total-zero charge sequence into neutral packets of at most D(Q) blocks. (D(Q) is the inherited zero-sum subsequence constant.) If each SCC of this chosen packet dependency graph contains at most L packets, Corollary 3a gives simultaneous event bounds

    slots ≤ L D(Q) s,       clocks ≤ L D(Q) w.

This is conditional on the actual packet-SCC bound L. Without it, bounded original SCC size and finite D(Q) do not suffice. All integer margins must already be neutral in each initial block, or an authenticated margin-neutral packetization must be supplied before applying finite-group extraction. No claim about arbitrary overlapping population columns is made.

Two simpler sufficient cases are worth distinguishing:

- If every original SCC block is neutral, execute the original blocks in requirement-first order, with exact bottlenecks max s(C), max w(C).
- If the condensation poset is an antichain, every selected block subset is relation-valid. Any neutral partition is acyclic and the ordinary Davenport bound survives (subject to exact-margin neutrality).

A full Boolean block cube is equivalent to an antichain condensation: if C<D, selecting D alone is forbidden. More weakly, it suffices to authenticate legality of every prefix of the particular neutral partition selected. “Downward closed” alone is not closure under arbitrary zero-sum subsets.

**Accessibility boundary.** The zero-defect ideal family Z(P,d), containing ∅ and P, is accessible by one-SCC deletion if and only if every SCC defect is zero. Indeed accessibility repeatedly deletes one SCC from P down to ∅; subtracting successive zero sums forces each deleted defect to be zero. The converse is immediate from ideal accessibility. Thus the inherited one-block antimatroid route applies directly to the original SCC coordinates exactly in the all-neutral-block case. Larger authenticated neutral blocks can still give a useful antimatroid after regrouping.

## 5. Hostile chains with fixed group, no zero increments, and exact populations

**Theorem 4 (unbounded precedence atom).** For every n≥2 there is an abstract binary-relation instance whose relation graph is a path, whose switch implication SCCs have size one, and whose charge group is C2×C2, with every single-slot endpoint charge increment nonzero, but with exactly two fully accepted states. They differ on all n slots, so every route between them needs one event of cost (n,∑c_i). For every even n the example simultaneously has one exact fixed population of n/2 ones.

**Construction and proof.** Let independent nonzero group elements be a,b. Use switches z_1,...,z_n and implications z_{i+1}⇒z_i. Relation-valid switch sets are exactly the n+1 prefixes. Define prefix charges

    t_0=t_n=0;
    t_k=a for odd k, b for even k, 1≤k<n,
    d_k=t_k−t_{k−1}.

Every d_k is nonzero: the endpoints jump between 0 and a nonzero element, and every internal jump is between distinct a,b. The charge of prefix k is t_k, which vanishes only at k=0,n. Thus only empty/full prefixes are accepted. The undirected relation graph is a path, explicitly admitting bags {i,i+1}; it has treewidth one. All implications go one way, so SCCs have size one. Both all-zero and all-one switches satisfy every local relation; no higher-arity relational predicate is used.

For even n, set original records' Boolean labels x=(0,1,0,1,...) and y=1−x. Preserve exactly n/2 ones. The source-to-target population increments are +1,−1,+1,−1,..., hence the endpoints share that exact integer margin. The switch prefix of length k has zero population defect exactly when k is even; this cannot add any accepted state to the already two-point charge fibre. Assign each original slot an additive charge function with difference d_k across its two endpoint records. For instance q_k(x_k)=0 and q_k(y_k)=d_k. The switch implications pull back to ordinary two-literal clauses in the original labels, with the same path scopes. Choose two distinct globally disjoint normalized resources per physical slot. There is no reuse, phase-as-capacity, or hidden helper reset. The direct event changes every slot once. QED.

Here D(C2×C2)=3: among any three terms, a zero term works, a repeated nonzero term gives a zero pair, or the three distinct nonzero terms sum to zero; two distinct nonzero terms show sharpness. Nevertheless n is arbitrary. Finite zero-sum theory finds neutral subsets, not necessarily legal ideals.

The smallest chain in this family exceeding D(Q) has n=4 and charges

    (a, a+b, a+b, a).

The neutral packets {1,4} and {2,3} each have size two, but their contraction has edges in both directions: original implications 2→1 and 4→3 go between the packets in opposite directions. The cycle-merging compiler merges all four slots. This is an explicit failure mechanism, not a contradiction of Session49's additive-only theorem, which excluded unverified relational prefixes.

An even simpler, zero-increment-allowing chain uses (a,0,...,0,−a) in any nontrivial finite group. It already proves the unboundedness. The construction above rules out blaming the phenomenon on zero-charge coordinates.

### 5.1 Equal-population paired lift (optional, all logical lengths)

For each logical switch z_i introduce two physical slots u_i=z_i and v_i=1−z_i, requiring u_i≠v_i by a Boolean binary relation. Put the implication chain on the u_i. Original labels are (u_i,v_i)=(0,1) at source and (1,0) at target. Every relation-valid state has exactly one occupied label per pair, so both one total population and all individual pair populations are preserved. The primal relation graph is a path on the u_i with one leaf v_i attached to each, again a tree. In switch coordinates each pair is one SCC of size two; the SCC condensation is the same chain. Only the two endpoints are charge-neutral, and the event costs 2n physical slots with their full clock sum.

Even each physical slot's charge increment can be nonzero: split each nonzero d_i=e_i+f_i using the other two nonzero elements of C2×C2, set q(u_i)=e_i u_i and q(v_i)=f_i v_i. The pair's source-to-target defect is e_i−f_i=e_i+f_i=d_i. This lift gives bounded **margin-neutral** original SCC blocks of size two and nonzero increments in fixed Q, yet unbounded required events. It directly refutes dropping the packet-cycle condition L from Corollary 3a's finite-group bound, without appealing to unpaired margin obstructions.

## 6. Scope and stopping boundary

These theorems classify one fixed binary endpoint interval. They do not prove that an arbitrary receiver is Boolean, that an arbitrary full predicate has arity two, or that a native current bank realizes the abstract charge labels. They do not add capacity by naming more kernel lifts. They do not infer an authenticated microscopic transition from two accepted arrivals. They do not establish a universal move bound from bounded separator width, arity, original SCC size, finite charge group, or exact total population.

The structural replacement is exact: retain the implication order **and** defects, then use zero-defect relative ideals or neutral packetization with acyclic dependency contraction. Small packets become a valid bound only after dependency cycles and literal physical costs are controlled. All-size claims rest on the proofs; finite checks only audit the explicitly recorded universes and hostiles.
