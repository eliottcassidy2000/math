# Non-reacquiring repairs: endpoint intervals, neutral trades, and native obstructions

24 October 2026. **PROVED, conditional all-size:** exact binary no-reacquisition routing criterion, neutral-difference partitions, and disjoint-block antimatroid sufficient condition with literal costs. **PROVED, native all-size supplied carrier:** paired prefix keys realize these binary relations without counting kernel lifts as capacity. **FINITE-EXACT:** literal source, path, partner, converse, reservoir and route controls. This is not a new arbitrary-order electrical bank or a proof that all receiver repairs satisfy strict temporal nonreuse.

## 1. Four different resource contracts

An arrival is a full accepted assembly, not a charge label. Keep source and target matrices, Hamiltonian paths, kernel lifts, normalized resources, inverse/converse records, original pins, population and every required current. An atomic event replaces full records between accepted arrivals. Unless separately supplied, no microscopic edge-by-edge implementation or intermediate-arrival semantics is asserted.

The following contracts are different.

1. **Capacity at arrivals:** no normalized unit is occupied twice at an arrival.
2. **Continuous occupancy:** the same unit remains assigned to its slot throughout a route. In the atomic record model this means the resource label never changes; physical continuity inside an unspecified implementation is not additionally claimed.
3. **No resource reacquisition:** the times at which each normalized unit is occupied form an interval. A released unit is never used again. Across-slot handoffs require explicit semantics; the theorems below use disjoint slot reservations and have no such handoffs.
4. **Per-slot full-record non-return:** once a slot leaves a particular full record, it never visits that record again. Continuing to hold a record is not a return. Distinct kernel lifts can share a normalized unit.

Continuous occupancy implies no resource reacquisition but not record non-return: a fixed-key phase may go 0,1,0. Conversely, record non-return need not imply no resource reacquisition: records (key a,phase 0), (key b,phase 0), (key a,phase 1) are three different records but reuse key a. Repeating an entire assembly is a still weaker prohibition than per-slot non-return. A capacity unit is never created by giving its kernel lifts different record names.

In the **binary disjoint-reservoir contract**, slot i has exactly two accepted full records r_i^0,r_i^1, their resources are distinct, and all 2n reserved resources are pairwise disjoint. Then no resource reacquisition and full-record non-return coincide. A route satisfies them precisely when every slot switches at most once. The original order-eight charged U/V carrier has this type of alternative reservation; the new native carrier in Section 6 independently authenticates another instance, with different guard scope.

If both records instead share one fixed key, the same binary criterion describes record non-return, but all routes already satisfy no resource reacquisition. Do not substitute that weaker fact for the binary disjoint-reservoir theorem.

## 2. Exact endpoint-interval criterion under arbitrary relational pins

Let F be any fully authenticated accepted subset of the binary product, after all relations, pins, exact charges, integer margins and capacity constraints. No factor may be silently omitted. Fix x,y in F and put S={i:x_i≠y_i}. For I⊆S let z_I agree with y on I and with x elsewhere. Retain

    L(x,y) = {I⊆S : z_I∈F}.

For each differing slot let c_i be its actual changed-owner-clock count between the two full records. Here c_i>0. The event I→J, I⊊J, has the exact cost pair

    (|J\I|, sum_{i∈J\I} c_i).                       (1)

Add any separate helper/visible occurrence work only if those occurrences were not already included among these coordinates. Double counting or forgetting a receiver coordinate is not permitted.

**Endpoint-interval theorem.** Under the binary disjoint-reservoir contract, strict no-reacquisition routes from x to y are exactly directed paths from ∅ to S in the DAG of retained subsets ordered by inclusion. The same statement holds for record non-return with binary fixed-key records. Delete edges violating either the slot or clock event budget. Reachability then decides strict routing exactly. A minimum-bottleneck dynamic program gives the exact minimum of either event cost; fixed-budget DAG shortest paths give the exact minimum number of events. The set of attainable simultaneous slot/clock budgets is the upward closure of the DAG paths' bottleneck cost pairs.

**Proof.** An endpoint-fixed slot cannot leave its initial record: it would have to return. A differing binary slot must switch from its initial record to its final record exactly once. Thus the set of completed switches increases strictly at every genuine event, and every arrival is z_I for a retained I. Conversely every increasing retained-subset path switches each slot at most once and contains only accepted arrivals. Disjoint alternative reservoirs convert this coordinate fact into global resource non-reacquisition. Literal occurrence counts add over disjoint slots, proving (1). QED.

Every such route changes exactly |S| slot occurrences in total, costs exactly sum_{i∈S}c_i owner-clock occurrences in total, and uses at most |S| genuine events. This total-work equality does not claim a small event bottleneck. A direct endpoint event always exists in the feasible-arrival model, with cost (|S|,sum c_i); an additional transition predicate, if required, must be inserted into the DAG and may remove that edge.

This is an exact finite criterion and an all-size characterization, not a polynomial-time bound: the endpoint interval has up to 2^|S| vertices. It preserves arbitrary full-record relations without projecting away their scopes.

The stated route starts at x with no earlier releases within its declared nonreuse contract. If it continues an existing history, carry that history forward and additionally forbid every acquisition of an already retired resource or return to an already departed record, as applicable. The same DAG with those changes removed remains exact, but a direct endpoint event may then be forbidden. Starting another endpoint query does not reset a user's continuing nonreuse constraint.

### Conditional receiver sections must be endpoint/history conditioned

Partition S into frame switches S_H and receiver switches S_R. A retained state is a pair (U,V), U⊆S_H, V⊆S_R. Its conditional receiver section is

    N_U = {V : z_{U∪V}∈F}.

Keep these actual sections, not just their charge supports. An edge from (U,V) to (U',V') must satisfy U⊆U', V⊆V', at least one strict, and the exact combined costs (1). This is exactly the same DAG, so it is an exact section compiler. With a prescribed increasing frame sequence, propagate all reachable V values, including allowed increasing within-frame moves before each crossing. A nonempty section at each frame is not sufficient; one compatible increasing witness chain must survive.

For a crossing between two such sections the exact endpoint-conditioned cost is the minimum over actual compatible pairs V⊆V', retaining the pair and charging frame and receiver changes together. But independently minimized crossings still cannot be composed unless their middle witnesses are connected in the required direction within that section. Ordinary undirected fibre connectivity does not establish this. Keeping the whole directed state/section relation, or an authenticated coherent monotone section, is sufficient.

In particular the old hole-difference existence rule can certify an accepted branch witness without certifying a temporal route to that witness. To make a no-reacquiring repair theorem, restrict its branch fibres to this endpoint/history interval and retain compatible witnesses and their costs. No general receiver refresh-and-restore construction meets the contract merely because each refresh is individually capacity valid.

## 3. Neutral-difference partitions when only additive guards remain

Suppose after disjoint record/reservation authentication the only restrictions are an additive guard g(x)=g(y) in an abelian group G, possibly including integer populations as Z coordinates. Let

    d_i = g_i(y_i)-g_i(x_i), i∈S.

Here G may be a complete authenticated paired dyadic code **only where its equivalence to the required physical guard was separately proved**. Otherwise use the original complete additive guard. Exact fixed integer populations must not be dropped.

For I⊆S, z_I is accepted exactly when sum_{i∈I}d_i=0. Therefore a strict route is exactly an ordered partition

    S = D_1 ⊔ ... ⊔ D_m, each D_j nonempty,
    sum_{i∈D_j}d_i=0.                                (2)

Its event costs are (|D_j|,sum_{i∈D_j}c_i). Conversely every such partition may be executed in any order and is strict. Thus the optimum bottleneck is the optimum zero-sum partition bottleneck, not a span or lattice test using trades outside S.

Every nonempty neutral block can be partitioned recursively into inclusion-minimal nonempty neutral blocks. Refinement never increases either event cost. Hence an optimum can be attained by a partition into neutral atoms. If additional relational predicates remain, their allowed prefix states must be checked, and refining a block may introduce forbidden arrivals. Neutral-atom refinement is not asserted then.

This proves a meaningful sufficient repair condition: a proposed hole compensation must be expressible as disjoint endpoint-aligned neutral trades with legal prefixes. Using a compensation coordinate and later restoring it is disallowed even if the final charge cancels exactly. Kernel lifts cannot be substituted for fresh compensation units.

### Uniform finite-group bound after exact population pairing

There is an all-size additive-only bound stronger than a generic lattice assertion. Suppose binary slots belong to disjoint population classes, and each class has an exact fixed V-count at both endpoints. Within each class, pair every endpoint switch 0→1 with one switch 1→0. These are disjoint two-slot packets preserving all integer margins. For overlapping margins, this pairing is valid only when each paired pair has the same complete margin column; arbitrary overlapping margin rows do not have that property automatically. Alternatively supply an authenticated margin-neutral packetization explicitly.

Let the remaining complete charge guard take values in a finite abelian group Q, and let D(Q) be the least integer such that every sequence of D(Q) elements has a nonempty zero-sum subsequence. The charge increments of the two-slot packets form a zero-sum sequence. Repeatedly select a zero-sum subsequence of at most D(Q) packets: when at least D(Q) remain use the definition on any D(Q) of them, and when fewer remain use the entire remaining zero-sum sequence. Remove each selected subsequence. This gives disjoint neutral endpoint-aligned events, each changing at most

    2 D(Q) slots.                                   (9)

Execute them in any order. Every population and charge guard remains exact; every differing physical slot switches once. The exact event clock cost is the sum of its actual changed-slot costs, bounded by the sum of the largest min(2D(Q),|S|) costs among the differing slots, or by 2D(Q)c_max. Total work is still precisely the endpoint work. The complete reservation must be disjoint as in Section 1, and there can be no other relational predicate unless closure under these intermediate trades was separately proved.

Witness extraction is constructive: on the chosen at-most-D(Q) packet charges, maintain a table of reachable group sums and one actual subset witness for each. Before inserting a charge g, a previously reached −g supplies a nonempty zero-sum subset after adjoining g; otherwise extend the table by g. The classical bound proves such a subset appears. Repeating on the remaining packet list produces actual disjoint events, not merely a charge label. The table has at most |Q| entries, so no enumeration of all assembly states is required for this additive-only algorithm.

For Q=Z/128×(Z/2)^4, the classical finite abelian p-group zero-sum theorem gives

    D(Q)=1+(128−1)+4(2−1)=132,

so (9) is 264 slots. This is a conditional inventory-size-independent bound, **not** a better bound for the 82-slot inventory, where direct endpoint replacement is already smaller, and not an authentication of a physical Q guard at other tournament orders.

For completeness, a self-contained proof of the group constant follows. If Q=∏_{j=1}^s Z/p^{a_j}, its group algebra over F_p is

    F_p[Q] ≅ F_p[t_1,...,t_s]/(t_j^{p^{a_j}}),

under t_j=X^{e_j}−1. Its augmentation ideal has nilpotence index D=1+∑_j(p^{a_j}−1): every monomial of degree D vanishes. Given a length-D zero-sumfree sequence g_1,...,g_D, the product ∏(1−X^{g_i}) lies in that vanished ideal power. Yet its coefficient at the identity is 1, since only the empty subsequence could sum to zero. Contradiction. Conversely p^{a_j}−1 copies of each independent basis generator e_j give a zero-sumfree sequence of length D−1, because each coordinate count is strictly below its order. Hence D(Q)=D exactly. This classical theorem is attributed to John E. Olson, “A combinatorial problem on finite Abelian groups, I,” Journal of Number Theory (1969), [publisher record](https://www.sciencedirect.com/science/article/pii/0022314X69900213); the proof here is self-contained, and the full publisher article was not consulted.

The auxiliary characteristic-p group algebra is only a zero-sum proof device. It does not replace actual nonnegative occupancy, turn a phase into capacity, or change the characteristic of any physical guard.

The factor 264 is sharp in the abstract additive theorem when arbitrarily many exact population classes are allowed. Take 132 disjoint two-slot classes, each constrained to exactly one V. Use the length-132 exchange-increment sequence consisting of 127 copies of e_0, one each of e_1,e_2,e_3,e_4, and e_0+e_1+e_2+e_3+e_4. It has total zero in Q and no proper nonempty zero-sum subsequence: without its last term, independent coordinate counts cannot vanish; with it, all four binary generators and all 127 cyclic generators are required. Give the initially-U slot in each pair its indicated charge increment, and the initially-V slot zero increment. Population constraints force each pair either to stay or to exchange fully. The zero-charge fibre therefore has exactly two states, differing on all 264 slots. This sharpness example is an **abstract authenticated-data model**, not a realization of those charge labels by native tournament currents and not sharpness with only the two core/helper populations.

This additive corollary does not repair arbitrary relational pins. The two-point relation in Section 6 can forbid every proper neutral event even when its charge is zero. The useful division is: authenticated additive-only guards give neutral partition routes and this uniform bound; additional relations require the endpoint DAG, a closed-trade proof, or a structural condition such as Section 4.

## 4. A structural all-size sufficient condition with exact costs

Partition selected physical slots into nonempty disjoint blocks E_1,...,E_m. Each block has a declared initial and final tuple, with **every slot in the block different**; other slots are fixed. In the family considered, a block is always entirely initial or entirely final. Write

    s_j=|E_j|, w_j=sum_{i∈E_j} c_i.

Supply a family A of allowed block subsets. Each subset must reconstruct an actual accepted full assembly with all original guards and capacities. A is an **antimatroid** if ∅∈A, it is closed under unions, and every nonempty member can lose one element and remain in A. Remove unused ground elements, so the union of A is {1,...,m}. No claim is made that arbitrary relationally constrained receivers have this property.

**Block-antimatroid theorem.** For any A,B∈A there is an accepted route switching each block in A△B exactly once and no other block. It satisfies record non-return, and no resource reacquisition whenever the two block states use disjoint globally reserved per-slot alternatives (or continuously fixed keys). It uses at most |A△B| events and has exact total work

    slot work = sum_{j∈A△B}s_j,
    clock work = sum_{j∈A△B}w_j.                      (3)

If the accepted family is exactly this synchronized block family, its exact endpoint event-budget feasibility is

    slot budget ≥ max_{j∈A△B}s_j,
    clock budget ≥ max_{j∈A△B}w_j,                  (4)

with maxima 0 for identical endpoints. Both budgets in (4) can be attained simultaneously. In particular the global thresholds are max_j s_j and max_j w_j. This is also the optimum among ordinary routes in this declared family, not merely among strict routes.

**Proof.** Accessibility gives a one-element-at-a-time chain ∅→B. Union every member with A. Union closure gives an accepted chain A→A∪B, changing each element of B\A once after deleting duplicate states. Likewise a chain ∅→A, unioned with B and reversed, gives A∪B→B, changing each element of A\B once. These sets of changes are disjoint. Each physical slot changes exactly once or never, proving non-return and (3). One-block events prove the upper bounds in (4). Conversely every differing block must change somewhere, and synchronization forces that event to pay at least s_j slots and w_j clocks. This proves both simultaneous lower bounds. QED.

Order ideals of any finite poset give a simple verifiable subclass: they are union-closed and remove a maximal element. Independent disjoint neutral trades give the Boolean-cube subclass. More generally antimatroids allow dependencies that are not all order-ideal dependencies. This is a sufficient structural criterion, not an equivalence or a claim about unverified native wiring.

For additive guards, if the entire antimatroid family is actually accepted, each block is guard-neutral: accessibility gives a chain from ∅ to the full ground set, and the one-block differences on that chain preserve the guard. One may instead establish neutrality first and then verify only the remaining relational antimatroid condition. Either way all guards remain authenticated.

As throughout the constructive statements, all accepted endpoint pairs permit atomic replacement before cost limits. Additional directed or microscopic transition restrictions must explicitly allow the constructed one-block events; they cannot be inferred from antimatroid membership alone.

## 5. The precise failure of undirected fibre lifting

Take binary states

    F={000,001,101,111,110}, x=000, y=110.            (5)

At one-slot budget, the induced accepted graph is exactly the five-state path in the displayed order. Every route from x to y must visit a state with third bit 1 and later restore third bit 0. The endpoint interval fixes third bit 0 and retains only x,y. Consequently the ordinary endpoint threshold is one slot, while the strict no-return/no-reacquisition threshold is two slots.

Project to the first bit as frame mode. The two fibres {000,001} and {101,111,110} are each internally one-slot connected, and their exact minimum boundary cost is one slot, attained by 001→101. Thus even **connected fibres plus exact witness-retaining boundary minima** suffice only for ordinary connectivity. They do not suffice for strict temporal routing. The restricted endpoint/history sections have only 000 and 110 and exact boundary cost two slots.

This example has no additive-charge obstruction: a constant guard accepts every state before the declared full-record relation. Section 6 makes every record and resource in (5) native, so the distinction is not merely a symbolic phase example.

## 6. Actual paired prefix keys realize the binary theory

Inherit the authenticated native carrier S_r=T_r→C3, q=8r, from the byte-pinned source archives. Its order-eight rigid seed has canonical code 29558785 and the endpoint-marked prefix paths

    P0=(3,0,6,7,1,5,4,2),
    P1=(3,1,6,7,0,5,4,2).

For key b∈{0,...,2^r−1} let the j-th prefix block use P_{(b>>j)&1}; attach the suffix q+k,q+(k+1 mod3),q+(k+2 mod3). Call the full record R_b(k). The inherited all-size carrier authenticates the source, K partner, inverse/converse square, endpoint-preserving gluing g=(4,0,2,3,1,5,6,7), and exactly one normalized four-word resource for each prefix key. The three kernel phases at one key are three full records, not three resources.

Choose n≥1 and r≥1 with 2n≤2^r. Slot i, indexed i=0,...,n−1, is supplied the six full records

    R_{2i}(k), R_{2i+1}(k), k∈F3.

Retain all six in the offered menu; impose the explicit phase pin at the actual suffix clock q:

    owner(q)=q, equivalently k=owner(q)−q=0.         (6)

It is **not** a raw-owner residue-zero pin when q is nonzero modulo3. A further arbitrary full-record relation Γ on the selected key bits produces exactly the accepted binary family Γ. The two keys differ only in their least significant bit, so precisely the first prefix block changes from P0 to P1. The paths differ at clocks 1 and 4 and nowhere else. Hence each binary switch costs exactly one slot and **two owner-clock occurrences**, independently of r.

The 2n alternative keys are distinct, so their normalized reservoirs are pairwise disjoint. Every arrival occupies n units/4n closed words, from a fixed reservation of 2n units/8n closed words. Kernel lifts add no capacity. Any immutable background must avoid the complete reservation. Target class and source class are constant across all these records, authenticated by the native gluing maps, not assigned by name.

Consequently every additional retained guard that truly factors through this fixed source/target class incidence is constant on the family, and must be checked once. A phase-sensitive, owner-sensitive, or different physical guard is not inferred constant. In particular this does **not** import the order-eight numerical aggregate-dyadic bank into order 8r+3, does not establish arbitrary all-order current invariance, and does not construct a single global Hamiltonian observer.

### Literal Boolean-predicate embedding, not universal repair

The selected key bit is read from an actual prefix clock without a field-coordinate convention:

    bit_i = owner_i(1),                              (8)

because P0 has owner 0 at clock 1 and P1 has owner 1 there. Thus any explicitly supplied Boolean predicate Φ, including a CNF formula, can be imposed as Φ(owner_0(1),...,owner_{n−1}(1)). Together with the offered paired-key menus and phase pin (6), its satisfying Boolean assignments correspond bijectively to accepted native full-record assemblies. No auxiliary records, phase multiplicity, or new capacity are introduced. For a formula input use r=ceil(log2(2n)); the source has 8r+3 owners, the 2n offered key sidecars and their paths have polynomial size, and copying the formula preserves its exact truth table. For n=0 use the trivial empty assignment separately. This is a parsimonious explicit Boolean-CSP embedding, not a claim that arbitrary predicates can be compiled from equality-only constraints or that they admit bounded-width separators or a common observer. It explains why arbitrary relational pins must retain their actual combinatorial difficulty rather than being declared repaired from charge coverage alone.

Converse is literal: transpose S_r and reverse every complete path. The same canonical phase is read at reflected clock 2; the pin in (6) becomes owner_reverse(2)=q. If instead using the converse head, its phase chart is shifted by 2. The two changed prefix clocks become N−1−1 and N−1−4 for N=q+3. Thus exact costs, normalized resources and no-reacquisition histories are preserved. The checker verifies K(S_r^op,R^reverse)=K(S_r,R)^op on actual owner labels.

**Native strict-repair separation.** With n=3,r=3 and Γ from (5), ordinary connectivity holds at one slot/two clocks. Every such route between the stated endpoints reacquires normalized key 4 after releasing it for key 5, and returns to the initial full record at slot i=2. The exact strict threshold is two slots/four clocks, achieved by direct 000→110. The unconditional boundary minimum between first-bit fibres is one slot/two clocks; the strict endpoint-conditioned minimum is two slots/four clocks.

**Native positive all-size family.** Partition any n chosen binary slots into blocks and impose the full-record relation that the key bits are synchronized within blocks and their active block set belongs to a supplied antimatroid. Section 4 applies to actual native full records with s_j=|E_j| and w_j=2|E_j|. It provides strict no-resource-reacquisition and per-slot record non-return with exact endpoint thresholds and total costs. This is a declared higher-arity relational family, not a claim that every such predicate follows from literal indexed-owner equalities.

**Native obstruction to universal bounded-locality repair under arbitrary relations.** Impose Γ={0^n,1^n}. Both arrivals are valid, have identical class incidence, and retain the same offered kernel lifts and disjoint alternative reservation. The only genuine event changes all n slots at exactly 2n clocks. Thus no bound independent of n can hold for arbitrary relational pins in this family. The all-size positive theorem is conditional on the proved structural property, not a blanket arbitrary-relation supply theorem.

## 7. Strict temporal routes are not stationary mixing processes

For finite per-slot menus M_i, every strict full-record non-return route changes slot i at most |M_i|−1 times. The total number of changed-slot occurrences, and hence the number of genuine events, is at most

    sum_i (|M_i|−1).                                (7)

Retain the current assembly and each slot's set of already visited full records. Every genuine strict move increases the rank sum_i |visited_i|. This history-augmented transition graph is a DAG. No entire assembly can recur on a strict route; that would return each changed slot to an earlier record. Allowing idle repetitions does not create genuine moves.

On a finite history-augmented Markov state space with transitions respecting strict non-return, any stationary distribution is supported where the probability of a genuine move is zero. Indeed the expected rank increment in stationarity is zero, while every genuine move has positive increment. Thus an indefinitely active reversible sampler on nontrivial record states cannot at the same time satisfy this strict contract. A resetting history starts a different contract; it does not prove non-return across resets.

Continuous fixed-resource occupancy is different: reversible phase sampling may repeatedly change records and still never release a normalized unit. Stationary or spectral statements for such a process concern that weaker resource contract. Likewise ordinary receiver refresh mixing cannot be re-labelled strict full-record or resource non-return when released alternatives are revisited. No new claim about the separately incoming weighted-projection results is made here; the point is the precise temporal hypothesis required before using any such claim.

## 8. Verification and scope ledger

Run `python -B check.py`, `python -O -B check.py`, `python -B verify_manifest.py`, and `python -B portable_check.py`. All checkers use explicit exceptions rather than removable assertions, the standard library only, no network, no repository imports, and no parent-packet execution. Source archives are exact byte copies under `.txt` names, with original paths and SHA256 pins in SOURCEPINS.json. MANIFEST.json closes every shipped file except itself; manifest-only replay copies only those files into a fresh temporary directory and reruns both Python modes.

Finite controls exhaust every accepted binary family through three slots, compare the endpoint DAG against an independent literal visited-record search for every accepted ordered endpoint pair and every nonuniform weighted budget, exhaust antimatroids through four ground elements and test all endpoint geodesics, and check neutral partition controls with integer and mod-three guards. Davenport controls exhaust small p-group sequences, check the exact 2048-dimensional dyadic augmentation calculation, count the two zero-sum subsets of the sharp length-132 sequence, and enumerate a smaller exact population-pair fibre. Native controls reconstruct source/target matrices and all offered paths, words, converse records, phase pins, disjoint reservations and exact costs at r=1,...,4, plus the hostile and positive routes. These finite controls supplement the proofs and inherited all-size carrier, rather than extrapolating all-size supply from samples.

Loss ledger: charge support alone loses compatible monotone witnesses; a history-free quotient loses the no-reacquisition constraint; a kernel phase label alone loses the distinction between record and capacity; abstract linear changes of variables need not preserve literal clock costs; constant class incidence does not authenticate an unrelated electrical bank; and no local record theorem supplies a global observer without additional proof.
