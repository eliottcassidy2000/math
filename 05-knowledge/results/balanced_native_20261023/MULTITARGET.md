# Directed target cycles followed by finite charge packets

**PROVED, conditional native extension.** The two-target packet rule extends to t targets in a common paired-charge fibre. For a finite abelian p-group G, every accepted endpoint return in a jointly compatible supplied inventory has a monotone native decomposition with atomic support at most

    t D(G^(t-1)) = t [1+(t-1)(D(G)-1)].                 (1)

Every differing slot changes once. Same-target path changes may be executed separately with support one. This is an upper bound; sharpness is not asserted for t>2. The targets, physical-row map, shared records and unused normalized capacities are part of the hypotheses, not outputs of the formula.

## Interface and construction

Let the t distinct target owners have one common paired-charge label in G. All source owners are outside this target set. Each labelled slot has a supplied old and new full path record, including actual pins, owner names, kernel lift and normalized resource. Suppose every endpoint mixture is independently compatible with the fixed shared record, locally guarded, and resource-injective against the unchanged remainder. The two endpoints have equal target margins and equal original physical-row current signature. All source weights remain attached to their original slots; no new stationary sampling law is inferred.

For each differing slot draw a directed edge from its old target to its new target, labelled by its actual source charge. Target-margin equality makes this a balanced directed multigraph. Decompose it into simple directed cycles, retaining every edge's actual slot index. Loops are already harmless same-target changes. A nonloop cycle uses at most t slots and preserves every integer target margin when applied atomically.

For a cycle C let delta_C(v) be the sum of incoming source labels at target v minus the sum of outgoing source labels. Each source row is safe because the target labels agree. The target rows require delta_C=0, and

    sum_(target v) delta_C(v)=0.

Consequently each cycle carries a charge in

    Gamma = ker[G^t -> G, (z_1,...,z_t) -> sum z_i]
          isomorphic to G^(t-1),

with the final coordinate reconstructed as the negative of the retained sum. This is an exact physical-row sidecar, not an unowned scalar cancellation. The total cycle charge is zero because the endpoint currents agree. Apply the Davenport extraction theorem to the sequence of cycle charges. A packet of at most D(Gamma) cycles has zero target-row charge and at most t D(Gamma) slots. Execute it, remove its cycles, and repeat. No slot is reused; every intermediate state has the original margins and guard. All-mixture compatibility supplies the actual shared owners and distinct ports at each arrival.

For a p-group the classical formula gives D(Gamma)=1+(t-1)(D(G)-1), proving (1). For a general finite abelian group the argument still gives t D(G^(t-1)); the p-group formula is not extrapolated to mixed prime groups. For t=1 no nontrivial target change exists and every supplied same-target record change is safe under the compatibility hypotheses.

For the original full dyadic code D(G)=132 this yields264 for t=2,789 for t=3, and1576 for t=4. The bound is independent of the number of supplied slots at fixed G,t. It need not be independent of tournament order, because its paired-charge group and available targets can change.

## Why cycles and group coefficients play different roles

Integer target margins cannot be replaced by congruences: a nonzero multiple of a modulus is still a wrong count of actual target paths. Directed cycles first impose these integer conservation laws exactly. Finite group extraction then cancels the remaining target-owned harmonic debts. Finally all-mixture compatibility turns the algebraic packet into actual paths. The three steps preserve different data and are needed in that order.

The concurrent `native_packet_decomposition_20261023/interface/ENDPOINT_FILTERS.md` handles the failure of independent mixing by a shared-record/resource implication poset. Its acyclic vacancy-chain hostile shows that even tiny charge group and singleton implication components can force unbounded packets. Thus (1) cannot survive by merely attaching a poset annotation while continuing to choose arbitrary cycle subsets. A proposed extension must extract charge-neutral *feasible filters*, or supply compensating excursions outside that endpoint interval.

## A literal three-target fibre, and the missing native coordinate

The complete original order-eight label table gives exactly

    {570,1965,2234}

at paired label (1,1,0,1;18). This is an actual three-target equal-charge fibre. However the complete common-clock domains of the current28-source supply contain no path to570. They contain the56 U/V paths and one extra path to1559. The arithmetic third target therefore supplies no new route under the retained pins. It is a precise target for a richer source inventory, not a licensed extra endpoint.

`multi_target.py` verifies this entire label fibre and the complete57-path inherited domain, then independently tests the cycle-plus-charge algorithm on seeded abstract inventories for groups C2+C4 and C3 and t=2..5. Every slot remains indexed; every intermediate target count and every target-owned charge is recomputed directly. These experiments are controls for the conditional proof, not native realizations of the arbitrary abstract labels.
