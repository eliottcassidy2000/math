# Literal C3 owner relations: holonomy, indivisible packets, and exact fibre bridges

24 October 2026. **PROVED, scoped all-size:** the equality compiler, flat-component packet classification, weighted-frustration boundary formula, and conditional neutral-fibre lifting theorem. **FINITE-EXACT:** independently enumerated owner-record controls in `check.py` and `receipt.json`. No global observer or unrestricted electrical routing claim is made.

## 1. Native carrier and exact scope

The pinned source is commit `c09b642a629b7d12b529b88c2697abaee78e73e7`. The inherited `NEUTRAL_FIBRE_LIFTING.md` §7 warns that a silent C3 kernel can force a large coupled event. `ONE_BIT_WORKSPACE.md` §3 supplies the actual carrier, and `inputs/CUBES.txt` supplies its endpoint-marked partner gluing. Exact source copies and their original paths are in `inputs/` and `Sourcearchives.txt`.

Let S_r=T_r -> C3, of order q+3, q=8r. Its rigid prefix blocks use

    P0=(3,0,6,7,1,5,4,2), P1=(3,1,6,7,0,5,4,2),
    seed canonical code 29558785,
    endpoint-marked partner map g=(4,0,2,3,1,5,6,7).

For slot i choose a fixed prefix key b_i in {0,1}^r, pin its entire prefix, and retain the complete three-record menu

    R_i(k)=P_(b_i) followed by (q+k, q+(k+1 mod3), q+(k+2 mod3)),
    k in Z/3.

All labels in this note are literal owners in the SAME source frame, not independently named owner copies identified by fiat. Aut(S_r)=C3 rotates the terminal cycle and fixes the prefix. The three k values use exactly ONE normalized K/converse resource b_i. Different b_i have disjoint four-word resource reservoirs. Consequently capacity is equivalent to distinct b_i, plus avoidance of any fixed occupied background resource. If b_i=b_j, these two slots cannot coexist, even with different phases. All phase assignments are capacity-valid when keys are distinct; no argument below increases capacity by counting lifts.

The endpoint-fixed maps g, extended by the identity, authenticate prefix-coordinate partner isomorphisms. Literal final C3 rotations are source automorphisms and, by K equivariance, actual partner isomorphisms for kernel changes. Each nontrivial kernel change alters exactly three owner-clock occurrences. Every supplied record has the same target isomorphism class. Thus every target-class-dependent per-slot character and every selected class-incidence/current *difference* is constant. This does not certify arbitrary owner-sensitive electrical conditions or transport an order-eight numerical harmonic bank to order q+3.

The first theorem permits any conjunction of literal equalities between indexed owner occurrences and pins to specified literal owners. It does not permit inequalities, disjunctions, arbitrary binary phase relations, or unverified target-side attachments. Atomic moves replace a packet of full records between accepted arrivals; they are not asserted to be individual edges in an older metagraph.

## 2. Exact equality compiler and holonomy theorem

At final clock q+t, slot i has owner q+(k_i+t mod3). Therefore

    owner_i(q+t)=owner_j(q+u)
      iff k_j-k_i=t-u in Z/3.                         (1)

A prefix/prefix equality is either constantly true or constantly false, because both prefixes are pinned. A prefix/final equality is constantly false: the owner sets are disjoint. A pin at a prefix clock is likewise a constant test; a final pin to q+v is k_i=v-t, and a final pin to a prefix owner is impossible. These rules are exhaustive, including equalities between occurrences in the same slot.

Build a directed voltage multigraph on slots with edge i->j labelled d=t-u for each equation (1). Its reverse has label -d. Parallel edges and loops are allowed. Delete true constant predicates; a false constant predicate rejects the assembly. Each unary phase pin is retained.

**Holonomy theorem.** Fix one root in each connected component C. The equations have a solution iff:

1. the signed voltage sum around every closed walk is zero;
2. after defining p_i as the voltage sum from the root to i, all pins in C specify the same root phase t_C=k_i-p_i.

When these conditions hold, every solution is uniquely

    k_i=t_C+p_i, i in C,                               (2)

with t_C fixed by pins or arbitrary in Z/3 if C is unpinned. Hence there are exactly 3^f accepted phase assignments, where f is the number of unpinned components. This count assumes the independent capacity and constant-predicate checks already passed.

**Proof.** Summing edge equations around any closed walk gives necessity of zero holonomy; substituting (2) into pins gives their compatibility. Conversely zero holonomy makes each p_i path-independent. Formula (2) then satisfies every edge, and precisely the stipulated pins. Any solution has this form because propagation along a root-to-i path determines k_i-t_C. QED.

Testing a spanning forest and each non-tree edge suffices. The compiler may return either the potentials and root domains, a nonzero cycle voltage, a contradictory pair of propagated pins, or the false constant/capacity predicate. Each is a finite checkable witness. On a forest without incompatible pins, there is never a holonomy obstruction.

**Realization converse.** Every finite directed Z/3 voltage multigraph and every list of phase pins is realized by actual native records with distinct resources: choose r with 2^r at least the number of slots, choose distinct b_i, realize i->j of voltage d by equality of clocks q+d at i and q at j, and realize k_i=v by pinning i's clock q to owner q+v. This realization uses no invented owner correspondence and retains the actual partners above.

The realizable full two-slot relations from equality-only constraints are precisely: the full 9 pairs if there is no nonconstant comparison; one of the three translation perfect matchings if consistent comparisons are present; or the empty relation if inconsistent. Pins intersect these relations in the stated unary way. Thus “arbitrary equality relations” does not mean arbitrary subsets of the nine phase pairs.

### Affine coordinates do not add native reflection holonomy

If slot i is merely reparametrized by k_i=a_i z_i+c_i, a_i in {1,-1}, then (1) becomes

    z_j=(a_i/a_j) z_i + (d+c_i-c_j)/a_j.

All edge slopes are vertex-gauge ratios; their product around a cycle is 1. Such reparametrization cannot create reflection holonomy. General AGL(1,3) edge systems admit a broader classification by common fixed points of root holonomies, including a reflection pinning one root value, but arbitrary edgewise affine maps are NOT authenticated by the native common-source equalities here. A claimed reflection cycle must supply extra non-native port data, rather than silently relabeling each incidence independently.

## 3. Mandatory packets and the charge-neutral factor rule

For any two accepted assignments k,k', subtraction of (1) gives

    k'_i-k_i = k'_j-k_j on every edge.

Thus their difference is constant throughout each component. A pinned component never changes. Every unpinned component either does not change at all or changes all of its slots by the same nonzero phase. Conversely each such component rotation is accepted. This exactly classifies mandatory coupled packets, not merely a lower bound.

Assign slot weights w_i>0 and define W_C=sum_(i in C) w_i. For changed-slot cost w_i=1; for actual owner-clock occurrence cost w_i=3. A packet changing components D has exact cost sum_(C in D) W_C. With no additional guard, the minimum budget joining two prescribed assignments is the largest W_C among components on which they differ (zero for identical assignments): change those components one at a time for the upper bound, and use indivisibility for the lower bound. Distinct C3 phase choices do not create vacancies or weaken that bound.

Now retain a complete additive guard in any abelian group H, including integer coordinates, by assigning the actual slot record charge A_i(k). Coincident physical owner rows must be combined as required by the original guard; H must retain the entire required current/margin vector, not an unproved scalar substitute. Define

    Q_C(t)=sum_(i in C) A_i(t+p_i).

The exact feasible root assignments satisfy sum_C Q_C(t_C)=the required total, with pinned roots retained. Call C charge-neutral when Q_C is constant on all three phases. This allows cancellations within C, but requires their exact proof. A target-class-only guard on the native carrier above is neutral even slot by slot.

For budget B retain every nonneutral unpinned component and every unpinned component with W_C>B as visible. Remaining neutral light components are independent C3 factors and may each rotate within budget B. The visible quotient has the exact charge equation and an edge when the sum of W_C over changed visible roots is at most B. Full B-connectivity is exactly connectivity in this quotient; each fibre is a product of light neutral factors, with a route changing each differing factor once. Neutral but heavy components cannot be erased.

This is an all-size conditional rule, not an assertion that arbitrary full current factors are neutral. It also does not solve the residual charged-root graph automatically.

### Exact spectrum of distances and budget components in the neutral family

For the fixed equality family with no additional phase-varying guard, choose any accepted reference assignment. Write m_C for the number of slots in an unpinned component C. The exact joint enumerator of changed slots and owner-clock occurrences is

    E(z,t)=product_(C free) (1+2 z^(m_C) t^(3m_C)).       (6)

A pinned component contributes 1. A frustrated or capacity-invalid assembly has no accepted reference and its unbased accepted-state enumerator is 0; formula (6) is only for a nonempty fibre. The factor follows because a free C has one unchanged phase and two alternative phases, each changing exactly its m_C slots and 3m_C clocks. Independence of free roots proves the product. A phase-varying joint charge guard generally destroys this product and must instead be retained as above.

For owner-clock budget B, the full packet graph has exactly

    3^(number of free C with 3m_C>B)

connected components. Their labels are the phases of the heavy components; light components each rotate within budget. Allowing extra multi-component packet edges cannot change a heavy component because all costs are nonnegative. Every connected component has 3^(number of light free C) states. The connectivity threshold is 3 max_C m_C, with maximum0 for a singleton fibre. The analogous slot-budget formula replaces 3m_C by m_C.

### Resource non-reuse is not full-record non-return

Every phase replacement in this lane keeps slot i on the SAME normalized resource b_i continuously. It neither releases b_i nor acquires another resource. Hence all our fixed-key routes satisfy the strict predicate “never reacquire a normalized resource after releasing it”: there is no release at all, and distinct slots never share a unit. This is not a claim of new capacity from the three kernel lifts.

A different predicate, “never return to a full record after leaving it,” concerns literal records R_i(k), not normalized units. In the fixed flat family, the direct endpoint route rotates each differing component exactly once. Each changed slot therefore visits only its initial and final records and never returns; unchanged slots simply hold their records. The route meets both predicates simultaneously. This is a positive all-size family with fixed continuously held resources, not a borrowed-resource workspace repair.

For the mode-dependent lifting theorem below, continuous resource occupancy still holds under its fixed-key hypothesis, but the route may revisit a full record while preparing later boundary witnesses. No full-record non-return theorem is claimed there. Nor do these statements address a stronger rule that declares mere continued occupancy at two observation times to be forbidden “reuse”; that is a different event/port semantics and must be specified separately.

## 4. Boundary distance between different equality fibres

The following extension covers supplied equality fibres with DIFFERENT relational graphs, not just one fixed graph. Let F and G be two nonempty solution sets on the same old slots, each given by (2), possibly with pinned roots. Let their components/potentials be C_F(i),p_i and C_G(i),q_i. Their resources and old prefix records are the same in both fibres.

Construct a bipartite **component overlay**: one left vertex for each F component, one right vertex for each G component. Slot i becomes an edge of weight w_i between its two component vertices, with equation

    y_(C_G(i)) - x_(C_F(i)) = p_i-q_i.                  (3)

Retain all pinned roots. Parallel edges are genuine distinct slot occurrences and their weights add. Its weighted frustration index phi is the minimum total weight of edges whose deletion leaves a solvable voltage-and-pin system. “Solvable” requires both zero cycle holonomy and compatible transported pins; cycle checks alone omit pinned-path obstructions.

**Boundary theorem.** The exact minimum record distance is

    min_(a in F,b in G) sum_(i:a_i!=b_i) w_i = phi.      (4)

A minimizing deletion certificate and any solution of its retained equations reconstruct minimizing full records in both fibres.

**Proof.** Given a,b, their component-root values violate (3) precisely at slots where a_i!=b_i. Delete those edges. The remaining system has these roots as a solution, so phi is at most the record distance. Conversely delete a minimum-weight feasible edge set and choose a solution of the remaining equations. Build a,b by the two potential formulas. They can disagree only on deleted edges, so their distance is at most phi. The inequalities give equality. With strictly positive weights every minimum deletion set can be chosen to coincide with the mismatch set. QED.

Therefore F intersects G iff the overlay itself is solvable. Every individual old slot can project onto all three phases in BOTH fibres while phi>0: the obstruction is in cycles of their combined relations, not in their unary phase sets. This is why the binary forced-interval boundary formula cannot simply be replaced by a per-slot ternary coverage check.

No efficient general boundary optimization claim is made. Exhaustive root assignments or edge deletions is always a finite exact algorithm. Any original voltage edge u->v of voltage d can be subdivided through a fresh right vertex e with constraints z_e-x_u=0 and z_e-x_v=-d: for fixed original roots the two-edge gadget has minimum disagreement 0 when the original equation holds and 1 otherwise. Thus arbitrary voltage disagreement objectives occur inside the bipartite boundary problem; bipartiteness alone does not remove the frustration issue.

## 5. Conditional all-size helper/fibre lifting

Supply finitely many authenticated complete helper modes h. For each mode, supply its actual old equality/pin predicates, fixed actual resources, and complete additive charge guard. Helpers have an exact distance d(h,g), counted separately from old slot occurrences. When predicates depend jointly on a helper and TWO old records, this is explicitly a higher-arity constraint; it is not covered by pretending to have only pairwise old/helper interactions.

For each mode, solve its equality system. Split any nonneutral free component and any free component of weight>B into separate visible root-phase choices. Check the full charge equation for each resulting state alpha. Reject impossible states. The remaining F_alpha is a product of neutral components of weight at most B, so it is B-connected. Its full phase choices, actual resources, mode-conditioned predicates and charge data remain as reconstructible sidecars.

Make a quotient vertex for each nonempty F_alpha together with its helper mode h_alpha. Compute delta(alpha,beta) by (4), including pinned visible roots. Join distinct quotient vertices precisely when

    d(h_alpha,h_beta)+delta(alpha,beta) <= B.            (5)

**Exact lifting theorem.** Two accepted full assemblies lie in the same B-budget atomic packet component iff their quotient vertices lie in the same component of this finite graph.

**Proof.** Any full packet has helper cost d plus old distance at least delta, so its projection is an edge or stays in one fibre. Conversely route inside F_alpha to the minimizing boundary witness, using one light neutral-component rotation per differing component. Apply the helper replacement and all boundary mismatches in one atomic packet; (5) is its exact cost and both arrivals are accepted. Repeat, then route within the final fibre to the specified endpoint. QED.

This theorem is exact for the supplied packet graph. It neither manufactures helper records nor gives a sequence of microscopic edges through forbidden intermediate assemblies. Resource injection is authenticated at every arrival; in this lane it follows from fixed distinct prefix keys. If resources vary across modes, they require an additional full capacity check at each supplied fibre and arrival and cannot be erased by this phase theorem.

## 6. Actual native controls that distinguish the obstructions

### 6.1 A frustrated triangle with locally total compatibility

Use S_2, q=16, and three slots with prefix keys 0,1,2. Impose

    owner_0(16)=owner_1(16),
    owner_1(16)=owner_2(16),
    owner_2(17)=owner_0(16).

These are k_1=k_0, k_2=k_1 and k_0=k_2+1. Every individual edge relation has three compatible two-slot pairs and full projection onto each slot's three phases. With the third slot unconstrained, each one-edge system has nine states. Every record has the same target class; every phase uses the same one unit at its slot; all three resource keys are distinct. Yet the three-edge system has no accepted assembly: its voltage sum is 1. Local charge coverage and edgewise matching do not imply a global section.

### 6.2 Two connected neutral fibres with a nonzero frustration boundary

For four old slots, define

    F: (k0,k1,k2,k3)=(a,a,b,b),
    G: (k0,k1,k2,k3)=(c,d,c,d+1).

Each has nine states, every slot has all three phases, and each component consists of two slots. Their overlay is one four-cycle with nonzero voltage. There is no common state, but changing one slot suffices to cross from an F record to a G record. The exact old boundary cost is one slot, or three owner-clock occurrences.

For a sharp helper bridge take TWO disjoint copies of this four-slot system, giving eight old slots. Use S_4, q=32, with old prefix keys 0,...,7 and a helper at prefix key 8. Its supplied local menu contains the actual full records R_8(0),R_8(1). Conditional on helper phase0 impose F in each copy; conditional on phase1 impose G in each copy. The conditional equalities are the declared complete higher-arity predicate. In literal clock form, equality k_j=k_i uses clocks32,32; k_j=k_i+1 uses clock33 at i and clock32 at j.

Each mode fibre has 9^2=81 actual old states and is connected by two-slot packets. Every old phase remains locally available in both modes. Their exact boundary distance is two slots, one from each frustrated overlay cycle. Any helper switch costs one additional slot. The full assembly has 162 states:

- Budget2: exactly two components of size81.
- Budget3: one component of size162.
- Sharp bridge cost: three changed slots, or nine owner-clock occurrences.

This is not the earlier all-equal large-component example: every neutral movable old component has size2, and each mode fibre is connected at budget2. Failure is entirely at the relational boundary between two nonempty, connected, charge-neutral fibres. Repeating the four-slot block L times gives boundary L, helper bridge L+1, and individual mode fibres still connected at budget2. For L>=2, ignoring boundary cost creates a false budget2 route. This is an all-size native conditional counterexample with 4L+1 distinct normalized resources, realizable whenever 2^r>=4L+1.

### 6.3 Resource-efficient parallel-overlay control

The same bridge obstruction needs fewer slots if a parallel-edge overlay is allowed. For each j=0,...,L-1 use one old pair (2j,2j+1). Mode0 requires k_(2j+1)=k_(2j); mode1 requires k_(2j+1)=k_(2j)+1. Their overlay contains two incompatible parallel equations per pair. Each mode has 3^L states, each free component has size2, and the exact boundary distance is L. A genuine helper switch therefore has sharp cost L+1 slots, or 3L+3 owner clocks, with just 2L+1 fixed distinct resources. Realize this in S_r whenever 2^r>=2L+1, retaining the same explicit helper-conditioned predicate qualification.

For L=2, S_3 has order27 and five resources. The full assembly has eighteen states: budget2 gives two components of size9, and budget3 joins all eighteen. The four-slot-per-block construction above is retained because its overlay is a simple cycle and the two mode partitions differ, not because four old slots are minimal for a positive boundary. Two old slots already give the smallest parallel-edge boundary example with full unary phase coverage and no pins.

The role of the helper is not inferred from a target class or a general wiring theorem. The mode-dependent equality predicate is supplied explicitly, and all target partners and normalized capacities remain independently authenticated.

## 7. Reproduction, loss ledger, and open scope

Run `python -B check.py`, `python -O -B check.py`, and `python -B verify_manifest.py` in this folder. The receipt records its exact finite universes. The compiler is checked against enumeration of literal path-owner equality, not its own equations. Boundary values are checked both by direct pairwise Hamming enumeration and independently by deleting overlay edges and solving their remaining constraints. Native tests reconstruct source matrices, complete pinned menus, source rotations, partner maps, normalized reservoirs, converse clock reflection, and exact clock costs. No Python assert is used, so optimization retains all checks.

Map: actual indexed owner comparisons -> voltage graph -> component root phases. Preserved: exact accepted assignments and packet cost, with fixed resources and full record sidecars. Forgotten: concrete path labels and root offsets unless those sidecars are retained. Simple-cycle rejection control: the three-slot voltage triangle. Boundary map: two affine equality fibres -> weighted component overlay. Preserved: exact minimum full-record disagreement. Cheapest lost nuance if only unary phase coverage is retained: the two-slot nonflat parallel-edge overlay; the four-slot control has a simple cycle and genuinely different component partitions.

Open and deliberately unclaimed: arbitrary phase-relation catalogues beyond equality/pins; a native realization of arbitrary edgewise affine reflections; compatibility of extra target-side attachments without marked transports; existence of charged helpers for a global electrical atlas; any global observer or unrestricted Hamiltonian-path routing theorem. This work resolves the relational phase carrier and its conditional neutral boundary, not those larger obligations.
