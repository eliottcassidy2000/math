# Simultaneous charge sections need disjoint full-record allocations

24 October 2026. **PROVED, conditional all-size:** cancellation-free occupancy compiler, the one-resource Hall specialization, sufficient hyperedge depletion criteria, and exact paid boundary lifting. **PROVED, scoped native all-size compiler:** supplied resource-set menus on the inherited rigid-prefix/C3 carrier, with retained kernel lifts, converse, class-dependent guards and literal costs. **FINITE-EXACT:** the small hypergraph censuses and actual native obstructions in `check.py`. This does not construct all-order charged native supply, transport the order-eight numerical dyadic bank, or give a global observer.

## 1. From a repaired hole to simultaneous source-pair demands

The inherited single-slot rule in `inputs/RECEIVER_BANKS.txt` says that a compatible binary repair of label lambda repairs a requested grade r precisely when lambda avoids H_r-H_(r-1), uniformly over the retained full boundary. That is a section-existence statement for one requested debt. Several demanded sections can each exist yet compete for the same normalized resources. Three kernel lifts of one tile remain one capacity unit.

Fix an actual complete separator tuple b, a finite set D of labelled demands, and a common normalized resource set U. Each demand i owns a disjoint set of local slots. A witness w in W_i(b) is the complete assignment of those slots: source matrices, actual paths and target partners, all kernel lifts, charge g_i(w) in a specified finite abelian group Q, exact integer grade vector a_i(w), every required indexed pin, and the entire normalized K/converse resource closure. Let R(w) subset U be its resources, after normalization. Reject any witness whose distinct owned slots already repeat a resource. Different lifts using the same resource set are different witnesses in the same capacity fibre, not additional units.

Immutable background resources are unavailable and are removed first. Shared boundary resources are either charged to one declared owner once, or are unavailable to every interior witness. They cannot be repeatedly charged as separate copies. The same once-only ownership rule applies to charge and clock occurrences. Shared owner labels need not be shared normalized resources; per-slot owner-clock costs need not equal costs counted on distinct global owners.

Every remaining factor is retained. For a general relation Phi on several demands, an allocation is a tuple (w_i) such that Phi(b,(w_i)) holds, all demanded charges/grades hold, and the R(w_i) are pairwise disjoint. Local nonemptiness alone does not imply this. The following independent product is valid without a final relational filter only if, at the fixed complete b, every noncapacity cross-factor is already enforced by the local menus. Call this the **boundary-separable case**. Arbitrary receiver-wide relations do not automatically satisfy it.

## 2. A nonnegative, commutative square-zero compiler

For each resource u use a commuting variable x_u with x_u^2=0. Precisely, the coefficient semiring consists of finitely supported maps from subsets of U to nonnegative integers, with multiplication

    [S] [T] = [S union T] if S intersection T is empty, and 0 otherwise.

This is N[x_u : u in U]/(x_u^2=0) in the indicated monomial-quotient sense. Its ordinary integer coefficients do not wrap modulo a characteristic. It is NOT exterior algebra: variables commute; there are no alternating signs and no determinant cancellation. Nilpotence enforces capacity, not antisymmetry or a claimed physical boson/fermion equivalence.

Write e_q e_t=e_(q+t) for the group-semiring labels of Q and use commuting grade variables z (Laurent variables if signed integer grades are permitted). Retain a separate formal tag t_(i,w) for every full witness. In the boundary-separable case form

    P = product_(i in D) sum_(w in W_i(b))
            t_(i,w) e_(g_i(w)) z^(a_i(w)) product_(u in R(w)) x_u.       (1)

**Exact occupancy theorem.** A term survives precisely when its selected full witnesses form a disjoint allocation. Its charge and grade are the exact once-owned sums. After extracting desired charge/grade, a nonzero coefficient exists if and only if an allocation exists. Retaining the tags reconstructs every actual witness, not merely its charge or resource support. After setting tags to one, coefficients count labelled full-record allocations, including kernel multiplicity; they do not count additional resource capacity, unmarked hosts or target classes.

Proof: distributivity selects exactly one witness per demand. A repeated resource makes that product zero; otherwise it survives with coefficient one and the claimed charge/grade. All coefficients are nonnegative, so no surviving selection can cancel another. QED.

For prescribed separate debts q_i and grades r_i, first restrict each menu to that exact charge/grade (or give each demand its own charge/grade variables). Total charge alone does not enforce separate debts. Exact grades are never replaced by unions over grades.

For an arbitrary relation Phi, retain each tagged tuple and discard it unless Phi holds. Equivalently, enumerate partial witness tuples, retain the separator data of every factor that is not yet closed, and check each factor when its full scope is present. Only then may histories with identical resource support, charge, grade and future-visible separator data be aggregated. This is an exact finite compiler, not a promised small-width or polynomial-time algorithm. Aggregating equal resource/charge monomials before checking a relation on discarded kernel lifts can manufacture witnesses or lose them. A relation connecting a new demand to an omitted interior coordinate must enlarge the retained state.

The construction is charge-coordinate agnostic. Q must already be authenticated as equivalent to the complete requested physical guard at the retained integer margins. A convenient quotient of a physical bank is not authenticated by (1).

## 3. Where Hall is exact, and where it is false

Assume the boundary-separable case and every surviving witness uses exactly one resource. Let N_i be the set of resources with at least one full witness at the requested charge/grade. Then an allocation exists **if and only if**

    |union_(i in J) N_i| >= |J| for every J subset D.                  (2)

Necessity follows from injection. For sufficiency use a maximum partial matching. If a demand is unmatched and there is no augmenting path, start alternating reachability from the unmatched demand. Its reached demands J have all neighbors reached; each reached resource is matched to a reached demand, with fewer reached resources than demands. This violates (2). Hence a full matching exists, and any retained witness over each matched resource reconstructs the allocation, because b already separates all other factors. This proves the usual Hall mechanism, with the required full-record hypothesis made explicit.

A violating J is a sharp obstruction certificate, not just a count of missing charges. If every N_i is the same N, at most |N| demands can be served even when every charge has arbitrarily many kernel witnesses.

For witnesses consuming several resources, neither neighborhood counts nor matching number of unions replace coloured hypergraph matching. Take two demands

    A = {{0,1}, {2,3}},      B = {{0,2}, {1,3}}.                       (3)

Every option has two resources. Each demand has two pairwise disjoint options. For every subset J of demands, its resource union has at least 2|J| resources; even the uncoloured union hypergraph has matching number at least |J|. Yet all four A/B option pairs intersect, so no disjoint allocation exists. Treating the two slots in each demand independently would also create forbidden mixed pairs. Four resources are minimal for the stated two-disjoint-options property. The native realization in Section 6 retains this exact obstruction.

One-resource Hall can also fail to imply feasibility if nonseparable cross-factors remain: three private keys with phases k_0=k_1, k_1=k_2, k_2=k_0+1 satisfy every resource condition but have no compatible tuple. That is a relational obstruction, not a capacity obstruction.

## 4. Useful sufficient criteria for full resource hyperedges

Remain at fixed b in the boundary-separable case and fixed requested charges/grades. Collapse each menu only for capacity testing to its distinct sets F_i={R(w):w in W_i}; retain the entire witness fibre over every set for relations, exact costs and reconstruction. Let k_i be the largest set size, with k_i=0 for an empty family. An empty family fails the criterion below. Empty resource sets, when genuinely available, consume nothing and are allowed.

For U_0 subset U define h_i(U_0)=|{E in F_i:E intersects U_0}|. For K>=0 let

    B_i(K)=max_(|U_0|<=K) h_i(U_0).

Fix any order of demands, and let K_i be the sum of preceding k_j. If

    |F_i| > B_i(K_i) for every i,                                  (4)

then **every** greedy sequence in that order can be completed: the occupied union before i has size at most K_i, so fewer than |F_i| options meet it. Choose any disjoint option and an actual full witness. Induction proves the claim.

An equivalent statement is that no set of at most K_i resources hits every option in F_i; when there is no empty option, its transversal number tau(F_i)>K_i. If an empty option exists, tau is infinity. This is a sufficient all-prefix condition, not necessary for the existence of one good allocation, and tau need not be easy to compute.

Two elementary, constructive sufficient replacements are:

- Supply more than K_i pairwise disjoint options in F_i. Every occupied resource hits at most one member of this supplied pack.
- Let d_i(u) count options containing u. The sum of the K_i largest d_i(u), clamped to the resource universe, is less than |F_i|. A union bound then proves (4). In particular |F_i|>K_i Delta_i suffices when all d_i(u)<=Delta_i.

For maximum witness size k and m demands, a pack of k(m-1)+1 options for every demand suffices in any order. It is not an optimal universal rainbow theorem. The strict inequality cannot in general be weakened to equality: (3) has k=m=2 and pack size k(m-1)=2, yet fails. This is a sharp boundary example for that elementary criterion, not a claim that the bound is optimal for every k,m.

These tests quantify over **whole normalized resource sets**, not numbers of raw paths or kernel lifts. Capacity loss from banning a resource means deleting every incident full witness. A single resource ban can erase many witnesses and many charges. A label-only hole-difference test supplies no such quantitative depletion estimate by itself.

If the criteria hold for every full boundary b in an alphabet D_0, any nonempty allowed relation Gamma subset D_0 admits an allocation by choosing an allowed b. Conversely robustness for every nonempty Gamma implies it for every singleton. This retains the inherited same-boundary quantifier; it says nothing about a new factor outside that boundary or about small-budget motion between its sections.

### A direct conditional bridge from forbidden labels to simultaneous repair

Suppose the requested receivers i have pairwise disjoint complete base-resource reservations, and are otherwise independent after a complete boundary b is fixed. Each supplied actual binary repair slot p has its two declared full records and a complete reservation bundle B_p for both options. Assume the B_p are pairwise disjoint for different candidates, avoid every base reservation and immutable background outside this candidate pool, and are assigned to at most one receiver. A bundle reserves both alternatives, but an actual binary slot occupies only its one selected normalized unit at any arrival; it does not occupy two units merely because two options are reserved. All slots, charge contributions and exact population margins have unique owners; reassignment only changes bookkeeping, not the number of occupied units. Unassigned slots retain specified frozen full outside records, with their actual population and charge contributions explicitly subtracted from the requested residual problem. An assigned slot is removed from that outside ownership and counted exactly once in its new receiver, not duplicated. Full cross-pins and guards must make all separately selected records glue.

Before applying the rule for a candidate p, subtract its actual U-record charge g_p(U) and baseline grade contribution from the requested receiver debt/grade, in addition to all frozen outside contributions. The difference label is lambda_p=g_p(V)-g_p(U), not g_p(V) unless U is already normalized to zero. The displayed r-1 shift refers to the V-count coordinate with U=0,V=1; retain every other exact margin separately. Uniform coverage over all residual charges makes this baseline translation harmless for the success edge, but reconstruction must use the correct candidate-specific residual. For receiver i let F_i be the union of its exact hole differences H_(b,r)-H_(b,r-1) over the requested boundaries and grades, after these candidate-specific coordinate requirements are applied (use F_(i,p) if the relevant grade domain depends on p). Make an edge i--p precisely when the actual two-record repair is available, all its required gluing hypotheses hold uniformly on that domain, and lambda_p avoids the applicable F_i (or F_(i,p)). Retain the full repair witnesses and source-specific costs, not just the labels.

**Uniform repair-assignment corollary.** There is an injective assignment of these declared successful repair candidates exactly when their bipartite candidate graph satisfies Hall. Under the disjoint-reservation hypotheses, any such assignment gives simultaneous compatible complete receivers, hence sections for every vector of requested debts at the specified grades. Proof: Hall assigns distinct p, the inherited hole-difference theorem supplies each requested full section, and the complete reservations plus the stated cross-factor hypotheses make their union feasible. No new normalized capacity is created.

The iff concerns the declared candidate assignment, not all conceivable multi-slot repairs or instantaneous allocations that might share parts of unused alternative reservations. Reserving both options is deliberately stronger than checking one arrival. If complete repair bundles overlap, they are hyperedges, and Hall on candidate names does not certify physical disjointness. If base-receiver stores overlap, this corollary is inapplicable even when every repair candidate differs.

A common stock of M candidates, with every demand losing at most e actual candidates and M-e at least the number of demands, satisfies Hall. Here e counts rejected physical candidate slots, not forbidden charge labels or raw lifts. The old eight order-eight repair candidates do not by themselves authenticate eight disjoint base receivers; that additional native supply is not supplied by this corollary.

## 5. Exact boundary work needs two occupancy alphabets

Let A_alpha and A_beta be two actual accepted allocation fibres on the same labelled owned slots. At each arrival enforce every charge, grade, relation and resource capacity. For full assignments X,Y use their literal nonnegative occurrence cost

    c(X,Y)=sum_(owned slot s) c_s(X_s,Y_s).

Use **separate** resource variables x_u^- and x_u^+ when enumerating pairs (X,Y). They enforce injection separately before and after the event. Using one alphabet would wrongly kill a legal pair that continuously retains a resource at both endpoints. Keep an additional variable v to exponent c(X,Y); use two cost variables when retaining slot support and clock cost simultaneously. Retain their paired exponents and apply both caps together; separate coordinatewise minima need not be jointly attained. Sum only over permitted endpoint pairs and any extra transition predicate. The minimum occurring exponent is exactly

    delta(alpha,beta)=min_(X in A_alpha,Y in A_beta) c(X,Y),           (5)

with infinity when no allowed pair exists, and coefficient sidecars recover an actual minimizing pair. Cost need not factor by demand until actual shared ownership is assigned; if it does, each demand-local witness pair contributes its own literal factor.

If helper work is owned separately, its actual cost is added, never minimized independently while discarding crossing witnesses. Assume any extra transition predicate is symmetric for the undirected statement that follows. At budget T, contracting A_alpha is exact for connectivity only when A_alpha is internally T-connected. Otherwise contract its T-components separately. Proof: every full event projects to an admissible quotient edge; to lift each quotient edge route inside its source component to the retained minimizing witness, perform that complete atomic event, and continue inside the next component. All interior preparation is real work and must obey T. Independent minima need not use the same middle record. For an asymmetric transition predicate use directed reachability, split fibres into their strongly connected components, and retain directional crossing witnesses; an undirected component contraction is not justified.

This is atomic feasible-arrival semantics, not a sequence of microscopic moves through partially replaced records. Two-alphabet capacity enforces no simultaneous reuse at either arrival. It does NOT enforce a temporal ban on reacquiring a resource after release; that requires route history or a special continuously held-resource construction. Nor does it imply per-slot full-record non-return or global-observer existence.

## 6. Native compiler, sharp hostile, and a paid six-state boundary

### 6.1 The all-size supplied-menu compiler

The historical native seed and endpoint-marked gluing are byte-pinned in `SOURCE_PINS.json`. Let S_r=T_r -> C3, q=8r, with rigid seed code 29558785 and

    P0=(3,0,6,7,1,5,4,2), P1=(3,1,6,7,0,5,4,2).

For a key b in {0,1}^r and phase k in F3 let R(b,k) concatenate the r selected seed paths and (q+k,q+k+1,q+k+2), with suffix arithmetic modulo 3. The common six pins in each prefix block leave exactly these 3*2^r records. Each key b is one normalized K/converse resource with four closed words and exactly three kernel lifts. Distinct keys have disjoint complete four-word closures. The full menus, not only charges, are retained.

To realize any supplied finite resource-set catalogue on at most 2^r keys, give each demand a fixed number of owned slots, order every permitted resource tuple, and retain precisely the corresponding full-record tuples. Variable-size hyperedges require a separate authenticated treatment of unused slots; the present compiler uses fixed size per demand. Arbitrary catalogues are explicit full-record predicates, not claimed to follow from unary owner pins alone. Equality of phases between slots is the literal shared-owner relation owner_i(q)=owner_j(q). A common actual phase k separates these equality factors. This supplies genuine native full records for the allocation problem, including multi-slot hyperedges.

For any two records, direct owner-clock disagreement is exactly

    d(R(b,k),R(c,l))=2 popcount(b xor c)+3 [k != l].                  (6)

Actual slot support is [b!=c or k!=l], not the number of changed key bits. The endpoint-fixed partner map g=(4,0,2,3,1,5,6,7) implements each prefix change; rotation of the suffix implements each phase change. Hence all records have the same target isomorphism class. Every per-slot guard depending only on source and target classes, including any dyadic such row, is constant. This is not an assertion about arbitrary owner-sensitive current rows, an authenticated numerical higher-order electrical atlas, or transfer of the original 82-slot C128 x (C2)^4 bank. Those remain separate native hypotheses.

Converse reverses source and path and sends clock j to q+2-j. It preserves (6), the full resource closure and the shared equality after reflected clocks are used. No raw-owner phase shift is silently dropped: k=owner(q)-q modulo 3 in this chart. Distinct prefix keys do not magically possess one global observer; their forced prefix successors can conflict. The compiler authenticates local records and relational allocation, not a common global observer or arbitrary target-side attachments.

### 6.2 Two demands, no allocation despite all local coverage

Use r=2 (carrier order 19), four owned slots, menus (3), and require the suffix-owner phase to agree across all four slots. Each demand has six full witnesses before fixing the shared phase, or two at each fixed phase. Every witness has two distinct normalized resources; its three phase lifts are one resource pair, not three independent pairs. Any A/B choice repeats a resource, so there is no accepted global allocation. This is a capacity obstruction even after the shared relation is solved. Every class-dependent charge is constant, so each locally requested charge/grade is fully witnessed; no all-charge surjectivity of a nontrivial numerical group is claimed.

The obstruction survives arbitrary additional kernel multiplicity on those same resource sets. It also defeats the stronger-looking uncoloured union-matching test in Section 3. This is a native two-record-demand example, not just a hypergraph assigned hypothetical charge labels. All demands use the same native source class and target class; arbitrary heterogeneous source/target-pair supply is not asserted.

For a same-carrier repair, start with the same hostile key catalogues inside an already authenticated S_3 inventory (keys 0 through 7), reserving the unused keys. Adding one previously unused key can succeed: add {0,4} to B and choose A={2,3}. The sufficient pack threshold need not hold. Alternatively adding {4,5} gives B three disjoint options; its pack size 3>2 guarantees extension after any A choice. These repairs use genuinely available authenticated keys at r>=3; they are not free kernel lifts. Enlarging S_2 itself to S_3 would change the carrier, so keys 0 through 3 across those two carriers are not silently identified as the same physical resources. With A as in (3), B augmented by {4,5}, and one common phase, exactly six full global records survive.

### 6.3 Six states expose the internal-fibre and boundary costs separately

Use S_3 (order 27), four old slots and one helper slot. In helper mode 0 the old key tuple is F=(0,1,4,5); in helper mode 1 it is G=(2,3,4,5). The old four slots have one common suffix phase k. The helper always uses key 6, with phase equal to its mode. Its declared menu contains only phases 0 and 1; the third mode is not silently discarded from an allegedly complete unfiltered menu. The full mode-dependent resource tuple predicate is supplied explicitly.

There are exactly six full accepted assemblies (mode,k). All five occupied resources are distinct in every arrival. Within one mode, changing k changes four old slots and costs 12 owner-clock occurrences. Between modes at the same k, two prefix choices change (cost 4), and the helper phase changes (cost 3): the exact bridge costs **3 slots and 7 clocks**. Between different modes and different k the cost is 5 slots and 19 clocks.

Therefore the clock-budget-7 graph has precisely three components of size 2, while at budget 12 all six states connect. A quotient with only two unsplit mode vertices has a 7-clock edge but falsely identifies the three different old phases inside each mode. The exact component-aware quotient does not. These values are computed from actual owner arrays, partners, normalized resources and retained full records, not just assigned weights.

At one fixed key allocation, the common-phase family has three states, and a direct change uses 4 old slots/12 clocks while continuously holding all four old resources. It neither reacquires any resource nor returns to a departed full record. Mode-varying routes have only per-arrival injection unless temporal history is separately enforced. There is no bounded-locality all-size connectivity theorem hidden in this example.

## 7. Loss ledger and verification

Map: full sections -> coloured resource hyperedges -> square-zero occupancy product. Preserved with witness tags: exact charge, exact grade, full compatibility, kernel lifts, normalized capacity and reconstructible actual arrivals. Lost by untagged resource projection: relational distinctions and literal costs. Lost by treating lifts as resources: capacity. Lost by uncoloured matching: demand identity. Lost by one before/after resource alphabet: legal continuous occupancy. Lost by contracting unconnected fibres: budget components.

`check.py` uses no Python asserts, checks normal and optimized execution, rebuilds native matrices and all selected normalized closures for r=1..4, compares algebra with literal tuple enumeration, exhausts stated small Hall/hypergraph catalogues, and checks the native obstruction and exact boundary graph. `verify_manifest.py` checks every frozen file. `reproduce.py` copies only the manifest-listed payload into a new temporary directory, reruns both interpreters there, and compares mathematical receipts. The resulting actual relocation receipt records its execution rather than claiming that relative paths alone prove portability.

All-size mathematical claims are supported by the proofs above and inherited native composition proof, not by extrapolation from r<=4. Missing obligations remain all-order charged native supply for the original complete current bank, arbitrary owner-sensitive guards, unrestricted all-size relational allocation efficiency, fixed-host global-observer realization, microscopic event decomposition, and temporal non-return for variable-resource borrowing routes.
