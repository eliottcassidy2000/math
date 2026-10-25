# Exact continuation states for source-owned clock capacity

2026-10-23, Session33 Niche lane. **PROVED, abstract finite theorem; FINITE-EXACT controls.** No physical LRC row, new native annulus census, uniform 153N inequality, arbitrary-source Collatz entry, or historical novelty is claimed.

## Inheritance, delta, and concept board

Read-only repository input: e5b42c67e5d2, clean at startup. AGENTS, startup packet, relevant frontier/guardrails, protocol and core-paper routing were read. Message processor was unavailable because this checkout has no `.machine-id`; none was invented. No Git mutation or publication was performed.

Closest mechanisms:

* `01-canon/theorems/THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile.md`: weighted Hall inequalities after deleting target-diagonal edges, actual ancestry versus marginal labels, minimal two-root swap.
* `05-knowledge/results/annulus_residual_flow_20261020k/README.md`: source-owned per-horizon capacity, legal handoffs, integral interval circulation, exact finite two-cut optimum.
* `05-knowledge/results/residual_flow_fibre_clocks_20261020k.md`: moving external credit out of the remainder, all actual cuts, and common-representative lifting failure.
* `05-knowledge/results/finite_horizon_budget_20261021a/README.md`: actual finite-arrival clock count and the invalid exchange of fixed-depth / all-depth quantifiers.
* `05-knowledge/results/tournament_future_boundary_20261020/README.md`: quotient histories by their responses to every still-permitted continuation.

The old Hall theorem is explicitly inherited, not claimed anew. The new synthesis is an exact characterization of the universally sufficient continuation state, its invertible mask coordinates, its realizability test, and the source-preserving lifting boundary. A bounded targeted search found no prior exact mask/Möbius statement; this is not a claim to historical priority.

Board: actual clock cut / negative source-pair owner / allowed obligation neighborhood / universal demand continuation / shared gate / source receipt. Hostiles: owner-word swap, arrival gaps, same singleton data with different pair capacity, shared gates with negative mask coefficient. Least-used sidecar: the exact allowed neighborhood of each original credit at each actual cut.

## 1. Model and the distinction between simultaneous reuse and later reuse

Fix a finite obligation label set Q and a finite set of actual event cells t. All source activation, expiration, eligibility and demand changes must occur on retained cell boundaries. Empty cells may be discarded; a final unbounded cell is allowed when all data stay constant thereafter.

At cut t, source-owned credits i have nonnegative capacities w_i(t), original immutable source-pair addresses, and eligibility sets A_i(t) subset Q. Zero capacity and empty eligibility never contribute usable credit. A demand vector d(t) in R_+^Q is paid by nonnegative x_ij(t) satisfying

    sum_i x_ij(t) = d_j(t),
    sum_j x_ij(t) <= w_i(t),
    x_ij(t) = 0 unless j in A_i(t).

There are no additional edge capacities, cumulative lifetime budgets, minimum holding times, indivisibility beyond integral units, or shared routing bottlenecks. Handoffs across cells are permitted. These are load-bearing hypotheses, not conveniences.

One source address is never spent twice at the SAME cut. It may pay different obligations at disjoint cuts because its negative indicator term exists at both cuts. This is the inherited annulus accounting convention. A one-shot consumed coupon has a different state space, and this theorem does not model it.

## 2. Exact finite service and compressed Hall obstruction

Suppress t. For S subset Q set

    r(S) = sum_{i: A_i intersects S} w_i.

Then all demand is paid iff d(S) <= r(S) for every S. The exact minimum unpaid mass is

    delta(d) = max_{S subset Q} [d(S)-r(S)].                 (1)

The empty set makes delta nonnegative. For integer capacities and demands, an integral payment achieving (1) exists.

Proof. Build the bipartite network source -> obligations (capacity d_j), obligation -> eligible owner (capacity infinity), owner -> sink (capacity w_i). A finite cut keeping exactly S of the obligations on the source side has minimum capacity d(Q\S)+r(S). Minimizing over S gives max service d(Q)-delta. Integral augmenting paths give the integral claim. This is the inherited weighted-Hall mechanism, now applied to a different explicitly typed service object.

For a fixed bank, let N(S) be its positive-capacity owners adjacent to S and define

    cl(S) = {j in Q: N({j}) subset N(S)}.

Then S subset cl(S), r(cl(S))=r(S), and d(cl(S))>=d(S). Thus it suffices to test closed sets cl(S)=S in (1). In particular, labels with identical owner neighborhoods enter an extremal witness together. This is an exact reduction, not a guaranteed polynomial-size test family.

## 3. The complete minimal continuation observable

Aggregate owners by their exact nonempty eligibility masks:

    c_A = sum_{i:A_i=A} w_i,             empty != A subset Q.

Then

    r(S) = sum_{A:A intersects S} c_A.                      (2)

Two banks have exactly the same feasible demand vectors, for EVERY nonnegative demand vector on the same labeled Q, iff their c_A agree for every nonempty A. Empty-mask credit is unobservable and deliberately omitted.

Proof of sufficiency: (2) and (1). For necessity, the feasible-demand region recovers r(S) as

    r(S) = max {d(S): d feasible and d(Q\S)=0}.              (3)

Indeed assign each owner touching S entirely to any one neighbor in S; this reaches r(S). Thus identical feasible regions give identical r. Let M=r(Q) and

    f(T) = M-r(Q\T) = sum_{empty != A subset T} c_A.

Boolean Möbius inversion uniquely recovers

    c_A = sum_{T subset A} (-1)^(|A|-|T|) [M-r(Q\T)].      (4)

This proves necessity.

For integral banks, the same equivalence holds with only integral demand continuations: the construction proving (3) is integral.

Scope of minimality: universal arbitrary demands are essential. For a fixed single native demand or bounded 0/1 demand alphabet, fewer distinctions can suffice; (3) can ask for a concentration larger than the native obligation permits. This theorem does not call its profile minimal for that smaller observer.

The state size is at most 2^|Q|-1 masks and often only the masks that actually occur. This can compress many source owners, but no unconditional complexity improvement is asserted.

## 4. Exact realizability diagnostic for a proposed compressed rank

Let r be any real-valued set function on Q with r(empty)=0. Compute c_A by (4). Then r is the rank of SOME independent nonnegative-capacity eligibility bank iff every c_A>=0.

Proof. Necessity is Section 3. Conversely, form one owner of capacity c_A with eligibility A for every positive coefficient. Inversion gives f(T)=sum_{A subset T} c_A. Taking T=Q gives sum_A c_A=M because r(empty)=0. Hence its rank on S equals M-f(Q\S)=r(S). Rational coefficients give rational owners; nonnegative integral coefficients give integral owners.

Thus ordinary monotonicity and submodularity of r are insufficient. Example on Q={1,2,3}:

    r(S)=min(|S|,2).

It is realized by a shared gate of total capacity 2 followed by three obligation branches each of capacity 1. But f(singleton)=0, f(pair)=1 and f(Q)=2, so

    c_{12}=c_{13}=c_{23}=1, c_{123}=2-3=-1.

No independent owner-mask bank has the same response to arbitrary demand vectors. The missing state is the shared gate (and branch constraints), not another phase label. A reachability graph alone cannot certify a capacitated flow network. One must either retain those bottlenecks or prove they are absent.

This is a falsifiable interface test: a negative coefficient rules out the proposed independent-credit interpretation even when the rank table is an entirely valid network capacity function.

## 5. Minimal hostile to total and singleton statistics

Let Q={1,2,3}, all capacities 1. Bank A has masks {1,2} and {3}; bank B has masks {1,3} and {2}. Both have usable total 2, and each singleton has rank 1. Yet for demand (1,1,0), bank A pays only 1 while bank B pays 2. The failed implication is “total capacity plus all individual service capacities determine simultaneous service.” The missing coordinate is a multi-label neighborhood overlap.

Three labels are minimal for a hostile preserving all singleton ranks AND total rank: on one or two labels these values already specify every nonempty subset rank, hence all feasible demands by (1). This minimality is about the label count, not every other notion of witness size.

Even simpler, with two labels and two credits, assigning both eligibility masks {1} versus one {1} and one {2} refutes total capacity alone.

## 6. Time gluing and lifting the compressed solution

The full independent-handoff continuation observable is the cell-indexed sequence c_A(t). Equality of these profiles is equivalent to equality of feasible demand histories under arbitrary cellwise demands. Necessity follows by concentrating a demand history on one cell and applying Section 3. Sufficiency follows by solving independently in each cell.

The joint distribution of a donor's entire future eligibility trace is NOT needed for this feasibility observable. One owner active on both cells and two separate owners active on one cell each have the same one-label capacity at each cell. Both pay the same histories with handoff. They differ if a fixed single owner is required or credit is consumed permanently. The continuation language determines the correct quotient.

To lift a coarse allocation, retain as a sidecar, for every cell and mask, the list of original owners and their capacities. A coarse mask node distributes total mass at most c_A among labels in A. Since every owner in its list is eligible for every label in A, split that allocation through the owners successively, never exceeding a capacity. Integral inputs can be split integrally. Repeat independently on each cell and concatenate adjacent pieces with the same owner and obligation. This yields a finite source-addressed receipt with no same-cut recycling.

This is a sufficient physical realizability criterion only when every retained owner is an authenticated original ledger term and eligibility really licenses the proposed scalar/semantic use. The theorem cannot construct absent later LRC target roles or turn root-label coincidences into common ancestry. It cannot discard cross-interactions not included in the paid bank.

## 7. Exact connection to the native 153N remainder

The inherited scalar annulus payment has unrestricted compatibility: every active negative interaction may pay any active positive interaction as a SCALAR bound, without physically transferring charge between terminals. Therefore only c_Q=b(t) is present and r(S)=b(t) for every nonempty S. Equation (1) collapses exactly to max(p(t)-b(t),0). Our theorem does not improve that already exact scalar bound.

Its use is in a stricter future certificate that demands actual ancestry, owner type, terminal, phase or other eligibility. Such restrictions can only increase the matching defect; they do not create negative credit. The exact mask profile says which additional state is necessary to assess the certificate and gives a concrete deficient set if it fails.

If the selected ledger bank is J=p-b in normalized equal units, any partial payment with unpaid mass delta gives J<=delta after unused negative slack is dropped. Restoring physical unit 3 gives I<=R+3delta. For unrestricted compatibility this is the existing bound R+3max(p-b,0); stricter eligibility gives a potentially weaker bound but a stronger source/semantic certificate. Moving external negative pairs into the bank requires subtracting them from the remainder exactly once, as in the inherited R' accounting. No payment quotient changes Delta E=A+I or A+M=153N.

For LRC, the actual later selector and common ancestry remain prerequisites. THM-2545 already uses the opposite operation, deleting the desired diagonal and looking for a Hall obstruction that forces a hit. The shared mechanism is preserving the support of legal continuations; one cannot identify an annulus payment edge with an LRC target edge merely because both are bipartite.

For the tournament continuation work, the map is from an owner bank to its response to all permitted future demands, just as a partial cycle history is mapped to its responses to future insertions. Preserved predicate: feasible continuation. Destroyed information: individual owner identity and cross-time trace correlations. Sidecar: original owner lists and actual clock cells. Cheapest decisive tests: a multi-label Hall cut and the signed coefficients (4). This is a genuine common quotient principle, not a numerical coincidence or a transferred prize theorem.

## 8. Exact artifact and stopping boundary

Run `python3 clocks/check.py` and `python3 -O clocks/check.py` from the workspace parent. Both use explicit checks, not Python assert. The checker compares an independent unit-owner augmenting matcher with every Hall subset over all 128 binary mask profiles on three labels and 8,192 demand cases (each demand coordinate 0..3), verifies inversion and closed-set compression, and checks both hostile controls. The receipt records the exact check count.

The finite universe is synthetic; it is NOT a native Collatz or LRC census. Larger real banks can use the theorem once actual eligibility and timing data are extracted. Extracting those data, proving a useful all-scale restriction, and establishing the uniform residual estimate remain separate open tasks. No theorem ID was reserved and no canonical status was changed.
