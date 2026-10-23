# Adaptive original electrical observers: four clocks can beat six

23 October 2026. **PROVED**: the finite adaptive-observer Bellman recurrences and their exact interpretation. **FINITE-EXACT**: optimal original-character, merged-target and literal-path observers on four complete original order-eight source fibres. No all-size wiring, network reconstruction from incomplete ports, or native spectral-gap theorem is asserted.

## 1. Actual inherited guard and the new coordinate question

The native electrical graph has converse-merged tournament classes as vertices. Its literal edge network is `G=B+2Q`; loops retain two stubs but contribute zero to the Laplacian. For an actual Hamiltonian path `P` of `T`, the partner `K_P(T)` keeps its seven path arcs and reverses every other arc. The inherited primitive potential `x` is a grounded mod256 harmonic character of this network, with root3527. Its original mod8 port phase is

    phi_x(T,P) = x_[T] - x_[K_P(T)] modulo8.

These are actual merged owner coordinates, not synthetic clock residues or the residue of the Hamiltonian count. Complete source hashes and original gauge are pinned in SOURCE_PINS.json at commit `0a36c32af8c89f6c4249a0ef929390432e556925`.

The inherited `native_guard_moves_farey_20261023/guards/README.md` proves that this phase is constant on every source automorphism orbit and unchanged by converse plus actual path reversal. Its static indexed-clock footprint is six on each of the four sources below. Here a query of clock `j` reveals the actual owner `P[j]`. The next clock can depend on prior answers. This changes the observation protocol, not the underlying path, source, electrical character or owner map.

**New finite conclusion.** The exact optimal worst-case adaptive costs are:

| Merged source | Ordinary class | Actual paths | Frozen phase | Merged target | Literal path |
|---:|---:|---:|---:|---:|---:|
|1238|1763|409|5|5|5|
|2242|3481|407|4|4|4|
|2246|3490|523|4|4|4|
|2315|3676|527|4|5|5|

Every nonadaptive cost in these twelve cases is exactly six; the independent checker recomputes this, rather than importing the earlier footprint result. In particular source2315 exhibits strict adaptive output compression: one original three-bit electrical character needs four queries, while the actual merged target needs five. At sources2242 and2246 even the whole literal path can be recovered with four adaptive clock queries.

This is finite nonuniform observation with a known source and a precomputed decision policy. Building the policy already uses the original target-owner atlas. It does not derive the phase or target from an unknown network.

## 2. Exact Bellman theorem and proof

Let `S` be a nonempty finite collection of literal paths of a fixed source and `F:S->Y` an observable. A clock query `j` partitions `S` into its nonempty answer classes

    S_(j,v) = {P in S : P[j]=v}.

Queries that are constant on the current state are useless and omitted. Define `D(S)=0` when `F` is constant on `S`. Otherwise

    D(S) = 1 + min_j max_(v:S_(j,v) nonempty) D(S_(j,v)).

**Theorem.** `D(S)` is exactly the least possible worst-case number of adaptive clock queries needed to determine `F(P)` for every `P in S`.

**Proof.** Every correct nonterminal decision tree chooses a first nonconstant clock. After answer `v`, its remaining tree must determine `F` on precisely `S_(j,v)`, so its depth is at least the displayed recurrence. Conversely choose a minimizing clock and attach optimal subtrees for each answer. Every child is a proper subset, so induction on `|S|` terminates. No unqueried label, target owner or phase may be used to choose the next query. The supplied trees make each decision only from the current answer class. QED.

A repeated query is constant because its answer has already been recorded. Thus this recursion cannot gain a fictitious advantage by querying a position twice. A leaf can terminate before identifying the path, but only if every actual path left has the same required output.

The minimax theorem does not assume that the source is rigid. All source automorphism lifts must be included in `S`, unless the observation problem explicitly quotients them. In this finite packet all four source automorphism groups are independently checked to be trivial, so no hidden kernel lifts have been discarded.

## 3. Average cost is a different objective

Fix the explicitly chosen measure: uniform over all actual Hamiltonian paths of the source. For each state define the minimum total external path length `E(S)`, not its average. The terminal value is zero, and otherwise

    E(S) = |S| + min_j sum_v E(S_(j,v)).

The proof is the same first-query induction: a query costs one for each path in `S`, and the remaining total cost is the sum over disjoint children. The optimal expectation at the root is `E(S)/|S|`.

To minimize average cost subject to a worst-case bound `b`, define `C(S,b)=0` at a homogeneous state; a nonhomogeneous state has `C(S,0)=infinity`. Otherwise

    C(S,b) = |S| + min_j sum_v C(S_(j,v),b-1),

where any choice with an infinite child is forbidden. This is not obtained by simply taking independently mean-optimal children: they can exceed the depth allowance.

For the original phase the exact results are:

| Source | Minimax depth | Best mean with no depth cap | Best mean at minimax depth |
|---:|---:|---:|---:|
|1238|5|1478/409|1478/409|
|2242|4|1423/407|1424/407|
|2246|4|1891/523|1905/523|
|2315|4|1900/527|1914/527|

At the last three sources every mean-optimal policy has worst-case depth at least five, and a depth-five mean-optimal policy is supplied by the recurrence. Hence the tiny one-query-saving tradeoff is real, not rounding: source2242 trades one unit of total cost over407 paths for reducing the worst case from five to four. The same receipt gives both objectives for merged target and full path.

## 4. Explicit policies, lower certificates and exact validation

receipt.json contains an explicit clock/answer tree for each of the twelve objectives. The displayed tree minimizes expected cost among all trees of minimum worst depth. Each terminal stores the common output and actual path mass. Every tree is replayed on all its paths, and its maximum depth and external path length are checked.

The lower bound is not an information-count heuristic. All reachable partial-observation states are solved by exact recursion. The state counts, shared by the three output objectives for these fixtures, are6188,5724,8411,8471. An independent tuple-of-path-indices recursion repeats the computation without importing the producer's bitset solver. It separately declares budget `D-1` infeasible by the depth-constrained recurrence. The two solvers agree exactly, including rational means.

All1,866 paths are enumerated by both permutation filtering and a distinct depth-first extension algorithm. For each record the checker constructs the literal partner arc by arc, checks its Hamilton path, and finds an actual isomorphism to the original canonical target tournament. The saved target isomorphism is verified on every pair. This gives a separate semantic check on the original forward-word owner lookup rather than trusting a residue alone.

Converse is checked literally: reverse `P`, transpose `T`, normalize again in the original chord gauge, and confirm the same source/target merged owners and phase. It transports clock `j` to `7-j`. The same optimal observer costs therefore apply to the transported converse problem, with every clock reflected. They do not authorize reusing unreflected physical pins.

These optimal clock trees are not hidden ternary trees. Their internal answer arities range from two to eight, and every published root has eight answers. A ternary cell symmetry in another representation cannot change that native owner-query alphabet by nomenclature alone.

## 5. Concurrent native-move recovery, not a second novelty claim

While this packet was computing, commit `a2015157620` published the same four-source native move census in `native_guard_moves_farey_20261023/native_moves/NATIVE_GUARDED_MOVES.md`. The reconfiguration part here is an independent recovery/audit of that incoming result. The adaptive-observer computation above is the new delta.

A move in this paragraph is a simultaneous replacement of one actual Hamiltonian path by another, in the same source, preserving the two indexed endpoint owners and this one frozen electrical phase. Its budget is the number of changed actual owner clocks. It is not a number of cell phases, an elementary metagraph edge, or a sequence of accepted one-clock changes.

For every complete endpoint/phase fibre, a minimum-bottleneck spanning tree gives the least budget making that fibre connected. The checker independently verifies connectivity at that threshold and disconnectedness one below. The original six-clock-barrier counts agree with the incoming census:20,19,35,31 for the four sources respectively. The complete histogram is retained.

A particularly small two-state certificate uses source1238, actual pins `P[0]=0`, `P[7]=3`, and original phase2. There are exactly two survivors:

    P = (0,4,2,1,6,7,5,3), target2289, forward mask1903700;
    Q = (0,6,7,2,1,5,4,3), target2009, forward mask1081112.

All six interior clocks change. Every budget at most five leaves the two vertices isolated; one budget-six atomic move connects them. Before imposing the phase, there is a budget-four actual route

    (0,4,2,1,6,7,5,3)   phase2, target2289
    (0,4,2,6,7,1,5,3)   phase0, target2501
    (0,6,7,4,2,1,5,3)   phase7, target2684
    (0,6,7,2,1,5,4,3)   phase2, target2009.

The rejected intermediate phases exhibit the precise effect of hard guarding. Treating this route as a legal guarded trace would silently allow forbidden intermediate states; the incoming auxiliary return-time theorem does not make that legal.

**One frozen character is not the full electrical bank.** This move does not claim to preserve every mod8 harmonic character. The independent electrical-duality lane supplies the separator `y=4 e_2289`: all edges incident2289 have even conductance, hence `L y=0 mod8`, but its two port phases at source1238 are4 and0. An aggregate or selected-character preservation statement must not be upgraded to per-port universal-bank preservation.

If the separately audited universal-bank result identifies its output partition with the target-owner partition on a fixed native source, the target column of section1 immediately gives the universal-bank adaptive costs too. This implication uses equality of partitions, not an extra decision-tree experiment; it remains dependent on that separate electrical theorem.

## 6. Scope and reproduction

The original atlas and frozen character are inherited inputs, not recomputed discoveries. Our actual input electrical receipt SHA256 is `417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27`. The prior guard receipt records a different whole-file checksum (`d8ef442...`) despite matching the character values consumed here. SOURCE_PINS.json binds these computations to the actual read bytes; no conclusion depends on treating the old checksum as current.

Commands, from this packet:

    python3 -B check.py
    python3 -O -B check.py
    python3 -B independent.py
    python3 -O -B independent.py

For a published layout without frozen input copies, pass `--repo PATH_TO_REPO` to check.py. It locates the three original files by their repository-relative paths and verifies their full byte hashes. Both normal/optimized replay pairs agree byte-for-byte. No floating arithmetic, optimization-sensitive assertion, external solver, or unpublished source is needed.

Open: bounded adaptive observers for other sources or growing order; a network-independent way to compute the electrical phase; uniform guards under arbitrary shared-port attachment; a complete native reconfiguration theorem; and all-size wiring. None follows from a small observed decision depth, the recovered six-clock move, an auxiliary coupling gap, or the presence of ternary cell symmetries elsewhere.
