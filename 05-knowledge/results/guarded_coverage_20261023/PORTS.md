# Heterogeneous ports: exponential tight cycles and constant-size workspace

23 October 2026. **PROVED:** native Gray-cycle obstruction and four-layer
stabilization. **PROVED / classical specialization:** the tight-Hall component
criterion. **OPEN:** arbitrary cross-slot equalities, full-atlas rewiring,
or a uniform weighted mixing bound for all assemblies.

## 1. Native menus with exponential necessary simultaneous support

Use the complete pinned cube X_r from [CUBES.md](CUBES.md), r>=2, and a cyclic
binary Gray ordering g_0,...,g_(m-1) of its m=2^r paths. There is one labelled
assignment slot for each Gray edge. Slot i has precisely the eligible paths

    E_i={g_i,g_(i+1)}, with indices modulo m.

These are **actual pin-defined menus**: retain the six common clocks in each
core; on every bit except the differing Gray coordinate, pin clock8j+1 to
its required owner0 or1 (offset by8j). Exactly the two displayed full paths
survive. All have the same actual source and partner class; hence every
target-dependent per-slot guard agrees. Impose unit capacity on each actual
path port, and retain all labelled slots. Shared constant pins are compatible.

The incidence graph is a chordless2m-cycle. It has exactly two perfect
matchings: all slots use their first endpoint, or all use their second.
Choosing either endpoint for one slot forces the choices around the cycle.
They differ in every one of the m slots. Consequently the sharp simultaneous
changed-slot budget is

    2^r = 2^(n/8), where n=8r is the source tournament order.

This extends Session44's native three-slot obstruction and the older
[graceful matching-support mechanism](../dynamics_labeling_entry_20261019/GRACEFUL_RELEASE_AND_MATCHING_SUPPORT.md).
The abstract chordless-cycle obstruction is already known there; the new
point is its pin-exact native realization with exponentially many distinct
ports in one fixed source. The incoming
[binary assembly theorem](../injective_native_blocks_20261023/assembly/BINARY_ASSEMBLY.md)
also identifies the same indivisible block through its literal implication
component. Its native family uses many distinct sources; this family uses
one source and many injectively owned paths. Neither result is a claim that
all unrestricted guarded fibres have just these two states.

**Failure anatomy.** Connectivity of each owner's path cube does not imply
connectivity of an injective multi-slot assembly. The missing coordinate is
vacancy. All eligible ports in this cycle are occupied. No single slot can
take its alternative without consuming another slot's actual resource.
An unrelated unused port outside all E_i does not repair this obstruction.

## 2. A one-block repair of the cycle

Append one genuine eight-owner seed block, so the fixed source during the
repair is T_(r+1). Keep its six common pins; leave its new bit unpinned.
Slot i now has the actual four-port menu E_i x {0,1}. Start and end in layer0.
All layer1 ports used below must be unused by the unchanged background.

To turn the first cycle matching into the second:

1. Move slot0 from(g_0,0) to(g_0,1).
2. In order m-1,m-2,...,1 move slot i to(g_(i+1),0).
3. Move slot0 to(g_1,1), then to(g_1,0).

Every step is one cube edge, changes one slot, changes two actual owner
clocks, and uses an unoccupied eligible port. The total is m+2 events and
2(m+2) changed slot-clock occurrences. All auxiliary bits return to zero.
This is an authenticated stabilization into a larger source, not a claim
that the frozen original source can move without new resources.

## 3. Four-layer stabilization for every connected local menu

The more general theorem does not need a tournament.

**Theorem.** Let V be an actual port set. Each labelled slot i has a connected
graph G_i on its eligible set E_i subset V; an edge is an authenticated legal
local move. Let a,b be two injective assignments with a_i,b_i in E_i. Enlarge
the actual port inventory to V x {00,01,11,10}, a four-cycle of layers, and
permit horizontal G_i moves and vertical layer edges. With no other changing
cross-slot constraints, there is an injective one-slot-at-a-time route from
(a,00) to(b,00) of exactly

    sum_i length(gamma_i) + 4s

events, for any supplied simple local routes gamma_i from a_i to b_i.

**Proof.** Process slots in any fixed order. Move the current slot from
(a_i,00) to(a_i,01), follow gamma_i in layer01, then move to(b_i,11).
Layer01 contains only that slot; unfinished slots remain in00; finished
slots occupy distinct b_j in11. There are no collisions. Once every slot
is in11, layer00 is empty. Move each slot from(b_i,11) through(b_i,10) to
(b_i,00). Distinct target ports ensure legality. All slots remain assigned
throughout, and every layer edge is a genuine resource change. QED.

**Native corollary.** For cube X_r and any fixed pin-defined subcube menu at
each slot, two additional seed blocks realize these four layers in source
T_(r+2). The local shortest route length is the Hamming distance d(a_i,b_i),
so every pair of injective assignments has a route of

    sum_i d(a_i,b_i)+4s <= s(r+4)

native events, changing exactly twice this many owner-clock occurrences.
The source has8r+16 owners, independently of the number of slots. Its actual
endpoints, all stipulated fixed pins, partner target class and converse
transports survive. The auxiliary layer bits are marked at the endpoints
and explicitly released during the return. Arbitrary connected, non-face
menus also work if authenticated local routes are supplied; their lengths
replace the Hamming bound.

The statement preserves every source slot and capacity-one resource. It
does **not** preserve arbitrary extra equalities between *changing* owners
of different slots. Such equalities can force a simultaneous block, exactly
as in the incoming binary assembly theorem. The four-layer reservoir must
be actually available, including inverse and converse ports. This is a
selected-inventory theorem, not permission to duplicate resources already
owned by the full atlas. A source measure on endpoint assignments can be
carried unchanged as marked demand; no distribution on temporary states or
polynomial spectral gap is inferred without a separate conductance bound.

## 4. Exact obstruction without stabilization: tight Hall sets

For a finite feasible eligibility graph L -> R, write Gamma(S) for the
eligible ports of S subset L. A one-slot vacancy move replaces its assigned
port by an unused eligible port, considered here as one legal endpoint move.
Let T be the union of all sets S with |Gamma(S)|=|S|.

**Theorem.** Two assignments are connected by vacancy moves iff their actual
assignments agree on T. Every connected component has a route diameter at
most2s^2+2s, where s=|L|. T can be found from any assignment by alternating
reachability from the unused ports: it is exactly the unreachable slots.

**Proof.** Hall's inequality holds because a saturating assignment exists.
Submodularity of |Gamma(S)| and Hall's lower bound show that the union and
intersection of two tight sets are tight. Hence T is tight. All its ports
are always occupied by slots in T. A first vacancy move changing any such
slot is impossible, so its entire restriction is invariant.

Delete T and Gamma(T). For any nonempty remaining U, its remaining
neighborhood has at least |U|+1 ports; otherwise T union U would be tight.
Equivalently, from any remaining slot a reverse alternating route reaches
a free port. If it did not, its closed reachable slot set would have as
neighbors exactly its matched ports, a forbidden tight set.

To connect two remaining assignments, process their alternating difference
paths by moving the free port along each path. The remaining differences
are alternating cycles. For a cycle of k slots, choose a shortest alternating
path from a free port to the cycle; its t internal slots lie outside the
cycle and are distinct. Move the vacancy through those t slots, move the
first cycle slot out, rotate the other k-1 slots, move that first slot to
its target port, and restore the external path. This costs2t+k+1<=2s+1
and changes no already repaired external slot at the end. At most s cycles
and at most s total path-slot moves give the displayed safe bound.

Finally start vacancy reachability at every free port, traversing an eligible
port-to-slot edge and then that slot's current matched port. No tight set is
reachable. Conversely the unreachable slots' neighborhood consists exactly
of their distinct matched ports, making that set tight. It equals T and
is independent of the starting assignment. QED.

This is a port-language specialization of classical matching reconfiguration,
not a novelty claim for the abstract criterion. In
[Bousquet–Hatanaka–Ito–Mühlenthaler, arXiv:1812.05419v1](https://arxiv.org/html/1812.05419v1),
Theorem18 and AppendixD use the frozen perfectly matched part of the
Edmonds–Gallai decomposition; AppendixA identifies that part by alternating
reachability. Their general framework also distinguishes reachability from
the harder shortest-route problem. Here all slots are saturated, so every
token jump is exactly a vacancy move at the same labelled slot.

If a native slot can only move along some edges of its eligible graph,
this abstract endpoint theorem alone does not authenticate a native route:
its intermediate ports can be occupied. The four-layer theorem supplies
that missing realization under its explicit reservoir hypothesis.

## 5. Verification and next questions

The cycle checker exhausts complete assignments for2<=r<=8 and checks each
event of the buffer repair against literal native clock pins. The four-layer
checker exhausts all663 feasible tuples of one to three two-dimensional
face menus and all7,912 ordered endpoint pairs, plus144 seeded higher-
dimensional pairs and fully occupied cubes through r=8. It records event
counts, not merely an abstract connectivity flag.

The Hall checker exhausts all eligibility graphs at dimensions1x2,2x2,
2x3,3x3,3x4. It compares subset enumeration of tight sets, independent vacancy
BFS, and complete assignment-graph components, including disconnected spare
ports. These tests support the proofs; they do not replace all-size arguments.

The next strongest questions are: can one extra native bit suffice for all
connected face menus; can a reservoir be authenticated in the *unchanged*
ambient source; and can simultaneous equality blocks admit a comparable
stabilization without losing their shared-owner constraints? No answer is
assumed from the successful two-bit construction.
