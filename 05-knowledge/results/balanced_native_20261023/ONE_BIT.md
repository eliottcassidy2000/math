# One-bit completion by resident displacement

**PROVED (scoped).** One auxiliary bit suffices for injective assignments with arbitrary connected local port graphs. This closes the OPEN one-bit question in `guarded_coverage_20261023/PORTS.md`; that earlier packet remains a frozen historical snapshot. It does not remove an electrical guard that is already disconnected after projection.

Concurrent reconciliation: commit `a793424e20`, `native_packet_decomposition_20261023/workspace/ONE_BIT_WORKSPACE.md`, independently proves this mechanism and the stronger version where full-record-to-resource maps have nontrivial fibres. That broader theorem, with its actual C3-kernel example, is the current retained-kernel interface. The present proof and BFS are an independent restricted derivation and check; no separate priority or stronger generality is claimed.

## Statement

Let V be a finite actual port set, and let s labelled slots have connected undirected graphs G_i on their eligible sets E_i subset V. Let a and b be injective eligible assignments. Choose a simple local path gamma_i from a_i to b_i, of length d_i. Enlarge the port inventory to V times {0,1}. Slot i may traverse an edge of G_i at a fixed layer, or change layer at its current base port. The auxiliary bit is unconstrained at intermediate times. No additional changing cross-slot equality is imposed.

Then (a,0) and (b,0) are joined by a sequence of legal one-slot moves with no repeated actual port in any state, of length

    ell <= 2 sum_i d_i + s.

Every horizontal gamma_i is used exactly once, in its given order. Each resident layer change is charged to a particular horizontal event; the final at-most-s layer resets have their own integer clocks. The bound is a guarantee, not an assertion of shortest length.

## Constructive proof

Call a slot finished after it has completed its local route. Finished slots stay at their distinct target base ports, although they may change layer. An unfinished slot not currently travelling remains at its original base port, again with either layer.

Start an arbitrary unfinished slot i along gamma_i. Before it enters the next base port v in layer l, inspect its resident. Away from the traveller's current base port there is at most one resident. If that resident occupies layer l, flip the resident to the vacant other layer; then move i horizontally into (v,l). If the resident is already in the other layer, or v is empty, move directly. The only possible double occupancy of a *base* port is now at the traveller's position. Actual ports (base,layer) remain distinct.

A simple route never returns to its starting vertex. Thus the possible double at the starting vertex, including the start of a displacement chain, cannot obstruct a later step. When i reaches b_i, mark it finished. If another slot j is at that base port, j is unfinished: a finished j would have b_j=b_i, contrary to target injectivity. Immediately make j the traveller. It has never moved horizontally, so it is still at a_j and can use its originally supplied gamma_j. Its first horizontal event removes the sole double at the old base port. Continue this displacement chain.

If the target base port is empty, the chain ends and all base occupancies are again single. Start any remaining unfinished slot. Each iteration finishes a new slot, so there are exactly s completions. A zero-length route is harmless: it can finish only at its unique assigned target, and cannot strand another finished slot. At termination the base assignment is b. Flip each remaining layer-one slot to layer zero; the distinct base targets make these resets legal.

There are exactly L=sum d_i horizontal moves, at most L resident flips, and at most s final flips. This proves the bound and all intermediate injectivity claims. No probabilistic or asymptotic assertion is involved.

## Native tournament specialization

Use the authenticated source T_r and complete pinned cube X_r from `guarded_coverage_20261023/CUBES.md`: r ordinal blocks of the rigid eight-owner source with canonical code 29558785. A cube bit selects one of the two paths

    (3,0,6,7,1,5,4,2), (3,1,6,7,0,5,4,2).

The same six original owner clocks in every block stay pinned. A face menu fixes additional cube coordinates. Append one genuine eight-owner block, obtaining the fixed enlarged source T_(r+1). Its last cube bit supplies the auxiliary layer. For s slots with injective face-menu endpoints, take local geodesics; d_i<=r gives

    ell <= s(2r+1),     changed-owner-clock work <= 2s(2r+1).

Every event changes exactly two actual owner clocks in one block. All source/converse/inverse normalized word orbits are distinct by the inherited cube theorem, so different layers are real capacity, not alternative kernel names for the same resource. The source is enlarged once before the route and then held fixed. Original fixed pins survive, endpoints return to auxiliary bit zero, full chosen paths remain explicit, and converse reverses the corresponding actual clocks. These rigid sources have trivial automorphism lifts; this is not an assertion that arbitrary sources have trivial kernels.

The partner-class constancy inherited from CUBES.md preserves every class-based guard along this route. The supplied block transporter does not fix every interior target owner; arbitrary external target attachments require their own compatibility sidecar. Likewise this theorem does not imply equality of the electrical atlases at orders 8r and 8r+8. It is an authenticated route on the enlarged fixed source, with that source's declared guards.

## Why this does not bypass the primitive electrical circuit

Session45's `incoming-integration/WORKSPACE_BOUNDARY.md` proves that if legal enlarged states project into a guarded base relation A, and each event changes at most b projected choices, they can connect only within the b-move components of A. For A={0^m,1^m}, b<m cannot suffice.

Here the projection deliberately permits two slots at one *base* port, while retaining distinct actual enlarged ports. The projected injectivity predicate is temporarily relaxed, with its exact repair stored by the extra bit and travelling-resident identity. Electrical charge is a different predicate: a neutral bit does not repair a nonzero physical-row charge. One-bit completion therefore solves the occupancy obstruction and leaves the authenticated dyadic circuit obstruction intact.

The earlier four-layer construction costs L+4s events; this construction costs at most 2L+s and uses half the enlarged port inventory. Neither bound dominates the other for all L. Memory and integer work are separate resources.

## Boundary-clock corollary

If a prescribed return kernel R on boundary assignments is stationary for the actual source law pi, compile each sampled endpoint pair by the construction. Observing only completed returns gives exactly R and therefore retains pi. If all routes have length at most a common integer H, add legal holding events at the endpoint to reach H; the marked boundary chain then has a fixed clock. This is not a uniform stationary law on intermediate states, and is not a congestion or spectral-gap estimate. Variable-duration sampling at arbitrary physical times would require its own occupation weights.

## Checks

`one_bit.py` exhausts all 663 feasible ordered tuples of cube-face menus in dimension two with one through three slots, and all 7,912 ordered endpoint pairs. It also tests every multiset of two connected local graphs on subsets of four ports: 64 local graphs, 2,076 feasible menus, 174,222 ordered endpoint pairs. Slot permutation makes the graph-multiset reduction exact for this two-slot test. It adds 128 seeded face cases and fully occupied opposite-corner assignments through dimension eight. Every assertion is an always-active exception test, including per-event legality and injectivity. `audit.py` separately checks the generated native event records and actual owner paths. The proof above supplies all sizes.
