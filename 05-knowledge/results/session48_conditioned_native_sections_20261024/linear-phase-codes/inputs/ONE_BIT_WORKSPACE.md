# One genuine resource bit routes arbitrary full-record menus

23 October 2026. **PROVED, scoped.** The construction below resolves the
one-bit completion conjecture stated in
`guarded_coverage_20261023/SYNTHESIS.md`, and extends it to local menus with
retained kernel fibres. It does not remove dynamic shared-owner or electrical
constraints. **FINITE-EXACT** controls are recorded by `check.py`.

## 1. Inheritance and the distinction that matters

The incoming `guarded_coverage_20261023/PORTS.md` supplies a four-layer route
and an actual native cube in `CUBES.md`; its one-extra-bit conjecture remains
open in that snapshot. Session45's
`incoming-integration/WORKSPACE_BOUNDARY.md` proves that workspace cannot
remove a barrier preserved by projection. Its
`injective-completion/THEOREM.md` distinguishes raw ports, normalized unit
tiles and full kernel decorations. Our earlier
`injective_native_blocks_20261023/assembly/BINARY_ASSEMBLY.md` handles
two-state dynamic compatibility by exact implication blocks. None of these
statements identifies a kernel lift with a fresh normalized resource.

The new mechanism is to distinguish one **moving slot** from passive slots.
The passive slots always have different base resources; the moving slot may
share a base resource with one of them in the other workspace layer. Passing
a passive slot requires switching its layer, not changing its base record.
When the mover reaches its destination, the displaced slot becomes the mover.

Concept board: actual resource / full record / passive injection / displacement
chain / shared-owner projection. Anchor: one-bit routing; niche: noninjective
kernel-to-resource maps; wildcard: the sharp dynamic-owner counterexample.

## 2. Full-record theorem, including nontrivial kernel fibres

For each labelled slot i let G_i be an undirected graph of authenticated full
records. A record includes its actual owners, path, target, full kernel label,
and every required local guard. Let q_i map records to one common set V of
physical resources. The maps q_i need not be injective: different full lifts
can consume the same normalized tile. Resource equality is literal equality
in V, after the specified physical-owner normalization.

Choose endpoint records a_i,b_i with pairwise distinct q_i(a_i), and likewise
pairwise distinct q_i(b_i). Supply a path gamma_i in G_i from a_i to b_i,
of length ell_i; write L=sum ell_i and s for the number of slots.

Assume a **genuine two-layer stabilization**:

1. Every record x has two authenticated versions (x,0),(x,1), consuming
   resources (q_i(x),0),(q_i(x),1). Distinct pairs are distinct resources;
   records with the same q and layer still compete for one unit.
2. Every local edge x--y is legal in each layer. Every vertical edge
   (x,0)--(x,1) is legal. Full records and kernel labels are transported as
   specified, rather than replaced by arbitrary representatives.
3. Apart from resource injection, all requirements hold for these local
   records and edges independently. Fixed consistent shared owners may be
   retained. An unresolved condition on changing owners in different slots,
   a joint kernel condition, or an aggregate charge equation is not covered
   merely by listing it as an additional annotation.
4. The immutable background occupies none of the resources used by the
   supplied routes in either layer. Equivalently, any background restrictions
   have already been removed while preserving assumptions 1--2.

**Theorem.** The endpoint assemblies ((a_i,0)) and ((b_i,0)) are joined by
one-slot legal moves. Each gamma_i is traversed exactly once in its given
order. There are L horizontal moves, at most L passive vertical moves, and
at most s final vertical resets. Thus the route uses at most 2L+s events.
All full records, actual resource identities and kernel fibres are retained.

No shortest-route claim is made. A horizontal move with q(x)=q(y) is a real
full-record change and is charged once; it needs no passive vertical move.
For nonunit event costs, retain the supplied horizontal edge costs and the
actual vertical costs of the particular slots switched. The event bound does
not identify those costs or discard integer clocks.

### Proof and constructive algorithm

Begin with all slots passive and choose any unfinished slot i as moving.
At this moment every passive slot occupies a different base resource. Traverse
gamma_i. For an edge x--y with different base resources, at most one passive
slot j is at q_i(y). If j occupies the mover's layer, switch j vertically.
Its other layer is empty: passive base resources are distinct and the mover
is at q_i(x), a different base. Move i horizontally. If j was already in
the other layer, simply move i. If q_i(x)=q_i(y), the same-layer resource
remains occupied by i and no other slot, so apply that horizontal edge
directly. The passive base injection survives every step.

On reaching b_i, mark i finished. If no passive slot is at q_i(b_i), the
current displacement chain ends; all base resources are now distinct again.
Otherwise exactly one passive slot j is there in the opposite layer.
Declare j moving and i passive. Passive injection is restored because j was
the only other slot at that base. Slot j is unfinished: a finished slot there
would have target resource q_j(b_j)=q_i(b_i), contradicting endpoint injection.
The full record of j is still a_j, because it has only undergone vertical
switches so far. Its supplied gamma_j is therefore ready to traverse.

Each handoff finishes a previously unfinished slot. No completed slot can be
displaced at a target, although completed passive slots may switch layers
when another route passes them. Consequently a displacement chain terminates
after finitely many handoffs. Continue from another unfinished slot if any.
After all slots finish, their base resources are the distinct target keys;
reset every layer1 slot vertically to0. Each horizontal edge triggers at
most one passive switch, and each slot needs at most one final reset. This
proves legality, reconstruction, termination and the claimed count. QED.

The partial map i->j defined by q_i(b_i)=q_j(a_j) has indegree and outdegree
at most one. Its path/cycle decomposition describes the handoffs, but is not
needed as a separate computation. A dependency cycle ends at the first
mover's vacated base. There is no requirement that intermediate base
assignments remain injective: that is exactly the obstruction the extra
resource bit removes. The full enlarged assignment is injective throughout.

The literal reference implementation scans passive slots and costs
O(sL+s^2) resource comparisons after the local paths are supplied. Maintaining
the passive base-to-slot dictionary instead gives O(L+s) dictionary
operations, plus full-record transport/guard-verification costs. Neither
bound includes discovering the local record graphs or their paths.

## 3. Actual tournament realization: one block, rather than two

Use the rigid eight-owner source of `guarded_coverage_20261023/CUBES.md`,
canonical code29558785, with

    P0=(3,0,6,7,1,5,4,2), P1=(3,1,6,7,0,5,4,2).

The partner isomorphism g=(4,0,2,3,1,5,6,7) fixes the first and last owners.
Ordinal composition gives source T_r of order8r. Pin the six common clocks
in every block; its complete pinned fibre is the cube {0,1}^r. Every cube
edge changes precisely the two indexed owners at clocks8j+1 and8j+4.
All paths have isomorphic partners. The source is rigid, and the full
source/converse/inverse word reservoir contains 4*2^r distinct words.
Thus these cube vertices are distinct normalized resources, not hidden
automorphism lifts. These facts are inherited; the checker independently
reconstructs the actual matrices and isomorphisms on its stated finite range.

For every slot choose a pin-defined face E_i, and injective endpoints a,b.
Append ONE more rigid block and release its two possible paths. This is the
literal product E_i x {0,1} in T_(r+1), of order8r+8. A shortest coordinate
route inside each face has ell_i=Hamming(a_i,b_i), so the theorem gives

    events <= 2*sum_i Hamming(a_i,b_i)+s <= s(2r+1),
    actual owner-clock work = 2*events.

Each event changes one slot. Every target class is preserved at every event,
so every target-class-dependent per-slot guard is preserved, and the selected
class-incidence/current difference is identically zero. This is stronger than
one modular check but does NOT authenticate arbitrary owner-dependent guards
or the native electrical atlas of an unrelated source extension. The actual
first/last owners, common pinned arcs and converse transports are retained.
On converse paths clock j becomes8(r+1)-1-j. All auxiliary bits return to0.

This proves the incoming one-bit conjecture in full. It also applies to any
connected subgraphs of this cube for which actual local routes are supplied;
the menus need not be faces. It changes the source by adjoining eight owners.
It does not find an unused buffer inside an unchanged native atlas.

One auxiliary bit is best possible uniformly in this product model: two
slots with the same one-edge menu {0,1}, initially (0,1) and finally (1,0),
have no legal one-slot move without an extra layer. This is already realized
by the seed's two pinned paths. The assertion is about auxiliary resource
bits, not a minimum number of added tournament owners among every conceivable
construction.

The theorem's noninjective q_i version allows full kernel fibres in the local
records, provided a genuine two-layer native realization is separately
authenticated. It does not assert that appending a block automatically
stabilizes an arbitrary target-changing or shared-kernel family. In particular
the incoming carry coordinates may be retained, but their many raw lifts
are not multiple normalized resources. The concrete rigid cube used here has
trivial source automorphism kernel.

### A concrete native kernel, with three full lifts per unit

There is also an all-size native example with a nontrivial kernel. Let
S_r=T_r -> C3, r>=1, where the final owners8r,8r+1,8r+2 form the directed
cycle. Its full records are

    R(b,k)=P_b followed by (8r+k,8r+(k+1 mod3),8r+(k+2 mod3)),
    b in {0,1}^r, k in Z/3.

Pin the same six clocks in each rigid prefix block, and leave the three final
owners unpinned. The ordinal decomposition forces each Hamiltonian path to
visit its blocks consecutively; the C3 has exactly these three Hamiltonian
paths. Thus the complete pinned family has3*2^r full records. All source
automorphisms preserve each strong component, fix the rigid prefix and rotate
the final C3, so Aut(S_r)=C3. In particular q(R(b,k))=b gives exactly THREE
full owner lifts per normalized source word, not three independent resources.

The same statement holds for the full normalized K/converse tile. Distinct
b cannot have equal source words, since that would give a source automorphism
changing a rigid prefix choice. S_r is not isomorphic to its converse: their
ordered strong-component profiles are (8,...,8,3) and (3,8,...,8).
Every partner is strong. Indeed its first eight-owner block is the inherited
strong partner; every owner in a later block points to the first owner of
that first block (none of these arcs is an exceptional boundary-path edge).
The retained global Hamiltonian path then reaches every owner from that first
owner. Thus a source word cannot equal any inverse-source word. Combining
these separations yields4*2^r different normalized words, or2^r unit tiles,
above3*2^r full records.

A cube-coordinate move uses the endpoint-fixed g in that prefix block,
extended by the identity elsewhere; it changes two actual owner clocks.
A kernel move k->k+1 uses the literal automorphism rotating the final C3,
extended by the identity on the prefix. Equivariance of K supplies the actual
partner isomorphism, and exactly THREE owner clocks change. The reverse
rotation is equally valid. Hence every record has the same target class,
while its actual final owner values and kernel coordinate are retained.
The prefix pins remain fixed. These final C3 owners are not simultaneously
declared fixed shared pins; extra conditions on them must restrict the local
routes explicitly.

Stabilize by inserting one rigid block immediately BEFORE C3:
S_r of order8r+3 becomes S_(r+1) of order8r+11. The owner embedding fixes
the original prefix and sends the old final owner8r+j to8r+8+j; this embedding
is part of the full record, not a silent reuse of local owner names. The new
block supplies the real bit f; full records are(b,f,k) with resource(b,f).
All three k-values still consume one unit at each(b,f). Any supplied path in
the full-record graph (cube edges and C3 rotations), with injective endpoint
base keys, therefore satisfies the one-bit theorem. A route has actual
owner-clock work2*events plus the number of horizontal C3 rotations, rather
than2*events. This is an actual retained-kernel realization, with all its
resource fibres explicitly distinguished.

## 4. Sharp boundary: a changing shared owner cannot be parked away

Take two slots in T_2, with cube coordinates (b0,b1), and face menus

    E0={00,01}, E1={10,11}  (binary strings write b1 b0).

Require that their actual owners at clock1 agree. In the first seed block
that owner is0 for b0=0 and1 for b0=1, so this is literally b0(slot0)=b0(slot1).
The complete accepted assembly consists of

    (00,10) and (01,11).

Both are injective; all four possible base paths use different normalized
resources. Both preserve every target-class guard. They differ in two slots
and four owner-clock occurrences. Append any number of neutral workspace
blocks, allowing their bits to vary independently, but retain the clock1
equality. Every legal state still projects to one of the two displayed
assemblies. One-slot moves cannot change that projection, so no one-slot
route exists, regardless of the number of auxiliary ports. A synchronous
two-slot base-bit flip succeeds, proving the sharp budget2 (four changed
owner clocks). Minimality here is in the number of slots: a genuine
cross-slot equality needs at least two. It is not a least-tournament-order
claim. Shared owner names do not assert a common-host amalgam; in this
example the two records do use the same actual source matrix.

This is a native instance of Session45's projection obstruction. The new
one-bit theorem repairs capacity, while shared changing values need actual
synchronous events or a separately justified modification of the constraint.
The distinction survives every modulus because it is equality of owners.

## 5. Map, loss ledger and reproducibility

Source: full authenticated local records plus their physical resource map.
Target: the enlarged injective assembly graph. Map: retain each full record
and add one real resource-layer bit. Preserved: slot identity, route order,
local guards, exact kernel labels, actual resource capacity and integer
event clocks. Deliberately relaxed after forgetting the bit: base injection.
Required sidecar: unused native layer resources and authenticated vertical
edges. Cheapest failure if omitted: the incoming tight Gray cycle, whose
base assembly has only two states and no one-slot route. Cheapest changing
owner failure: the two-slot example above.

Run `python -B check.py` and `python -O -B check.py` in this folder. Both compare
against the frozen receipt; `--write` deliberately regenerates it. Text-source
hashes normalize CRLF to LF for portable clones. The
script independently verifies each emitted state, tests arbitrary record
graphs including noninjective kernel maps, reconstructs native cubes, and
checks the shared-owner hostile. Its JSON states the finite universe; the
all-size claim rests on the proof above. No Collatz, LRC or unrestricted
tournament wiring conclusion is asserted.
