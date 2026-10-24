# From local partner hosts to a single global observer

Status: PROVED conditional all-size exact compiler. The incoming local-host
controller does not require this stronger condition and is not contradicted.

## Inheritance and missing predicate

Incoming `compatible_charged_helpers_20261023/controller/PARTNER_GLUING.md`
correctly distinguishes a coherent global partner tournament from a partner
obtained by complementing nonpath arcs of ONE global Hamiltonian observer.
A local source/path record is Hamiltonian on its own induced cell. It need not
be a contiguous segment of a global Hamiltonian path. This note gives an exact
criterion and a concrete failure for the incoming shared-owner example.

Fix a FULL labelled source tournament S on a finite owner set V. Each cell's
actual local path declares every internal pair either required-consecutive
or required-nonconsecutive. Source orientations are already coherent with S.
Retain all such flags; any contradictory flag on the same pair makes the
instance impossible. We seek a Hamiltonian path P of S such that K_P(S)
restricts to every prescribed local partner. This is equivalent to satisfying
all these adjacency flags, not merely preserving the local subsequence orders.

## Exact directed-block compiler

Orient each required-consecutive pair according to S, and let F be their
union. A global observer exists only if:

1. no pair is both required and forbidden;
2. every owner has indegree at most1 and outdegree at most1 in F;
3. F has no directed cycle.

These conditions say F is a vertex-disjoint union of directed paths, including
isolated owners. They are necessary because F must be a subset of one directed
Hamiltonian path. Each directed path component must occur as a contiguous
block, in its uniquely prescribed orientation.

Contract those blocks. Make a directed transition B->C exactly when the last
owner of B points in S to the first owner of C AND that pair is not forbidden.
Then:

    accepted global observers <-> Hamiltonian paths of this block digraph.

The bijection expands each block in its required order. Every global path
respects F so determines such an ordering, and every admitted transition is
precisely a legal unblocked source arc. Internal pairs that are not neighboring
within a block never become adjacent; internal positive edges were tested
against the forbidden list. Thus no negative flag is lost.

This is an exact reduction, not an efficient algorithm for arbitrary instances.
The block digraph need NOT be a tournament: different endpoint pairs may allow
both directions, or neither. Applying tournament Hamiltonian-path existence
without checking this new domain would discard the obstruction.

A useful sufficient condition, when there are no extra absolute-time pins,
is that this block digraph be semicomplete: each distinct block pair has at
least one directed transition. The usual insertion proof gives a Hamiltonian
path in a semicomplete digraph: if a vertex can precede the first, insert it;
otherwise follow the path until the first position it can precede, using
semicompleteness to obtain the preceding incoming arc; if none, append it.
This proof is supplied here so no external graph-class theorem is assumed.

## Integer clocks survive contraction as block lengths

If an additional GLOBAL pin says owner v must occupy clock c (zero based),
and v has offset j within its block, then that block must start at c-j.
All global pins in one block must give the same start; out-of-range starts
are impossible. A block's length is its number of actual owners. In a block
ordering its start is the sum of lengths of preceding blocks. Thus the exact
clock-retaining reduced problem is a Hamiltonian block ordering with these
prescribed cumulative-length values. A semicomplete transition digraph alone
no longer suffices when these additional global clocks are imposed.

Local cell clocks in the incoming packet are NOT automatically global clocks.
Converting them without an explicitly chosen embedding of time coordinates
would add a constraint the native controller never promised.

## The actual shared-owner witness cannot be one global observer

In the incoming example the core176 owners are0..7, and helper368 local owners
embed by rho=(8,9,10,2,6,11,12,13). Its shared-owner equality admits exactly
(core U,helper V) or (core V,helper U).

For (U,V), the core path contains 2->3, while the embedded helper path
contains 2->12. Both would have to be global observer edges, impossible
because a path has only one successor at owner2.

For (V,U), the core contains 6->3 and the helper contains 6->12, giving the
same impossibility at owner6. All other host-completion arcs are irrelevant.

The existing local source and partner host constructions remain valid, as
does their original-order-eight electrical controller. This result only rules
out the stronger claim that their partner host is induced by a single global
observer agreeing with every local consecutive-pair flag. The incoming text
explicitly declines that claim. It identifies a distinct coordinate to retain
if a future all-size wiring theorem asks for global-observer realizability.

## Connections and scope

The earlier linear-forest expansion counts path-compatible arc subsets.
Here the SAME elementary linear-forest condition is a realizability gate for
actual required observer arcs, followed by an exact endpoint-transition graph
and length-clock constraints. This is a genuine operation on named arcs,
not an analogy inferred from matching numbers. Disjoint resources, kernel
lifts, electrical harmonics and converse remain separate requirements.
Converse reverses S, each block and each block ordering, and sends global
clock c to |V|-1-c; it preserves feasibility under those explicit changes.
