# Arithmetic square-sum contexts separate actual boundary pairings

Status: **PROVED, elementary synthesis of two inherited mechanisms**;
same-session audit and exact positive/hostile controls. No priority claim.

The abstract separator is already proved in
`../session33_continuation_states_20261023/ports/PORT_CONNECTION_ALGEBRA.md`,
Theorem 2. The fixed-label arithmetic bridge is already proved in
`../induced_square_forests_subdivisions_20261019/INDUCED_SQUARE_FORESTS.md`
and its independently reviewed `SUBDIVISION_MOD7_CONSTRUCTION.md`.
The new connection here is that the separator can be implemented by
**native positive-integer square-sum contexts while retaining the source
labels**. We do not rediscover either abstract theorem or assert interval
Hamiltonian coverage from a sparse embedding.

## The actual-state theorem

Let S be a finite set of distinct positive integers. Let F and F' be two
spanning linear forests in its induced square-sum graph, each consisting of
r>=2 nontrivial paths. Assume their degree-one vertex sets are the SAME
2r named ports B, and all other vertices have degree two. Their paths induce
perfect matchings P and Q on B. Suppose P!=Q.

There is a finite set of new positive integer labels, disjoint from S, such
that, in the induced square-sum graph on the combined label set:

1. F extends to a Hamiltonian cycle using r private length-three connecting
   paths, with all source vertices and every edge of F preserved.
2. F' cannot extend to a Hamiltonian cycle while preserving every edge of F'.
3. Each added vertex has degree two in the full induced graph. There are no
   new-to-old or new-to-new edges other than those on the designated paths.

The total number of added vertices is exactly **2r**. The labels may be large.
There is no claim that the combined set is an integer interval, or that a
universal bounded-height context exists.

### Step 1: a separating matching

Choose ab in Q but not P. The P-paths containing a and b are distinct.
Orient and order all P-paths so the first starts at a and the last ends at b.
Join successive path ends by r-1 abstract connector edges R. Then F union R
is one spanning path with endpoints a,b. Put N=R union {ab}; N is a perfect
matching on the ports and F union N is one spanning cycle.

For F', its a-to-b path together with the connector ab forms a closed
component. At least one other F' component remains, since r>=2. Hence
F' union N is a disconnected 2-factor. This is the inherited abstract
separator, with the last edge restored to make a cycle test.

### Step 2: realize every connector without changing the old labels

For a requested connector A--B and a current label bound H, choose r0 to be
the first integer greater than 2H+4 with parity A+B+1. Set

    v = r0^2-B,
    s = (A+B+1-r0^2)/2,
    w = s^2-A.

Then A--w--v--B is a square-sum path, since

    A+w=s^2,   w+v=(s-1)^2,   v+B=r0^2.

Here s is an integer, H<v<w, and the new labels are distinct. For any old
y!=B, v+y=r0^2+(y-B), with nonzero deviation smaller than the adjacent-square
gap 2r0-1. It is therefore nonsquare. Likewise for y!=A, w+y=s^2+(y-A)
is nonsquare since 0<|y-A|<H<2|s|-1. These statements also check that w has
no unintended edge to B and v none to A. The one new-new pair is the intended
middle edge. All old labels and their full induced adjacency remain unchanged.

Repeat this construction for the r connectors of N, using the updated H
each time. Later steps avoid every earlier internal vertex, so every added
vertex has final degree exactly two. It is legitimate for A and B already
to be adjacent in the old graph: the new length-three path is a distinct
route, not a duplicate edge in a simple graph.

### Step 3: no alternate route can repair the rejected state

A Hamiltonian cycle in the enlarged graph must use both incident edges at
every new degree-two vertex. Thus it must use all connector paths. Retaining
F or F' then saturates every old vertex at degree two. No unused old edge can
be inserted. The only candidate extending the selected state is its union
with the connector paths. Suppressing their degree-two interiors recovers
the connected cycle for P or the disconnected 2-factor for Q. This proves
all three assertions.

For marked Hamiltonian PATH extension one can instead omit the final ab
connector and prescribe endpoints a,b; it uses 2(r-1) new vertices. The
endpoint marking is necessary: unmarked path queries allow new endpoint
choices and are a different predicate.

## What compression is still possible?

For actual states on the same S and B, different endpoint pairings cannot be
identified by any deterministic code that must answer all such fixed-source
square-sum cycle-extension tests. The continuation family is now arithmetic,
not merely an unrestricted abstract graph gadget.

This does NOT assert that all (2r-1)!! pairings occur as spanning states on
every S, nor that this lower bound applies to one restricted Gerbicz catalog.
It also does not forbid signed linear compression of weighted sums. The
Session33 distinction between distinct response rows and matrix rank remains
essential. Exact counts, parity counts, and existence are different observers.

## Explicit hostile and tests

Use S={1,3,35,46}, whose induced square graph is a four-cycle. Take

    F  = {{1,3},{35,46}},
    F' = {{1,35},{3,46}}.

They have identical named ports, degrees, component counts and no isolated
vertices. Apply the native connector construction for N=F'. The first state
closes to an eight-vertex Hamiltonian cycle; the second closes to two
four-cycles. Exact labels and all induced edges are retained in
`native_contexts.json`. This uses the square-sum four-cycle whose sharp
birth at 46 was proved in `../crossroads_family_20260926_squares.md`.

`verify_native_contexts.py` checks every pair of labels for unintended edges,
not only the six intended new square equations. It also tests all relabelled
and reversed versions of the two connector choices and square-scaled copies.
The implementation uses explicit exceptions and passes normally and with -O.
The all-size theorem is the proof above, not the finite test range.

## Bridge to the current tournament target

The same general research move is useful in the order-ten proof: retain
distinguished-component multiplicities until the complement is accounted
for. The two arguments do not identify square-sum graphs with tournaments.
There the target is a weighted determinant-gap inequality; here it is
extension of a selected spanning state. In both cases prematurely replacing
joint incidence by marginal counts loses a term that later operations see.

The next useful square-sum question is the opposite one: for the ACTUAL
finite nice-pair transition catalogs, which of these native contexts remain
legal? Compute response equivalence under that restricted context family.
Any smaller state space must be justified by its restricted operations, not
by an invariant that the universal separator has already disproved.
