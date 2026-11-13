# Consecutive star buffers for graceful forest completion

**PROVED scoped: G1–G3 below.** The abstract tree may change. This is a
completion and exchange theorem, not a proof that every prescribed tree is
graceful. No external priority claim is made.

## G1. Sharp bounded entry

Let F be a forest on the fixed labels 0,...,m. Its edge differences are
distinct and are exactly {1,...,m}\C, where C is nonempty. Write r=|C| and
D=max C. Adjoin labels m+1,...,m+D−1 and freeze all their edges to 0.
Then F extends to a graceful tree on 0,...,m+D−1, retaining every old
label and edge. The r missing differences can be installed in any order,
each by choosing a crossing edge of that difference.

Proof. The enlarged frozen forest still has r+1 components. For a missing
d, consider all difference-d edges together with the frozen forest.
The difference-d graph consists of residue chains modulo d. If d<D,
the D−1 consecutive new labels meet every residue class; their frozen
spokes join every chain to 0. For d=D, the new labels meet every class
except m modulo D. If D=m, that class already contains 0. If D<m,
color m is frozen, and its unique possible edge is (0,m); it joins the
remaining class to 0. Thus every free-color contact graph is connected.

At any intermediate forest with more than one component, every remaining
color therefore has an edge crossing its components. Add such an edge.
No cycle is created, the color is new, and after r additions the graph
is a spanning tree. Its differences are precisely 1,...,m+D−1. QED.

The bound D−1 is sharp for every m=D≥3 under this prescribed star padding.
Take C={m} and

    F={(0,d):1≤d≤m−2} ∪ {(1,m)}.

The only isolated old vertex is m−1. With b<m−1 new star leaves, every
available color-m edge (i,i+m), 0≤i≤b, lies inside the other component.
At b=m−1 the edge (m−1,2m−1) is the unique completion. This is a sharp
bound for this extension operation, not for arbitrary relabeling or padding.

## G2. Exchange diameter after entry

At the G1 buffer size, and at every larger size with the same frozen-star
convention, any two completions can be joined by at most r moves, each
changing at most two free-color edges. This follows from the inherited
[connected color-contact theorem](../graceful_connected_color_contacts_20261019/CONNECTED_COLOR_CONTACTS.md).
All its contact hypotheses were proved in G1, so this is an application,
not a new general rainbow-tree exchange theorem.

The mechanism installs a target edge. If its fundamental cycle contains
the old edge of that color, exchange them. Otherwise remove a non-target
edge of another free color from that cycle, remove the old edge of the
target color, and use connected contacts to replace the other color across
the remaining cut. Previously installed target edges remain untouched.

For the inherited order-12 frozen-color obstruction (13 labels, free
colors 6,9,11), the original five completions have components of sizes 2,3.
Adding just label13 with frozen edge (0,13) gives18 completions in one
component. A literal two-move bridge, specifying starts of the three edges
of colors (6,9,11), is

    (2,0,0) → (2,4,1) → (4,2,1).

The middle state uses edge (4,13) of color9. The endpoints belonged to
different old components. This is an actual repair using a new label.

## G3. Exact eventual stabilization

Let the maximum label n grow, always freezing the new spokes (0,x), x>m.
For every n≥m+D the completion sets, with their one/two-color adjacency,
are literally identified by deleting the additional forced leaves.

Proof. A free edge joining two new vertices would form a triangle with
their two frozen spokes, hence cannot occur. Every free edge has an old
endpoint ≤m and difference ≤D, so its larger endpoint is ≤m+D.
All subsequent labels are therefore forced leaves. QED.

The threshold m+D is uniformly sharp. For C={m}, freeze
F={(i,m):1≤i≤m−1}. At n=m+b, 0≤b≤m, all b+1 possible color-m edges
connect the two frozen components, so there are exactly b+1 completions.
The last new completion appears at n=2m=m+D. Thereafter there are m+1.

This produces a useful class of infinite objects with finite tails: the
entire later exchange metagraph is a fixed finite graph with forced leaves
attached to every state. The preserved data are labels, frozen edges and
free differences; prescribed abstract tree type is the destroyed data.

## Exact controls

`graceful_padding.py` checks every621 partial-color forest through m=5,
the two sharp families for3≤m≤30, and the inherited obstruction through
maximum label24. Its completion counts for n=12,...,24 are

    5,18,26,34,40,50,63,70,82,82,89,89,89.

All n≥13 in this particular example have connected two-color fibers.
The general theorem guarantees connectivity by n=22 and stabilization by
n=23; an example may stabilize earlier. Enumeration checks support the
proof and do not supply its all-order quantifier.
