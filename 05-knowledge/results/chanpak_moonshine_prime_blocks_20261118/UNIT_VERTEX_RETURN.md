# A one-vertex addition cannot supply a renewable +2 step

**All-order, stronger than the transitive-clone restriction.**
If a tournament T and its vertex deletion T−v have Hamilton counts
differing by exactly two, then their pair is either 1→3 or 3→5.
In particular a target-first selector at larger heights cannot
simply delete one vertex and reduce the count by exactly two.
Bigger or variable drops, and changed old arcs, are not excluded.

## Proof

By the inherited odd-cycle expansion (THM-026),

    H(T)-H(T-v) = 2 sum_{odd directed cycles C containing v} H(T-V(C)).

All complements have positive H, taking H(empty)=1. Hence
there is exactly one odd directed cycle containing v, and its
complement is transitive. If any cycle contains v, there is a
triangle containing v: follow the cycle from the first out-neighbor
of v to the last in-neighbor of v, and take the first transition
from an out-neighbor to an in-neighbor. So the unique odd cycle
is a triangle v→a→b→v.

Partition the other vertices into O beaten by v and I beating v.
Then a∈O and b∈I. The only arc from O to I is a→b, by triangle
uniqueness. There cannot be a two-edge path inside O ending at a,
as it gives a directed 5-cycle through v and b. Hence a has at
most one predecessor t in O: two predecessors would form such a
path in one direction. If t exists, t beats every other vertex
of O (other than a, which it already beats), for the same reason.
Dually b has at most one successor s in I, and if s exists every
other vertex of I beats it. Both t and s cannot exist, since
v→t→a→b→s→v would be a directed 5-cycle.

Therefore the active core is {v,a,b}, or that set plus exactly
one of t,s. All other I vertices dominate the entire core and
all other O vertices lose to the entire core. The remaining I
also dominates the remaining O. Both external blocks are
transitive as subsets of the triangle complement. So exterior
padding changes neither H(T) nor H(T−v). If neither t nor s
exists, deleting v from the core leaves a two-vertex transitive
tournament. If exactly one exists, deleting v leaves a directed
triangle (t→a→b→t or a→b→s→a). The assumed +2 then gives
the pairs 1→3 and 3→5. Both are actual, on three and four vertices.

This is the marked renewal obstruction suggested by the random
Lonely Runner owner-permutation and mod-6 acceleration sessions
recorded in UNIT_CLONE_CLASSIFICATION.md. It turns the inherited
unit-deletion certificate into a small-core classification.
