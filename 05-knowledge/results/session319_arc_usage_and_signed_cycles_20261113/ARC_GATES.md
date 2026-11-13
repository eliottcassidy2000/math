# Dark arcs, cycle cuts, and a reusable two-vertex gate

A Hamilton arc weight W_T(u,v) counts spanning directed paths using u→v.
A zero is an arc of the tournament, not an absent adjacency. We use
*dark* and *nontraceable* synonymously. The classical exact structure
theorem is [Busch–Jacobson–Reid](https://doi.org/10.1002/jgt.20173);
see also its [primary restatement](https://combinatorialpress.com/article/jcmcc/Volume%20063/vol-063-paper%201.pdf).
The repo's [guard reconstruction](../session70_attachment_guard_alignment_20261026/session70-guards/REPORT.md)
already proves that the **directed** positive-support graph determines
the full tournament. We do not claim to discover nontraceable arcs or
reprove that classification. The point here is a quantitative local
test and a repeatable, exact weight-preserving construction.

## A quantitative cut test on any Hamilton cycle

Suppose a spanning cycle of T has cyclic word u,X,v,Y, and consider
the chord u→v. The Y side is nonempty if this is not an arc of the
cycle, since otherwise v→u would lie on that cycle. The X side is
also nonempty for a chord. The chord closes the cycle C'=u,v,Y
on the complement of X. Write a(p),b(p) for the start and end of a
Hamilton path p of T[X]. Then
```
W_T(u,v) >= sum_p ( #{z in C', z != v : b(p) -> z}
                  +#{z in C', z != u : z -> a(p)} ).
```
In the first case append the cycle cut to start at z, in the second
prepend it cut to end at z. The excluded endpoints are precisely the
cuts that would remove u→v. The outputs are all distinct: X occupies
the prefix or suffix, its local path and the cut position are recoverable,
and a nonempty proper prefix and suffix of the same length cannot have
the same vertex set.

In particular if T[X] is strong (or a singleton), every x in X can
be the start and the end of local Hamilton paths. For each x in X
and z in Y, whichever orientation exists supplies at least one
summand. Thus W_T(u,v)>=|X|·|Y|. More finely the bound is at least
|Y| times the sum over x of min(local starts at x, local ends at x).
It follows that a dark chord must skip a decomposable segment on
**every** spanning cycle. A chord skipping just one vertex forward
cannot be dark. No dark arc exists in a strong tournament through
order four, since two intervening vertices and at least one return
vertex are required.

## The two-vertex gate

In a tournament T with marked arc u→v, adjoin x,y by making u the
sole in-neighbor of x, v the sole out-neighbor of y, and orient x→y.
These rules agree on the shared pair. Any spanning path containing
u→v in the enlarged graph must start at x and end at y: using the
unique predecessor of x or successor of y would conflict with u→v.
Deleting x,y therefore identifies its marked paths bijectively with
the spanning paths of T containing u→v strictly internally (u not
first and v not last). Conversely the two outer arcs exist precisely
for those internal paths.

If T is strong, the enlarged graph is strong by the incoming and
outgoing links to x,y. After one gate every marked occurrence is
internal, so further gates preserve the **same** arc weight exactly,
even though the total Hamiltonian count need not stay fixed. Thus
the gate is an idempotent projector on the four-cell marked response
(start/end flags); on the graph it can be repeated without bound.
Over a directed triangle it gives a dark arc at every odd order
>=5. Over a cyclic wrap w→a→TT_m→w, the marked arc w→a has
exactly 2^m−2 internal paths: pick any nonempty proper subset of
the chain before w, and place the remainder after a, with both
parts increasing. Further gates preserve this count as well.

For a dark chord skipping **exactly two** vertices x,y in a spanning
cycle, the gate is not merely sufficient; it is forced. The cut
inequality implies x beats every vertex of the complementary cycle
except u, and every vertex of that cycle except v beats y. Since
the original cycle contains u→x, x→y, y→v, those are precisely the
gate orientations. Thus the zero question there is reduced exactly
to absence of an internal occurrence of u→v in the smaller
complementary tournament. This is a local inverse theorem at the
shortest possible dark gap, not a claim to classify longer gaps.
