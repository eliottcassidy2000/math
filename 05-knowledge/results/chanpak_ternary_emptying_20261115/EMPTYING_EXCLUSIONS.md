# Empty the single outlet: the exclusions 7 and 21 and the 63 witness

This is an alternative to the [Busch seam certificate](../chanpak_window_descent_20261115/SMALL_EXCLUSIONS.md), not a novelty claim for the exclusions. It replaces 144 seam completions with a rooted deletion recursion, starting from order four; the companion cutoff removes the need for a stronger extremal bound. The same recursion explains a known actual 63 witness. The general two-wing formula below makes the relation between its empty core and the Moon–Busch 5-blocks explicit. It does not establish spectral completeness.

## The one-outlet emptying rule

Let H(W) count directed Hamilton paths, with H(empty)=1. For a vertex u of W define

~~~text
R_u(W;z) = sum over paths L starting at u of H(W − V(L)) z^|V(L)|.
~~~

Write E_u and O_u for its even and odd coefficients summed, respectively; S_u=E_u+O_u. These are weighted sums, not numbers of full paths from u. Every directed odd cycle through a newly added vertex v that beats only u opens on removal of v into an even-order path out of u. Conversely every such even-order path closes at v because its final vertex is not u and therefore beats v. By the [odd-cycle packing expansion](https://arxiv.org/abs/2307.05569), comparing packings with and without a cycle through v gives the exact rule

~~~text
H(W+v)=H(W)+2E_u(W).
~~~

Moreover W+v is strong exactly when u lies in the first SCC of W. Indeed all other vertices beat v, and u reaches the entire W exactly when it is in the first SCC. This explains how a rooted restriction helps the exclusion rather than merely giving an exact formula.

The two parities can be computed without constructing new tournament classes. Split a rooted path by its next vertex:

~~~text
R_u(W;z) = z * [H(W−u) + sum over u→w of R_w(W−u;z)].
~~~

Here a root singleton has constant weight H(W-u) before multiplication by z. Thus E_u(W)=sum O_w(W-u), and O_u(W)=H(W-u)+sum E_w(W-u), with the sums over successors of u. H itself is computed by the usual partition by first vertex; the replay also directly computes the weighted rooted polynomials by grouping by the path vertex set.

For W=A→B with B nonempty and root u in A, a path from A has an optional B-part, so

~~~text
E_u(W) = S_u(A) Q_2(B)/2,
Q_2(B) = sum over subsets X of B of H(B[X]) H(B−X).
~~~

To justify the factor 1/2, the signed sum of the optional B-part vanishes. Expand both Hamilton counts in the signed two-slot convolution sum_X (−1)^|X| H(B[X]) H(B−X) by odd-cycle packings. Every selected odd cycle and every leftover vertex is an odd-size block that can lie in either slot, giving signed factor zero. In a nonempty tournament at least one such block exists. This is an integer cancellation; it is not a congruence argument.

## Small strong cores, via deletion rather than cross-seam cases

The condensation of any tournament is transitive, so H multiplies over its strong components. It is sufficient to exclude strong counts 7 and 21: the only possible nontrivial factorization of 21 needs a factor 7. For the large-order cutoff one may use the companion [elementary order-seven bound](SEVEN_CUTOFF.md) in place of the sharper Moon–Busch bound. The strong tournaments of order three and four are C3 and the unique near-transitive four, with H=3 and 5 respectively.

At both order five and six, either there is a vertex of outdegree one (after taking the converse if necessary) or every in- and outdegree is at least two. In the latter case order five is the unique regular C5 tournament: starting with a vertex r, its predecessors and successors each form an ordered pair, and degree two forces the four cross arcs. Concretely, write successors a→b, predecessors c→d. Vertex b must beat c,d; c must lose to a, forcing a→c; then d→a. This yields the circulant with steps 1 and 2. It has five triangles by the degree formula and two directed 5-cycles (the steps round a Hamilton cycle must be all 1 or all 2), giving H=15. At order six the degrees are three 2s and three 3s, giving eight triangles. Each single deletion has min in- and outdegree at least one, so its first and last SCCs cannot be singletons; if it were not strong it would need at least six vertices, not five. Hence every single deletion has at least one directed 5-cycle by the order-five routes below, giving H≥1+2(8+6)=29.

In the low-degree case remove the one-outlet vertex, then use the first-SCC criterion. For order five the root has exactly seven tests before the converse: attach to a source of TT4, to the source of source→C3, to one cyclic root of C3→sink, or to one of the four roots of the unique strong four. Add the unique regular case. These tests follow directly from the complete order-four classification (acyclic, one cyclic SCC and a singleton in either order, or strong). In the table the entries are outgoing bit masks, indexed from zero. The seven low-degree choices come first.

| Constructed strong five | H | E palette over its roots and those of its converse |
|:--|--:|:--|
| (14,28,24,16,1) | 9 | 3,4,7,12 |
| (14,20,24,18,1) | 15 | 6,8,9 |
| (10,28,25,16,1) | 9 | 3,4,5,10,11 |
| (6,28,24,17,1) | 11 | 4,6,7,9,13 |
| (22,12,24,17,2) | 13 | 5,6,8,12 |
| (22,28,8,17,4) | 9 | 3,4,5,10,11 |
| (22,28,24,1,8) | 9 | 3,4,7,12 |
| (6,12,24,17,3) | 15 | 9 |

Duplicates are allowed; the rooted recurrence, not a general orientation enumeration, gives the entries. For the 5-cycle check just close a path starting at the selected root of the 4-host through v: in the first three hosts use 0123, in the fourth host use successive rotations of 0123. The regular C5 has the step-1 cycle. All converses have the reversed cycles. Thus order five never gives 7 or 21. For a one-outlet extension of a **strong** five, H+2E=21 would require E=6,5,4,3 when H=9,11,13,15 respectively. Each is missing from its row union.

Only decomposable five-vertex hosts remain. If their first SCC is a singleton, then its order-four tail contributes the pairs (H,Q_2)=(1,16),(3,24),(5,32), and the new counts are 17,27,37. (On order four write Q_2=16+8c_3 in the two-slot expansion.) If the first SCC is C3 and the tail is TT2, use S_u(C3)=3 and Q_2(TT2)=4 for new H=3+12=15. If the first SCC is the unique strong four (outgoing masks (6,12,8,1)) and the tail a singleton, the S palette is {5,6,7,9}, so the possible new counts are 15,17,19,23. Thus no decomposable host reaches 21 either.

The finite table is local; the unbounded exclusion comes from the product reduction and the strong-order cutoff, not from an assumed growth pattern.

## What opens at 63: two wings with a common junction

For a,b≥2 take transitive rows A' and B' of respective lengths a,b, with A' dominating B'. Add a junction r that beats only the first vertex in A' and loses only to the last vertex in B'. Each wing together with r is nearly transitive. All cycles go through r, and every complement is transitive. Within a wing, an odd cycle has one forced endpoint and an odd subset of the remaining vertices; a crossing cycle has two forced endpoints and an even subset of the other vertices. Their counts are (2^(a−2)), (2^(b−2)), (2^(a+b−3)). Thus with X=2^(a−1), Y=2^(b−1),

~~~text
H(W)=(X+1)(Y+1).
~~~

The weighted path polynomial rooted at r, marking only vertices *after* r, is

~~~text
R_r(W;z)/z=z(1+z)^{a+b-1}+(1+z)^b-z.
~~~

If the path enters A' it must start at its first vertex, with unrestricted subsequent subset of A' and B'; otherwise any B' subset except the singleton last vertex is allowed. Consequently E_r is the odd part in the displayed polynomial, (XY+Y-1). A new one-outlet vertex beating only r gives

~~~text
H(W+v)=(X+1)(3Y+1)-2.
~~~

Taking a=b=3 gives 63 = 25+2·19, **not** 3² times a 7-path component. Deleting the outlet from the [known 63 witness](../chanpak_two_path_session/ONE_OUTLET.md) with rows (4,5,152,3,11,31,63,123) leaves exactly this two-wing core: reindexed rows (2,76,1,5,15,31,61), r=1, A'=(6,5,4), B'=(3,2,0).

Each individual near-transitive wing of lengths 2,3,4 gives a block with 3,5,9 paths respectively. The two-wing product consequently also displays the familiar 15=3·5 and 45=5·9 via the same common junction (and 5 is a single strong-four block). The exposed root stores information that the product alone misses. For (a,b)=(2,3),(3,2) the deleted cores both have H=15, yet adding the same type of one-outlet vertex gives 37 and 33. Reversing the order of the wings is an actual arithmetic control, not just a re-expression of 15.

The inverse for this two-wing family is literal but sparse: factor H+2 and ask whether its factors are 2^(a−1)+1 and 3·2^(b−1)+1, with a,b≥2. The rule wins at 63 and leaves most counts unresolved. A general construction must find physically available weighted root responses, not assume orbit-wise coverage from 7 or 21 being excluded.

Reproduce the tables and the two-wing numerics with `python replay_empty_spectrum.py` (same receipt under `-O`).
