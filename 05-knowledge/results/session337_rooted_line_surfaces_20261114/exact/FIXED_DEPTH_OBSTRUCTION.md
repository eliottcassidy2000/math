# No fixed-depth transitive cone covers the two-hub spectrum

Let $X>1$, put $L=\lfloor\log_2 X\rfloor+4$, and let $D_r(X)$ be the set of Hamiltonian counts at most $X$ obtained by a **two-hub endpoint** within $r$ tiling-line moves from its transitive class. The order can vary arbitrarily; so can the intermediate tournaments, colors of the lines, and endpoints' choice of hubs. Then
$$
 |D_r(X)|\le (L+1)^3\sum_{j=0}^{\min(L,12r)}\binom Lj3^j.
$$
In particular every fixed-depth cone gives only polylogarithmically many counts. More strongly, **for all sufficiently large $X$, the cone of depth $\lfloor L/192\rfloor$ still misses $\Omega(\sqrt X)$ counts at most $X$ already realized by strong one-hub tournaments**. No coverage conjecture is needed.

## Disagreement degree is an obstruction after projection

Lift the walk to fixed vertex names, transporting each path under an isomorphism, and reversing it if needed for the converse merge. The $i$-th move reverses every pair except the undirected path edge set $E_i$. Thus the endpoint agrees with its reference transitive order (the initial one or its reverse according to the walk's parity) except on $\triangle_i E_i$. As each $E_i$ has maximum degree two, at most $2r$ disagreement edges touch any vertex.

Remove the two hubs. The endpoint spine is transitive, as is the reference order restricted to its remaining vertices. The disagreement edges within the spine are exactly the inversions of their two total orders. Each vertex's displacement from its reference position is the difference of its inversion counts in the two directions, in magnitude at most $2r$.

For a hub, its reference out-neighbors in the reference spine are exactly a suffix after some cutoff $t$. Transfer that cutoff index to the endpoint spine. Only vertices within $2r$ reference positions on either side of the cutoff can be on a different side, so at most $4r$ membership bits change when replacing the reordered reference suffix by a true suffix. At most $2r$ additional arcs incident to the hub disagree with the reference. Thus its endpoint bits differ from a suffix in at most $6r$ positions. The two-hub spine word differs from the combination of **two simple cutoff columns** at at most $12r$ positions.

This count survives removing neutral padding: strip all singleton strong components except the two named hubs, as in the [S336 target inverse](https://github.com/eliottcassidy2000/math/blob/8d2d70178f9d17b191b93dcaffa6c1a8e25dbcb6/05-knowledge/results/session336_overlap_inverse_and_intercyclic_20261114/invert/TARGET_OVERLAP_INVERSE.md). The count is unchanged, the reference remains transitive, disagreement degree cannot increase, and the stripped spine length is at most $L$. This does **not** assert that the clipped walk is a valid metagraph walk: only its endpoint comparison is used.

At length $m\le L$, the two cutoffs have $(m+1)^2$ choices. A word with $j$ mismatches has at most $\binom{m}{j}3^j$ choices after fixing the cutoffs. This bounds the number of possible counts even if many words coincide, yielding the displayed inequality upon summing over lengths.

## Quantitative escape from the rooted cone

When $r\le L/192$, every summand obeys $j\le L/16$. Multiplying by $(1/16)^{j-L/16}\ge1$ before summing over $j$ gives
$$
 \sum_{j\le L/16}\binom Lj3^j\le (19/16)^L\,2^{L/4}.
$$
As $19^4<2\cdot16^4$, this is exponentially smaller than $2^{L/2}$, even with the polynomial factor. So the number of counts in this cone is $o(\sqrt X)$. The family of one-hub words $u0^{k+1}$, with $u$ ranging over the $k$-bit prefixes starting in 1, has $2^{k-1}$ distinct whole counts below $2^{2k+1}$, by the binary padding construction proved in [the second-branch note](ROOT_SECOND_BRANCH.md). This supplies the asserted $\Omega(\sqrt X)$ known attained counts outside the cone.

The result concerns the minimum depth of a **two-hub witness**. It does not rule out shallower witnesses needing more hubs, or a genuinely global count selector that lets depth grow. Its practical message is that a bounded neighborhood of the unique first blue line is too thin even when the terminal overlap may be chosen freely.
