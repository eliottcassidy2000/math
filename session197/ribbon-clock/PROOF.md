# Exact backward span of the binary cloning clock

Session 197. The two-port construction worker supplied the clock definition below; the path count and reordering bound are proved independently here. This is a graph-family obstruction, not an obstruction to realizing its counts by other tournaments.

## Definition

Let `C_n`, `n>=3`, consist of two vertices `r,a` and a transitive set `B` of size `n-2`, with

`r -> a`, `a -> B`, and `B -> r`.

Adding a new last vertex of the transitive set `B` doubles the old Hamiltonian-path count minus one. It is the fixed binary-affine clock obtained by retaining the same root `r` indefinitely.

For a supplied linear ordering of a tournament, its backward span is the maximum distance between the endpoints of a backward edge (zero if no edge is backward). The minimum backward span minimizes this over **all** vertex orderings.

## Path count

**Proposition.** `H(C_n)=2^(n-2)+1`.

If `r` is not last in a Hamiltonian path, it must be immediately followed by `a`, since `a` is its only outneighbor. Split `B` arbitrarily into the vertices before `r` and those after `a`. Each part must appear in its unique transitive order, and every such split is valid. This gives `2^(n-2)` paths, including the empty-left and empty-right cases.

If `r` is last, `a` must be first because its only inneighbor is `r`. There is exactly one additional path, namely `a`, all of `B` in transitive order, then `r`.

Hence `H(C_(n+1))=2H(C_n)-1`.

## Reordering cannot keep the span bounded

**Theorem.** The minimum backward span of `C_n` is exactly

`ceil((n-1)/3)`.

Consider an arbitrary ordering and let the positions of `a,r` be `p,q`. Write its backward span as `w`. The graph contains a directed triangle, so `w>=1`.

If `p<q`, three forced backward gaps give

- `w>=p-1`: when `p>1`, the first vertex lies in `B` and is beaten by `a`;
- `w>=q-p`: the edge `r -> a`;
- `w>=n-q`: when `q<n`, the last vertex lies in `B` and beats `r`.

The zero-length boundary cases are harmless. Summing gives `3w>=n-1`.

If `q<p`, let `e,l` be the first and last positions occupied by `B`. Because `B` has `n-2` vertices in distinct positions, `l-e>=n-3`. The edges `a -> B` and `B -> r` give

`w>=p-e`, and `w>=l-q`,

where a negative right-hand side is simply vacuous. Therefore

`2w >= (p-q)+(l-e) >= 1+(n-3)=n-2`.

Since `w>=1`, it follows again that `n-1<=2w+1<=3w`. This proves the lower bound without assuming that the vertices of `B` keep their transitive order.

For the upper bound, set `w=ceil((n-1)/3)`. Choose nonnegative gaps `g_1,g_3` and a positive gap `g_2`, all at most `w`, with

`g_1+g_2+g_3=n-1`.

Such a choice exists because `n>=3` and `3w>=n-1`. Put `a` at position `g_1+1` and `r` at position `g_1+g_2+1`; fill all other positions with `B` in its transitive order. The only backward edges are from `a` to earlier members of `B`, the edge `r -> a`, and edges from later members of `B` to `r`. Their lengths are at most `g_1,g_2,g_3`, respectively. This attains the lower bound.

## Exact scope of the obstruction

At a fixed positive integer width `w>=1`, this exact graph family can continue only through order `3w+1`, and its largest count is `2^(3w-1)+1`. At width three this means order at most ten and count at most 257.

The obstruction survives arbitrary reordering, but it does **not** survive replacing the graph by a different graph with the same count. For example:

- `C_11` has count 513 and minimum backward span four.
- The distinct span-2 word `022333332` has count 513 on nine vertices.

The latter word comes from the independently verified bounded-output closure in `../ribbon/closure-2-1001.json`; this packet also checks it directly by Held–Karp dynamic programming. Thus a count-preserving compiler can sometimes succeed even when a reordering-only compiler is impossible. No uniform such compiler is established here.

## Checks

Run `python session197/ribbon-clock/verify.py` from the repository root. Python 3 only; assertions must remain enabled.

The verifier enumerates every ordering of `C_n` through `n=8`, checks the exact minimum span, independently counts Hamiltonian paths through `n=12`, checks the balanced construction through order 200, and directly verifies the contrasting count-513 witness.
