# Literal transitive cloning cannot stay at fixed backward span

Session 197. This note concerns replacing one fixed vertex by a growing transitive module. It does not exclude a different, count-preserving graph construction.

## Definitions and statement

For a nonempty tournament `T`, let `H(T)` count directed Hamiltonian paths. Let `bw(T)` be the minimum, over all linear orders of its vertices, of the largest distance between the endpoints of a backward edge; the maximum over no backward edges is zero.

Fix `v` in `T` and an integer `m>=1`. Write `T_m=T[v <- Tr_m]` for the tournament replacing `v` by transitive vertices `b_1 -> ... -> b_m`, all having exactly the former orientations of `v` toward every other old vertex.

**Theorem (complete dichotomy).**

1. The following are equivalent: `v` belongs to a directed triangle; `v` belongs to a strong component with at least two vertices. In a tournament a nontrivial strong component actually has at least three vertices.
2. If `v` belongs to no directed triangle, then `H(T_m)=H(T)` for every `m>=1`.
3. If `v` belongs to a directed triangle, then, for every `m>=1`,

   `bw(T_m) >= ceil((m+1)/3)` and `H(T_m) >= 2^m+1`.

In particular, a literal transitive-cloning family with nonconstant path count eventually leaves every bounded-span class, even if each graph may be arbitrarily reordered. At a fixed positive width `w`, its clone block has size at most `3w-1`. The bound is sharp as a bound depending only on `m`: take `T` to be a directed triangle. At width zero no triangular case is possible.

No assertion is made that an individual cloning step always increases `H`, nor that alternative tournaments with the same count have large backward span.

## A vertex is in a triangle exactly when its component is nontrivial

Put `I={x:x->v}` and `O={x:v->x}`. If `v` is in no triangle, every cross edge is `I->O`: an edge `o->i` would give `v->o->i->v`. Thus the whole tournament has the ordered decomposition

`T[I] -> {v} -> T[O]`.

No vertex other than `v` can both reach and be reached from `v`, so its strong component is a singleton. Conversely, the three vertices of any directed triangle are mutually reachable, and therefore lie in a common nontrivial strong component. These two arguments prove both directions for arbitrary order and without assumptions on other components.

Every Hamiltonian path in an ordered join `A->B` consists of a Hamiltonian path of `A` followed by one of `B`. Therefore replacing the middle singleton in the displayed decomposition by a transitive block changes no path count. Here empty sides contribute the empty-product factor one; this is only a product convention, not a change to any implementation's empty-input path-count convention. This covers a singleton `T`, sources, sinks, and singleton components between larger components.

## The forced induced clock and its span

Choose a triangle `v->x->y->v`. In `T_m`, the set `{x,y,b_1,...,b_m}` induces exactly the clock with

`r=x`, `a=y`, `B={b_1,...,b_m}`,

and edges `r->a->B->r`, with `B` transitive. Call this graph `C_(m+2)`.

Minimum backward span is monotone under taking induced subtournaments: restrict any ambient order to the selected vertices. Its backward edge orientations are unchanged, and deleting other vertices only decreases the distances. Minimize over the ambient order to obtain `bw(T_m)>=bw(C_(m+2))`.

For completeness, `bw(C_n)=ceil((n-1)/3)` for every `n>=3`. Consider any order, with `a,r` at positions `p,q`, and span `w>=1`.

- If `p<q`, the forced backward gaps are at least `p-1`, `q-p`, and `n-q`: respectively `a` to an earlier member of `B`, `r->a`, and a later member of `B` to `r`. Boundary gaps zero need no edge. Their sum gives `3w>=n-1`.
- If `q<p`, let `e,l` be the earliest and latest positions of `B`. There are `n-2` such positions, so `l-e>=n-3`. The orientations `a->B->r` give `w>=p-e` and `w>=l-q` (a negative right-hand side is vacuous). Consequently `2w>=(p-q)+(l-e)>=n-2`. Together with `w>=1`, this implies `3w>=n-1`.

For the matching upper bound, choose three gaps `g_1,g_2,g_3<=ceil((n-1)/3)`, with `g_1,g_3>=0`, `g_2>=1`, and sum `n-1`. Put `a,r` at positions `g_1+1,g_1+g_2+1` and fill the other positions with `B` in transitive order. All backward edges have length at most the relevant gap. Apply the identity at `n=m+2`.

## The path-count lower bound

The clock has exactly `H(C_(m+2))=2^m+1` paths. When `r` is not last, it must be followed by its sole outneighbor `a`; choose any subset of `B` to precede `r`, with that subset and its complement each in transitive order. These give `2^m` paths. When `r` is last, the unique additional path is `a,B,r`.

Path count is monotone under induced extension. To see this directly, insert a new vertex `z` into an existing path `p_1,...,p_s`. If `z->p_1`, prepend it. Otherwise insert it before the first `p_i` with `z->p_i`, or append it if none exists. The predecessor of the insertion points to `z`. This deterministic extension is injective because deleting `z` recovers the original path. Repeating the operation proves the claimed induced monotonicity and hence `H(T_m)>=2^m+1`.

## Same-count replacement remains possible

The bound is about the cloned graph, not its scalar count. Two directly checkable strong witnesses illustrate the distinction. The word convention is oldest-first: in a span-2 word, bit zero reverses the edge to the vertex two positions earlier, and bit one reverses the edge to the immediately preceding vertex; nonexistent initial positions are ignored.

- Word `022333332` has order 9, backward span at most 2, and path count 513. The clock with that count is `C_11`, whose minimum span is 4.
- Word `003231331133133` has order 15, backward span at most 2, and path count 65,537. The clock with that count is `C_18`, whose minimum span is 6.

Both are checked by graph-level subset dynamic programming, separately from the frontier search that found the second word. No all-target construction, minimal-order claim, or minimum-span claim for either alternative witness is needed here.

## Verification and provenance

Run `python session197/ribbon-cloning/verify.py` from the scratch root. The verifier is standalone Python 3 and uses explicit checks rather than removable assertions. It exhausts small tournaments for the triangle/component equivalence and constant-count case, checks literal induced clock realizations, exhausts all orders of small cloned graphs for the lower bound, checks sharpness of the clocks, and recounts both alternative witnesses independently.

The clock identity agrees with the previously independently reviewed `../ribbon-clock/PROOF.md` (SHA-256 `8e7a7fca634179d970795b08fc6a79a21390a2c037d50fa05dbb1c82fd7c8a01`). Its proof is reproduced above, so the present mathematical argument is self-contained. The broader cloning dichotomy is a consequence, not a claim of priority or novelty.
