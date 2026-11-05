# A returning three-port ternary clock

## 1. Exact scope and result

A tournament is encoded by rows: bit `j` of row `i` is one exactly when `i -> j`. Let

`T = [2,12,145,149,99,143,175,19]`, with distinguished vertex `r=7`.

Replace `r` by a transitive chain of `m >= 1` vertices, each having precisely the same incidences with vertices `0,...,6` as `r`. Call the resulting actual tournament `T_m`. Then

**H(T_m) = 60 * 3^m - 1.**

Consequently `H(T_(m+1)) = 3 H(T_m) + 2` for every `m >= 1`. The first six values are `179,539,1619,4859,14579,43739`. All these tournaments are strong.

This is an actual returning physical interface: expose the same three old vertices `0,1,4` throughout. Add one vertex which beats exactly those three vertices and loses to every other old vertex, and retain all three old ports. This sends the actual graph `T_m` to the actual graph `T_(m+1)`, with exactly the same marked interface.

This is a THREE-port construction. It is not a claim about the earlier two-port language, arbitrary input tournaments, both ternary digits from one state, or the full Hamiltonian-path spectrum. Section 5 proves that a pure transitive-clone `3H+4` companion is impossible, including with more ports. Other physical operations are not ruled out.

## 2. Literal graph and physical operation

Label the clone chain `7,...,m+6` in forward order and put `B = sum(2^j : 7 <= j <= m+6)`. The first seven rows are

`[2,12,17|B,21|B,99,15|B,47|B]`.

Clone `j` has row `19 | sum(2^k : j < k <= m+6)`. These are disjoint bit sums, so bitwise OR and addition coincide here. This formula specifies every edge, without identifying an equal-response vertex with a possibly different physical vertex.

For ordered ports `[0,1,4]`, write `R(mask,drop)` for a rightward birth: the new vertex beats the selected old ports and loses to every other old vertex, then one entry is removed from the ordered list of the three old ports followed by the newborn. Removing a port only forgets its exposure; it does not delete a graph vertex. The returning word is the single operation

`R(7,3)`.

The newborn beats precisely `0,1,4`, all old clones point to it, and the retained ordered ports are still `[0,1,4]`. This proves graph-level iteration for arbitrary length.

The seed is strong, as direct reachability verifies. Alternatively, all vertices reach the directed cycle `7 -> 1 -> 2 -> 7`; from it one reaches `3` through `1`, `4` through `2`, `5,6` through `4`, and `0` through `2`. Every nonempty blowup of a vertex in a strong tournament is strong: directed paths in the quotient lift, and paths between vertices of the same block can go around a directed quotient cycle. Thus every `T_m` is strong.

## 3. Run decomposition and exact finite coefficient certificate

An `r`-skeleton is a directed word containing every vertex other than `r` exactly once and `r` exactly `k >= 1` times, with no consecutive occurrences of `r`. Let `a_k` count these words. Contracting each maximal clone-chain run in a Hamiltonian path of `T_m` gives exactly one such skeleton. Conversely distribute the `m` distinct clone vertices into its `k` ordered nonempty runs, placing each run in increasing chain order. There are `k! S(m,k)` distributions, where `S(m,k)` is a Stirling number of the second kind. Therefore

`H(T_m) = sum(k >= 1) a_k k! S(m,k)`.

Here the out-neighbors of `r` are `O={0,1,4}` and the in-neighbors are `I={2,3,5,6}`. The only edges from `O` to `I` are

`1 -> 2, 1 -> 3, 4 -> 5, 4 -> 6`.

Between successive copies of `r`, each nonempty internal fragment starts in `O` and ends in `I`. Along that fragment at least one edge crosses from `O` to `I`. Different fragments have disjoint vertices, hence selecting one such crossing from each gives a matching. This crossing graph has matching number exactly two: its edges use only left endpoints `1,4`, and `1 -> 2, 4 -> 5` are disjoint. Hence `a_k=0` for `k >= 4`.

The finite coefficient certificate is

`a_1=179, a_2=180, a_3=60`.

The supplied verifier obtains this twice, independently:

1. Literal enumeration: permute the seven ordinary vertices and choose a nonempty subset of their eight gaps in which to insert `r`. Test every adjacent directed edge. This enumerates every skeleton exactly once and also tests every possible `k=4,...,8`.
2. Subset dynamic programming on visited ordinary vertices, last vertex, and number of occurrences of `r`, allowing a transition to `r` only from an in-neighbor and a transition from `r` only to an out-neighbor.

The complete DP vector, including the zero-occurrence coordinate, is `[59,179,180,60,0,0,0,0,0]`. The zero-occurrence entry is simply `H(T-r)` and is not part of the cloning formula.

Since `1!S(m,1)=1`, `2!S(m,2)=2^m-2`, and `3!S(m,3)=3^m-3*2^m+3`, substitution gives

`179 + 180(2^m-2) + 60(3^m-3*2^m+3) = 60*3^m-1`.

The theorem for all `m` follows from this exact decomposition, not extrapolation from finite numerical values. Independent ordinary Hamiltonian-path subset counts for `m=1,...,6` are additional controls.

## 4. Why this literal append needs three ports

The returning newborn has exactly three old out-neighbors and `m+4 >= 5` old in-neighbors. A rightward birth with at most two exposed ports has at most two old out-neighbors. A leftward birth with at most two exposed ports has at most two old in-neighbors. Therefore this particular one-vertex extension cannot be either kind of two-port birth, even after choosing any different pair of old vertices as ports or relabeling the graph.

This is a distinction between exact physical operations. It does not exclude a multi-vertex two-port construction with the same scalar count, a different graph representation, or a more general replacement rule.

## 5. All-order obstruction to a pure-clone +4 companion

For an arbitrary finite tournament and root, `a_k` is divisible by `(k-1)!`. In an `r`-skeleton there are `k-1` nonempty internal fragments between consecutive occurrences of `r`. Permute these fragments while leaving the initial and final fragments fixed. Every resulting word remains directed, because every internal fragment starts at an out-neighbor and ends at an in-neighbor. The action is free: the fragments use disjoint nonempty sets of ordinary vertex labels. In particular, `a_3` is even.

Suppose some fixed rooted tournament had transitive-clone counts satisfying `H_(m+1)=3H_m+4` for all sufficiently large positive `m`. Then `H_m=C*3^m-2` on that tail. The run formula is a finite linear combination of `j^m`. If `K` is the largest index with `a_K != 0`, its coefficient at `K^m` is exactly the positive integer `a_K`. Linear independence of distinct positive exponential sequences forces `K=3` and `C=a_3`. (If `K<3`, the same comparison is impossible.) But `a_3` is even, making every `a_3*3^m-2` even, contrary to the classical oddness of a tournament's Hamiltonian-path count.

Thus no pure transitive-root-cloning clock, at any order or number of ports, can return `3H+4`, even eventually. In the specifically proposed coefficient form, `a_2=3a_3` and `a_1=3a_3-2` fail for the same parity reason. This obstruction does not concern non-clone gates, state-dependent words, or count-preserving rebasing.

## 6. Reproduction and limits

Run `python3 verify.py` and `python3 -O verify.py` from this directory. Both enforce all checks with explicit exceptions, not Python assertions. Their JSON receipts must agree. `core.py` and `verify.py` are self-contained and use only the Python standard library.

The proof does not need the discovery search or its bank of unlabeled tournaments. The explicit eight-vertex graph and its exhaustively checkable skeleton counts are the entire finite witness. No minimality claim about order, number of ports among all alternative constructions, or prime distribution of this sequence is made.
