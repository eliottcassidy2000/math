# Induced gluing with 19 and a small-jump strong-core bound

**19 row.** For prime $q\ge11$, a tournament with $19+q-1$ Hamiltonian paths can contain induced witnesses of counts 19 and $q$ exactly for $q\in\{17,19\}$. More generally, an induced +18 jump from a prime $q\ge19$ exists exactly for $q=19,23$; the $q=23$ hosts cannot retain an induced 19. These are statements about induced retention, not about what counts can be realized after rewiring.

**General bound.** Suppose the source has exactly one nontrivial strong component and path count $q$. If a host retains it induced and has $q+2r$ paths with $0<r<q$, then the host has exactly one nontrivial strong component, on at most $r+3$ vertices. The same conclusion applies to a prime-count source automatically. All induced counts in the host survive passage to that component, in the sense that for every induced vertex set its intersection with the component has the same count, with the empty count defined as 1. In particular, for induced gluing with odd $3\le p\le q$ and prime $q$, a host of count $p+q-1$ exists with both witnesses if and only if one exists on at most $(p+5)/2$ vertices.

## Proof of the general bound

Insert the other vertices in any order. The [positive odd-cycle-packing formula](https://arxiv.org/abs/2307.05569) implies induced monotonicity, and a neutral insertion creates no directed triangle. So a neutral insertion adds just another ordinal singleton without enlarging the core. At a strict insertion let $C$ be the current core and $2e$ the increase. Uniform contact to $C$ keeps it a module and, if strict, creates a directed triangle in the quotient. All quotient Hamilton paths then lift, giving at least $3H(C)$ paths and spending at least $H(C)$ in half-gap. This is impossible because $e\le r<q\le H(C)$.

Mixed contact instead gives the [one-sided exterior block normal form](../seventeen_row/README.md): since $e<H(C)$, after dropping neutral pads the new core consists precisely of $C$, the new vertex, and a single block of $j\ge0$ active pads. Each active pad lies on a rooted triangle meeting $C$; mixed contact makes $C$ plus the new vertex strong. Let $d$ be the local half-increase with no pads. The [triangle budget](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md) gives $d\ge\tau(C)\ge |C|-2$. The $j$ different rooted triangles show $e-d\ge j$. Consequently the new core has
$$
|C|+1+j\le e+3\le r+3.
$$
Apply at the last strict step. Ordinal singleton components stay ordinal and singleton on restriction, so deleting them changes no induced count. This completes the proof, including for an arbitrary number of strict steps. The additive constant three cannot be lowered uniformly: split one vertex of a directed triangle into an ordered pair. The resulting strong tournament has two directed triangles and no longer odd cycles, so the count increases from 3 to 5 and the core order is $1+3=4$.

Two useful refinements are immediate. Use the actual last half-increase $e$, not $r$, when earlier steps have spent some budget. And for $j\ge1$ the exact exterior formula gives $e-d=2^{j-1}O+(2^{j-1}-1)E\ge 2^{j-1}$, since the mixed contact has $O\ge1$ and $E\ge0$. Thus the new order is at most $e+3-(2^{j-1}-j)$. The formula also gives an exact finite successor rule: at each current core, choose a mixed cut, then no pad block or one active block with a feasible dyadic offset. Neutral padding need never be enumerated. In particular, a small jump cannot be obtained by delaying many individually cheap strict insertions. If the initial core order is $n_0$ and there are $k$ strict insertions, the $i$-th core has at least $n_0+i$ vertices. Summing the per-step costs gives $r\ge k(n_0-3)+k(k+1)/2$. The last step alone also spends at least the final core order minus three.

## Exact replay for +18

For prime $q\ge19$, at the first strict insertion the old core has at most nine triangles and order at most eleven. This uses the same triangle budget and strong order bound. Enumerate Hamilton-path orderings for the cores, fixing the consecutive forward arcs and rejecting prefixes with more than nine directed triangles. This loses no eligible isomorphism type. It leaves 504,088 encodings at order eleven, 36,684 strong. The checker first uses the sum of rooted triangle weights as a lower bound to prune the mixed cuts; every surviving local half-gap is computed using exact weighted even paths. The complete eligible table is:

| $q$ | $d=6$ | $d=7$ | $d=8$ | $d=9$ | Exact cases with a separate induced 19 |
|---|---:|---:|---:|---:|---:|
| 19 | 38 | 76 | 38 | 76 | 76 |
| 23 | 0 | 46 | 46 | 184 | 0 |
| 29 | 0 | 0 | 58 | 0 | 0 |

The counts are fixed-Hamilton-path encodings, not isomorphism classes. The first strict step costs at least 12, leaving at most 6, so there cannot be another strict step by the [all-order exclusions of +2,+4,+6 from counts above 9](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md). Hence the first step must give all 18. For each $d<9$, the two exact exterior ratios $(9-d+E)/(O+E)$ fail the nonnegative power-of-two test, so padding never closes a partial hit. At $d=9$, direct subset DP checks the exact augmented count and every induced subset. This proves the +18 statement and the $q\ge19$ part of the row. The [11](../eleven_row/README.md), [13](../thirteen_row/README.md) and [17](../seventeen_row/README.md) rows settle smaller primes by symmetry.

The seven-vertex self-pair has outgoing bit masks $(62,124,56,48,96,64,13)$, vertex 0 in the low bit. Its full count is 37; the two input masks 119 and 123 each have count 19, and their intersection 115 has count 9. Run **node check_nineteen.js** to reproduce the table and this literal witness.

The strongest consequence for the larger project is the scope of this obstruction: the finite core bound holds at all orders of the host, not just for previously generated tournaments. It rules out certain ways of retaining the old graph. It does not rule out preserving the old continuation data while replacing or rewiring the old graph, which is what the bounded fresh deposit requires.

## The shield threshold and ternary digits

There is a simple structural change at the edge of this regime. If a host retains the old strong core $C$ and has any second nontrivial strong component, ordinal multiplication forces its count to be at least $3q$. At equality, the component containing $C$ has count $q$ and the other has count 3. A strong extension of $C$ with the same count is impossible: adding the extra vertices with no increase creates no directed triangle, so they are singleton components. The count-3 strong component is a directed triangle (a strong tournament of order $s$ has at least $s-2$ triangles, each contributing two paths on top of the empty packing). Therefore the only host with *separate* nontrivial components at the threshold, after deleting singleton pads, is an ordinal sum of the original core and a triangle, in either order.

For the two nearby targets $3q+2$ and $3q+4$ with $q\ge3$, there can be no second nontrivial component at all. The product of outside component counts would be an odd integer at least three, and less than five since $3q+4<5q$. Hence the count would be divisible by three, a contradiction. In the induced setting a literal triangle shield can implement the zero ternary digit, but the nonzero digits cannot come from a disconnected shield. They require interaction with the strong core (or an operation that no longer preserves the old core as induced). This does not by itself forbid a repeatable fresh deposit.
