# Prime induced jumps of sixteen and the 17 gluing row

**Result.** For prime $q\ge17$, there is an induced increase from $q$ to $q+16$ exactly for $q\in\{17,19,23,29\}$. In every such increase, after deleting the host's singleton strong components, the host is the old strong core with one mixed-contact vertex added. All four possibilities admit an induced 17-path witness distinct from the induced $q$-path witness. The [11-row](../eleven_row/README.md) and [13-row](../thirteen_row/README.md) results exclude the smaller $q=11,13$ from the induced 17 gluing row by symmetry. Thus the complete 17 row of induced $p+q-1$ has those four partners. Nothing here rules out gluing with rewiring.

## Exterior block normal form

Let the old tournament be a strong core $C$ with $q=H(C)$, and singleton components in its ordinal order. Insert $x$ with mixed contact to $C$. Write $r$ for the full half-increase and $d$ for the half-increase in $C+x$ alone. Assume $r<q$.

Use the [positive odd-cycle-packing formula](https://arxiv.org/abs/2307.05569), in the [rooted-insertion form](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md). A triangle through $x$ and two pads contributes $q$, and is impossible. Call a pad before $C$ active when $x$ beats it, and a pad after $C$ active when it beats $x$. There cannot be active pads on both sides. On the left an active pad cannot precede an inactive pad, and on the right an inactive pad cannot precede an active pad. Each forbidden arrangement would give a root triangle through two pads. After dropping outer inactive singleton components, the old exterior is therefore either empty or a contiguous block of $j\ge1$ active pads on one side.

For a prefix block, let $O,E$ be the sums over odd- and even-sized nonempty directed paths in $C$ ending at a vertex beating $x$. Each path on vertex set $S$ is weighted by $H(C-S)$. For a suffix block, use paths starting at a vertex beaten by $x$. Every exterior rooted odd cycle uses an increasing sequence of pads from that block, and one such core path, with the two vertex counts of the same parity. There are $2^{j-1}$ odd and $2^{j-1}-1$ nonempty even choices of pad subsets. Thus
$$
r-d=2^{j-1}(O+E)-E.
$$
In particular a positive residual $s$ is filled on that side if and only if $(s+E)/(O+E)$ equals $2^k$ for an integer $k\ge0$, provided the resulting full $r$ is below $q$. The reverse construction just adjoins this active block in ordinal order and orients all its edges to $x$ accordingly. This is an exact inverse criterion in its regime, not just a triangle lower bound. If the exterior is empty, the residual is zero.

## Reduction and replay

Before the first strict insertion the unique nontrivial strong component $C$ of the prime witness is unchanged. Uniform contact from $x$ would preserve $C$ as a module and lift every path from a quotient with a new directed triangle. This costs at least $2q$ and is impossible here. The triangle budget forces at most eight triangles in $C$, and the strong order bound gives at most ten vertices. Relabel every core along a Hamilton path. Requiring $i\to i+1$ and recursively pruning any prefix once it has nine triangles retains every eligible isomorphism class. There are 100,657 such fixed-path encodings at order ten, of which 10,300 are strong.

For each prime $q\ge17$ compute the local half-gap $d$ for each mixed contact, by counting even-sized paths from an outneighbor of $x$ to an inneighbor, with the same complement weights as above. These paths are exactly the possible internal rooted odd cycles with $x$ deleted. The complete list with $d\le8$ is:

| $q$ | $d=5$ | $d=6$ | $d=7$ | $d=8$ | Exact cases with a separate induced 17 |
|---|---:|---:|---:|---:|---:|
| 17 | 51 | 68 | 34 | 119 | 85 |
| 19 | 0 | 38 | 76 | 38 | 38 |
| 23 | 0 | 0 | 46 | 46 | 46 |
| 29 | 0 | 0 | 0 | 58 | 58 |

The entries count encodings, not isomorphism classes. The first strict step costs at least ten, leaving at most six. The earlier all-order short-gap exclusions rule out a second strict step since there is no positive induced increase of 2, 4 or 6 from count at least 11. Therefore the strict insertion must give exactly 16. For every encoded contact with $d<8$, test independently on each side whether $(8-d+E)/(O+E)$ is a nonnegative integer power of two. None passes. When $d=8$, the budget excludes an active exterior, and direct subset DP verifies the whole and induced counts. All later neutral additions preserve the strong core. This proves both classifications.

Run the dependency-free checker with **node check_seventeen.js**. Literal witnesses are given below; a row mask encodes the outgoing neighbors and vertex 0 is the low bit. The old $q$-witness omits the last vertex. The indicated mask selects a different induced 17-witness.

| $q$ | Host paths | Outgoing masks | Mask for 17 |
|---|---:|---|---:|
| 17 | 33 | 94,124,120,112,96,1,32 | 111 |
| 19 | 35 | 30,124,24,48,32,5,61 | 111 |
| 23 | 39 | 94,28,56,48,32,3,62 | 125 |
| 29 | 45 | 190,252,248,16,32,72,153,56 | 123 |

The structural takeaway for the bounded fresh deposit is that, below the core cost, extra exterior length gives at most two dyadic offset sequences per fixed mixed contact. It cannot behave as a freely programmable additive window. To use exterior variation for a prescribed target, either change the two endpoint profiles, or leave the one-core/single-insertion regime; the theorem says nothing against those routes.
