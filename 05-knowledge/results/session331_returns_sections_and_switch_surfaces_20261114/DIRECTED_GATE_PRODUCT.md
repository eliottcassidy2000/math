# Directed switches and an exact product of endpoint counts

For a tournament \(A\) and marked entry, exit \(s,t\), let \(p(A;s,t)\) be the number of spanning directed paths from \(s\) to \(t\) (one for a singleton). This is not its unmarked Hamiltonian count.

Arrange marked tournaments \(A_1,\ldots,A_k\) in a line. Orient all cross arcs from later blocks to earlier blocks except for the single bridge \(t_i\to s_{i+1}\) across each adjacent interface. The complete graph is a tournament. A path from \(s_1\) to \(t_k\) must cross each cut from left to right exactly one more time than in the other direction. The only forward arc on that cut is the designated bridge. Thus the path never goes backwards across a cut, and visits every block contiguously at its prescribed ports. Conversely the local paths concatenate across the bridges. The exact fixed-port count is therefore
\[
p(A_1;s_1,t_1)\cdots p(A_k;s_k,t_k).
\]
For nontrivial strong interior blocks and singleton buffers at both ends, the construction is strong: each block can use its next bridge, and later blocks have return arcs to earlier ones (in particular from the last singleton to the first). Cycles **containing the marked closing arc** of the last singleton to the first are in bijection with these fixed-port paths. Cycles not containing that arc are not claimed to obey the factorization.

## A directed Gauss–Bonnet control becomes a literal cell

The [square-switch surface calculation](SQUARE_SWITCH_SURFACES.md) distinguishes disjoint local rewiring attempts by interlacement. For a directed Hamiltonian path, close it by an untouched formal return; swapping two successors on the directed route always cuts it into two circuits. A pair of interlaced successor-swaps may reconnect it. This test is not a count of all paths on an arbitrary tournament.

Here there is a stricter physical cell where the two available fixed-port paths are **exactly** the two choices in the signed calculation. Start with the path \(0,1,\dots,7\), add arcs \(0\to5,\,4\to1,\,2\to7,\,6\to3\), and orient all unspecified arcs from larger to smaller. This gives a strong tournament with rows \((34,4,137,19,39,78,159,59)\). The two switched pairs are the cuts \((0,4)\) and \((2,6)\); each alone splits, together the path is \(0,5,6,3,4,1,2,7\). The direct replay enumerates every path with endpoints \(0,7\), obtaining just those two. The gate product therefore makes an exact \(2^m\) marked-path family at arbitrary depth, not an unverified assumption that no additional completions appear across blocks.

## One total, different factors

The same \(4+7=5+6\) split also makes exact products of marked path counts. Here \(6=1\cdot6=2\cdot3\), with buffers making both assembled tournaments strong at order thirteen.

| core order | rows | (entry, exit) | fixed-port count |
|---:|---|---|---:|
|4| \(2,4,9,3\) | \((0,3)\) | 1 |
|7| \(2,4,25,19,35,79,31\) | \((6,0)\) | 6 |
|5| \(2,8,3,21,7\) | \((4,0)\) | 2 |
|6| \(2,4,25,19,35,15\) | \((5,0)\) | 3 |

The unmarked Hamiltonian path totals are **10,641** and **16,883**. The marked equality really is specific to what is observed. There is no inference that 6 is an unrestricted Hamiltonian path count, nor that this gate product alone proves Chan–Pak.
