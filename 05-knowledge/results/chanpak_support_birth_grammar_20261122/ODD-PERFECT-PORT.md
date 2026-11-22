# A strong all-odd-exponent lower rail, and why one apex cannot straddle

**PROVED.** For every odd \(p\ge5\), a fixed eight-vertex source with three growing transitive slots has a strong realization of \(E_p-1\), where \(E_p=2^{p-1}(2^p-1)\). The order is \(2p+2\). Unlike the known ordinal-product lower witness, the construction is strongly connected and has one fixed compressed-path grammar. It does not by itself give the upper neighbor or an arc joining the two counts. The one-apex obstruction below holds at every exponent \(p\ge3\) for that architecture.

## A three-port strong grammar for the lower endpoint

Use outgoing masks \((2,8,3,37,15,151,63,95)\) on vertices \(0,\ldots,7\); choose ports \((6,4,1)\) in that order. Replace them by transitive tournaments of widths
\[
                 p,\quad (p-3)/2,\quad(p-3)/2.
\]
Every edge between different blocks follows its edge in the source. The source is strong and so is every inflation: reachability between blocks follows source paths, and reachability within each block follows a source cycle through its port. With five unmarked vertices the total order is \(2p+2\).

**The bounded-run reduction.** For any port \(z\), each return between consecutive appearances in a spanning path requires an arc from an outneighbor of the port to an inneighbor. Different returns use disjoint such edges in the inflated graph. Check the return edges in the eight-vertex source, treating endpoints at other marked ports as freely replicable and ordinary endpoints as singletons. For each port, any two potential returns have at least one shared ordinary endpoint. Thus no mark can occur in more than two runs at any inflation size. The full small return test and exact table are in the standalone replay. This uses the [two-run criterion](../chanpak_passive_pubsub_20261117/reader/PROOF.md), here with the other marked slots inflated too.

Collapse each maximal run in a marked block to its mark. Unmarked vertices occur once; a mark appearing once has a unique filling, while a mark appearing twice has \(X_i-2\) fillings for \(X_i=2^{\text{width}}\). Indeed, choose a nonempty proper subset for the first run and put the complement in the second, each in increasing order. The finite compressed path grammar is recursion on the last symbol and the remaining unmarked mask, with each port multiplicity capped at two. It gives:

| Run counts at \((6,4,1)\) | Skeletons |
|---|---:|
| \(111,112,121,122\) | \(45,33,15,12\) |
| \(211,212,221,222\) | \(15,11,5,4\) |

Write the three slot powers as \((X,Y,Z)\). Expanding the table, the exact count is
\[
H=-1-X-3Y-3XY+3Z+3XZ+4YZ+4XYZ.
\]
On the selected widths \(Z=Y\), the linear and mixed antisymmetric terms cancel. Also \(X=8Y^2\), so
\[
H=4XY^2+4Y^2-X-1
 =32Y^4-4Y^2-1
 =2^{p-1}(2^p-1)-1.
\]
The new ingredient is the pair of independently splittable half-width slots: their balanced diagonal cancels the intermediate modes but keeps the \(2^{2p-1}\) mode. This is a scalar lower rail, not a claim that the source can be completed or flipped to make the upper endpoint. For Mersenne exponent \(2\) an all-exponent straddle is impossible outright, since the upper target is the excluded value \(7\). The unrestricted question starts at exponent \(3\).

## Excluding the fixed one-apex single-contact architecture

Let a strong tournament consist of a transitive spine and one apex. A strong-to-strong reversal of an apex contact gives two normalized words \(U0V,U1V\), where both flanks are nonempty. Let \(u\) be the binary value of \(U\), \(w\) the value of the reverse complement of \(V\), and \(F=1+p\) the one-apex binary count. The [single-contact identity](../chanpak_apex_flip_local_global_20261122/README.md) says
\[
L=2uw+u+[u\text{ odd}]+F(u)+F(w)-1,\quad
L'-L=2(\lceil w/2\rceil-\lceil u/2\rceil).
\]
Negative orientation is reduced to positive by reflecting and complementing. For unit gap, \(1\le w-u\le3\).

Set \(T=2^{p-1}\), so that \(E_p-1=2T^2-T-1\). If \(u\ge T\) the product already exceeds the target. If \(u\le T-3\) then \(w\le T\); the bound \(1\le F(i)\le i+1\) gives \(L\le2T^2-3T-4\), below the target. Thus only \(u=T-2\) and \(u=T-1\) need be checked. The exact ceiling equation now leaves four pairs, and the binary recurrence gives:

| \((u,w)\) | \(L\) |
|---|---:|
| \((T-2,T-1)\) | \(2T^2-9T/2+3\) |
| \((T-2,T)\) | \(2T^2-3T/2-1\) |
| \((T-1,T+1)\) | \(2T^2+3T/2-1\) |
| \((T-1,T+2)\) | \(2T^2+15T/4-1\) |

The displayed formulas apply for \(T\ge8\), where none equals \(2T^2-T-1\). At \(T=4\) the four exact lowers are \(17,25,37,45\), again missing the target 27. Thus no strong perfect straddle can use only one apex contact in a fixed transitive deletion, regardless of the number of bits in the flanks or their apparent residue flexibility.

There is no escape by reversing an ordinary arc while keeping the **same** apex deletion transitive in both hosts. A single transitive-preserving reversal must be between adjacent spine vertices, so the tag word either stays fixed or swaps adjacent \(01\) to \(10\). The weighted inversion formula has pair weight one at distance one, and \(2^{d-2}\) at distance \(d\ge2\), for the half-count \((H-1)/2\). The swap adds its own one and moves each earlier 1 and each later 0 farther from its partner. The distance weights are nondecreasing and strictly increase beyond distance two. Both words must start in 1 and end in 0, so a unit gap forces the swapped pair to follow the first vertex and precede the last. The only possibility is \(1010\leftrightarrow1100\), with counts \(9,11\). A fixed transitive deletion therefore cannot support any of the desired strong perfect seams by any single arc.

The companion Python replay gives the exact skeleton table and independent small inflated literal counts, together with the four one-apex branches. A common strong one-arc grammar for both signs at every odd exponent remains open here.
