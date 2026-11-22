# A reflected rank gives a second strong lower witness

**PROVED.** For every odd \(p\ge5\), the strong three-port realization of \(E_p-1\) has a nonisomorphic strong neighbor of the same height. Reverse the arc from singleton 7 to the middle of port 6. More precisely, write \(A=2^{p-1}\) and give this port ranks \(0,\ldots,p-1\) in its transitive order. Reversing its rank-\(r\) contact from vertex 7 changes the count by
\[
            (A-1)\big(2^r-2^{p-1-r}\big),\qquad 0<r<p-1.
\]
The middle is the unique neutral contact in that entire port. Its two adjacent ranks have opposite changes of absolute value \(3(A-1)2^{(p-3)/2}\).

**FINITE-EXACT.** At \(p=61\) the connected component accessible from either witness by count-preserving single-arc flips consists of exactly these two isomorphism classes. No member has an arc leading to \(E_{61}+1\).

## The signed derivative under rank reflection

Use source rows \((2,8,3,37,15,151,63,95)\) and inflate ports \((6,4,1)\) by transitive blocks of respective widths \((m,u,v)\), with \(m\ge2,\ u,v\ge1\). In the perfect family \(m=p,\ u=v=(p-3)/2\). [Strong lower rail](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md)

Only paths using the reversed arc change. Contract its two endpoints to an ordered atom and take the count for the new order minus the count for the old. Before or after this atom a run from port 6 may be required to use ranks strictly **below** the selected rank, or strictly **above** it, respectively. In each contraction only one of these constraints can occur. Every return to a given port away from the atom requires a return edge of the original source. At most one return interval can cross the atom, and at most one can avoid it by the source's shared-endpoint criterion. Thus three runs per replicated port suffice for this finite enumeration.

Let \(y=2^u-2,\ z=2^v-2\), the counts of ways to fill two runs in the other two ports. In the difference of the ordered-atom skeletons all nonzero entries have either one or two runs in port 6 and one rank constraint. The following two columns determine all sixteen signed entries. For each row the coefficients for (one run above, one below, two above, two below) are \((-2b,d,-b,b)\).

| Runs in ports \(4,1\) | \(b\) | \(d\) |
|---|---:|---:|
| \(1,1\) | 15 | 10 |
| \(1,2\) | 11 | 8 |
| \(2,1\) | 5 | 5 |
| \(2,2\) | 4 | 4 |

For an interior rank, a single run cannot put all \(m-1\) remaining ranks on its required side. With two runs and \(L\) eligible ranks in the distinguished run, the number is \(2^L-1\), since a rank on the other side guarantees the second run is nonempty. Set \(B=15+11z+5y+4yz,\ D=10+8z+5y+4yz\). Summing the signed coefficients gives the full rank profile
\[
 \Delta_0=-2^{m-1}B,\qquad
 \Delta_r=B(2^r-2^{m-1-r})\ (0<r<m-1),\qquad
 \Delta_{m-1}=(2^{m-1}-2)B+D.
\]
At an extreme, the eligible side containing every remaining rank fills two runs in \(2^{m-1}-2\) ways and one run in one way. The other side cannot fill its distinguished run. This proves the extreme formulas and excludes a neutral extreme. The interior antisymmetry persists even when \(u\ne v\). On the perfect family put \(S=2^{(p-3)/2}\). Then \(B=4S^2-1=A-1\), giving the stated derivative. The table follows from the finite recurrence appending each permissible next source symbol or the atom, with each ordinary vertex used once and the port multiplicities capped at three. [Standalone recurrence and direct small checks](middle_reflection.py)

The midpoint reversal preserves strongness because the original arc \(7\to r\) can be replaced by \(7\to 0_{\rm port}\to r\). Before the reversal, the outdegrees of vertex 7 and this middle vertex are \(2p\) and \((3p+1)/2\); reversing the arc moves one unit from the former to the latter. Their difference is at least two for \(p\ge5\), so the two degree multisets differ.

## The entire component at exponent 61

In source-block order on vertices \(0,\ldots,123\), the midpoint arc has labels \((92,123)\). Examining all \(\binom{124}{2}=7626\) reversals in **each** of the two states gives the following signed count changes. The last two columns are the closest strictly negative and positive values to zero.

| State | Neutral arcs | Changes of \(+2\) | Closest negative | Closest positive |
|---|---:|---:|---:|---:|
| Original | 118 | 0 | −18,889,465,896,294,208,798,720 | 47,223,664,846,288,638,115,840 |
| Middle contact reversed | 116 | 0 | −618,970,018,201,538,258,838,487,040 | 1,237,940,043,897,066,293,863,383,042 |

For the original state, 116 neutral flips swap consecutive vertices inside one of its three transitive modules. For the midpoint state, 114 do so inside the four remaining modules of widths \(29,29,30,30\). In either state one additional neutral flip interchanges source twins 0 and 2 (labels 0 and 30), and the last neutral flip joins the two states. Every listed adjacent swap relabels the same tournament, and the degree multisets distinguish the two remaining states. This exhausts the reachable component up to isomorphism. It leaves open other components at the same Hamilton-path count.

The all-arc count uses a two-path decomposition. For a physical arc from vertex \(i\) to vertex \(j\), split any path using it at that arc. This gives two disjoint, nonempty directed paths covering the whole graph, the first ending at \(i\) and the second starting at \(j\). Conversely they concatenate uniquely. The gain on reversal is the opposite boundary count minus this one, since paths avoiding the edge survive. On the quotient by transitive blocks, enumerate the two paths by their ordinary vertex masks, terminal symbols and numbers of visits to each block. The ordinary masks are complementary. All block labels are allocated to their successive visits, with each visit nonempty and internally increasing. If a boundary rank is specified, its run ends there or starts there as required. Inclusion–exclusion on the unpinned runs counts the allocations exactly. In the original cyclic port quotient the shared-endpoint criterion bounds visits per path by two. When the port subgraph is acyclic, a path cannot revisit a port between consecutive ordinary vertices, so the recurrence terminates without a visit cap. [Endpoint count and rank filling](rail_arcs.py)

After each neutral reversal, consecutive twins in transitive blocks may be exchanged by relabelling. For the other neutral reversals the component replay uses canonical adjacency rows. It refines vertex cells by their counts of outgoing neighbors in each cell, then minimizes the encoding over all permutations within surviving cells. This refinement is invariant under relabelling, and minimization retains the entire oriented graph. A breadth-first search on these encodings therefore explores every reachable isomorphism class.

The same component description holds **FINITE-EXACT for every odd \(7\le p\le127\)**. An exhaustive contracted-arc calculation covered 1,422,520 signed flips across 122 host states. The original has \(2p-4\) neutral flips and the middle-reversed state has \(2p-6\). The known relabelings and the transition to the partner use all those neutral flips. There are no \(+2\) exits. At \(p=5\) the connected component has four isomorphism classes forming a square, with respectively \(7,7,5,5\) neutral arcs and no \(+2\) exits.

For the rank formula run `python3 -B middle_reflection.py`. For the whole component run `python3 -B rail_plateau.py --exponent 61`. [All-arc component replay](rail_plateau.py)

To reproduce the finite progression with the component implementation, use the following from this directory.

~~~sh
for p in $(seq 5 2 127)
do
    python3 -B rail_plateau.py --exponent "$p"
done
~~~
