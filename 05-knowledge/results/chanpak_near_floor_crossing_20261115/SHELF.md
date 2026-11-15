# Crossing can have negative defect at unbounded strong order

For a tournament \(T\) and a subset \(X\), let \(T+x\) add one vertex beating exactly \(X\). Write \(H\) for the Hamilton path count, \(\delta=H(T+x)-2H(T)\), and call \(X\) crossing if some Hamilton path uses at least two \(X\to X^c\) arcs. [HYP-9732](../../hypotheses/HYP-9732-crossing-lemma-a-double-crossing-cut-has-defect-at-least-three-so-unit-cuts-are-monotone.md) predicted \(\delta\ge3\) for crossing cuts. The accompanying 11-vertex witness already refutes even the unit-cut consequence. This note gives an unbounded, strong family and an exact capacity reformulation.

## One growing transitive block

Take the 12-vertex outgoing masks (LSB vertex zero)
\[
(4030,4024,4090,1008,4064,4032,3971,3840,3584,3072,2056,8),
\qquad X=\{3,5,7,8,9,10,11\}.
\]
Replace the transitive module on vertices \(7,8,9\) by a transitive block \(B\) of size \(r\ge1\), all marked. Retain all the orientations outside and to that block. The outside types are early \(0,\ldots,6\) and late \(10,11\). Early vertices beat \(B\), and \(B\) beats late vertices. The only late-to-early arcs are \(10\to3\) and \(11\to3\). This description is unchanged in the extension if the new vertex joins the early side: it beats the block and both late vertices.

After identifying the block with a single vertex when \(r=1\), both the source and the extension are strongly connected (direct reachability check). Replacing a vertex of a strong tournament by a nonempty transitive block preserves strong connectivity: lift any outside route, and to connect any two vertices in the block use a directed cycle through the original vertex. The source path \(0,2,1,3,4,5,6,B,10,11\), with \(B\) in ascending order, is Hamiltonian and crosses the cut at \(3\to4\) and \(5\to6\), for all \(r\).

A Hamilton path visits \(B\) in at most two contiguous runs. Each return after leaving it needs a late-to-early edge, and every such edge consumes vertex \(3\). Compress each run to a slot in the outside path. One run has exactly one expansion. Two ordered nonempty runs have exactly \(2^r-2\) expansions, since one chooses a nonempty proper subset for the first run and each subset then has its unique internal order. Thus both counts are affine in \(2^r-2\), with coefficients determined by \(r=1,2\). Literal DFS gives \((H,H^+)=(535,1101),(1069,2157)\), respectively. Hence
\[
 H(T_r)=267\cdot2^r+1,\qquad
 H(T_r+x)=528\cdot2^r+45,\qquad
 \delta(T_r,X_r)=43-6\cdot2^r.
\]
For all \(r\ge3\) this is negative, with \(\delta/H\to-2/89\). The construction has exactly three arcs oriented across the cut in the crossing direction, independent of \(r\). [Replay](verify_shelf.py) independently enumerates literal Hamilton paths for \(1\le r\le5\) and checks reachability.

## Directed cut capacities and an exact minimization rule

For each arc \(w\to u\) in the old tournament define \(J_{wu}\) as its usage among Hamilton paths of \(T\). Let \(K_{u,w}\) count ordered spanning two-path covers whose first path ends at \(u\), whose second starts at \(w\); both paths are nonempty. Give \(w\to u\) capacity \(c_{wu}=J_{wu}+K_{u,w}\).

A Hamilton path with \(j\) crossings \(X\to X^c\) contributes \(j+1\) valid insertion positions, by telescoping its changes of side and counting permissible endpoints. A new Hamilton path not from a valid insertion has a broken junction \(u,x,w\) with \(u\in X^c,w\in X,w\to u\). Removing the new vertex gives exactly a cover counted by \(K_{u,w}\). Thus
\[
 \delta(T,X)=-H(T)+\sum_{w\in X,\ u\notin X,\ w\to u}c_{wu}.
\]
All capacities are even: every cut has even capacity by Rédei parity applied to both counts, including the empty/full cuts; subtract the cut of \(\{w,u\}\) from the sum of the singleton cuts to isolate \(c_{wu}\). This is the nonnegative arc-supported version of the cut identity; it is not pairwise identification with the odd-cycle-path weights, which may also be nonzero on reversed arcs.

To minimize defect over all crossing cuts of a fixed tournament, enumerate the pairs of disjoint arcs \(w_1\to u_1,w_2\to u_2\) that co-occur in at least one Hamilton path. For each, force \(w_1,w_2\) on the source side and \(u_1,u_2\) on the sink side, and take a minimum directed cut in the network above. The least of those capacities minus \(H(T)\) is exactly the least crossing defect. Each crossing path supplies such a pair; conversely every cut produced from a co-usable pair is crossing. Co-use can be tested by contracting each required arc to a two-vertex directed word, allowing transitions precisely when the exit of one word beats the entry of the next; ordinary subset DP works on this directed transition graph. This avoids searching all cut masks once the capacities and co-use table are available. Minimum defect does not, by itself, decide existence at an exact prescribed defect such as \(+1\).
