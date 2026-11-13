# Radial growth in the single-arc-flip metagraph

Here the edges mean **one arc flip**, not complement lines. Merge tournament isomorphism classes with their converses. The distance in this quotient is the minimum labelled Hamming distance over all identifications and a possible converse: lift a quotient walk one edge at a time, and project a shortest labelled walk for the other inequality.

## A stable embedding, and its limit

For $k\ge0$ set $P_k(T)=TT_k\to T\to TT_k$. This commutes with converse up to isomorphism. It is an isometric embedding, and its image is **geodesically convex**: every shortest path between two image vertices stays in the image.

To prove it, take any identification between two padded tournaments and let $s,t$ be their outermost sources. If $s$ is not mapped to $t$, swap its image with the image of the preimage $x$ of $t$. The mismatch on $sx$ disappears. For every other vertex the two affected arcs compare binary pairs of the form $(1,b)$ and $(c,1)$ before, versus $(1,b)$ and $(1,c)$ after. The distance does not rise because $|b-c|\le (1-b)+(1-c)$. Thus every minimizing identification pairs the sources. Delete them and repeat, then apply the sink version of the same argument. This works in both possible converse frames. Hence a minimizing identification aligns the entire padding. A shortest walk in the labelled cube cannot flip any edge on which its endpoint labels agree. Lifting any shortest quotient walk therefore shows that it stays padded, proving convexity as well as isometry.

In particular adding four sources and four sinks is an exact $n\mapsto n+8$ stabilization for this marked metric, preserving the Hamiltonian count as well. This is **not** Bott periodicity, nor does it say the whole larger graph repeats. It does not even preserve complement-line wiring: if $k>0$, the entire image is nonstrong, and [THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) makes nonstrong nodes independent for complement lines. The analogy with [marked Bott work](../collatz_bott_marked_periodicity_20261007d.md) is the difference between an exact transported operation and an invariant that forgets the marked interface. See also [Bott's original stable homotopy paper](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/bott4.pdf) for the genuinely topological statement.

## An exact first radial shell and a full radial product law

A neighbor of the transitive tournament is obtained by reversing an arc across a gap spanning $m\ge3$ positions in a total order. The interval becomes the strong tournament $S_m$ with just its first-last arc reversed; all earlier and later vertices are singleton SCCs. The $m=2$ flip merely exchanges two adjacent entries of a total order and gives no new node. The interval is strong because its forward path and reversed endpoint arc give a directed cycle. It is self-converse by reversing the interval order. Thus the neighbors are parametrized exactly by $m$ and the unordered pair of outside lengths $\{a,b\}$ with $a+b=n-m$. Their number and ordinary generating function are
$$
\left\lfloor (n-1)^2/4\right\rfloor,\qquad
\frac{z^3}{(1-z)^2(1-z^2)} .
$$
This is a wiggly support degree, not the complement-line degree. It gives the exact values 16 and 81 at orders 9 and 19.

More generally, let $S(z,u)$ count ordinary strong nonempty atoms by order and distance from transitive, and $S_{\rm sc}$ count their self-converse subset. Distance from transitive equals the minimum backward-arc count over total orders. It is additive across SCCs: any total order pays at least each internal minimum, and concatenating optimal orders in SCC order pays no cross cost. Therefore the ordinary class series is $(1-S(z,u))^{-1}$. A self-converse SCC word consists of an arbitrary half-word, an empty or self-converse strong center, and the reversed converse half. Its series is $(1+S_{\rm sc}(z,u))/(1-S(z^2,u^2))$. Averaging those two expressions gives the merged radial generating function. It is exact in terms of the strong atom series; it does not claim those atom series are closed-form or recover edges between radial layers.
