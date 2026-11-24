# Two cycles on fifty shifted heights, and four ways to join scales

**PROVED with finite certificates.** For all integers \(k\ge2\), all 50 integers
\[
             h_{v,s}=4^k v-2+s,\qquad 1\le v\le25,\quad s\in\{-1,+1\},
\]
are exact Hamilton-path counts of strongly connected tournaments. On these heights put an edge when \((h+2)+(h'+2)\) is a square. The resulting graph has exactly **two** undirected Hamilton cycles at every such scale. At \(k=1\), the labels corresponding to the two genuine Chan–Pak holes are \(7+2=9\) and \(21+2=23\). Deleting them leaves a graph with three leaves and no Hamilton path.

For every \(r\ge1\) the same initial 25 vertices give a graph on \(50r\) different attained heights with at least \(2^r4^{r-1}\) Hamilton cycles under the shifted square-sum test. The construction only changes the exponent along 26 prescribed dyadic families. Each count individually has a strong realization; one can also pad every witness to one common order at most \(2k+6(r-1)+14\). The edges of this graph test its integer labels. A marked tournament arc is extra information, addressed below for several actual unit flips.

## The exact integer carrier

Write \(Q_n\) for the square-sum graph on \(1,\ldots,n\). The preceding [onset and port calculation](../../../session211/square-sum-port-bridge/PROOF.md) gives zero spanning paths at 24 and ten at 25, all with endpoint 18. Upon insertion of 25 the two individual ports are 11 and 24. Retaining their spanning two-path covers explains the ten paths despite the zero at size 24. The [two-sheet lift at 36](../square36_signed_collatz_20261123/LIFT.md) preserves this sort of information as square edges of numerical heights. To apply the two-unit shift literally, give vertex \(v\) the two odd labels \(h+2=4^kv+s\). The height which must be realized now has offset \(-3\) or \(-1\) from the dyadic center.

Use the [one-return path identity](../chanpak_absolute_liouville_leech_20261122/LEECH-LABELS.md), previously compiled at the eight pure-A Niemeier labels. Given a tournament \(D\), contact set \(I\) and a transitive block \(X\) of width \(w\ge1\), add a vertex \(b\) with \(b\to X\to D\), \(I\to b\), \(b\to D\setminus I\). Count ordered spanning two-path covers \((P,S)\) of \(D\) with \(P\) ending in \(I\), permitting an empty \(S\). Denote their number by \(G\) and the number also having \(S\) start in \(I\) by \(J\). The ordinary height of the assembled tournament is
\[
                   G2^w+H(D)-J .
\]
Each vertex of every selected \(D\) can reach \(I\), making the assembled host strong. The finite table gives the slopes \(G\) of actual marked cores for each constant term. The [replay](square25.py) records all orientations and contact sets, verifies the path-cover coefficients by subset recurrence, and counts the full hosts separately at width two.

| Odd part \(a\) | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(G\) at offset −3 |128|768|10240|448|9|88|1664|15|17|19|336|368|200|
| \(G\) at offset −1 |1024|6|5|224|9|176|26|15|34|608|336|46|1600|

Each entry is its column label times a power of two at most \(2^{11}\). Given \(v\), choose the column indexed by its odd part and the row indexed by \(s-2\). Read the width from the target by \(2^w=4^kv/G\). For \(k\ge6\) this gives a positive width throughout. For \(2\le k\le5\), the same rule works whenever the width is positive. The remaining 38 distinct integers have literal strong witnesses of order at most ten in the replay. In particular the restriction \(k\ge2\) is sharp among positive integral exponents. For \(k\ge6\) the rule also realizes every \(4^kv-2\pm1\) for arbitrary positive \(v\) with odd part at most 25, extending the finite label interval to an infinite set of coefficients.

## How the holes destroy the first graph

The two labels on the same sheet have sum \(2\bmod 4\). Opposite sheets satisfy the shifted square test just when \(v+w\) is a square, since their sum is \(4^k(v+w)\). The copies of the same \(v\) are linked just at \(v=2,8,18\), where \(2v\) is square. Thus the full graph is independent of \(k\). The [exhaustive graph replay](square25.py) counts unoriented Hamilton cycles by fixing one start and orienting by its two neighbors. For sheet sizes 23, 24, 25, 26 the totals are 1, 0, 2, 4. For size 25 both cycles are fixed by sheet exchange; their same-label edges are, respectively, \(\{2,18\}\) and \(\{8,18\}\). These result from the two ordinary spanning paths at 25 ending at 18 whose other endpoint is 2 or 8. The exhaustive count rules out additional asymmetric cycles here.

At \(k=1\) the heights are exactly the first 50 positive odds. By the [all-order hole theorem and the small-height atlas, THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md), precisely two heights in this interval are absent. They are \((v,s)=(2,+1),(6,-1)\). In the graph of the other 48, heights 39, 75 and 89 have degree one (their base sheet addresses are \((10,+1),(19,+1),(23,-1)\)). This proves failure of a spanning path at the first scale. At each larger scale both addresses have strong witnesses and the two cycles return. The missing absolute heights explain this local graph defect; square-sum rules alone do not exclude 7 and 21 as tournament counts.

## Three one-arc selectors

For odd parts 9, 15 and 23, the two offsets can also come from reversing a single fixed arc of the core with the same contact set. The following core rows are outgoing masks; vertex indices start at zero.

| Odd part | Common slope | Contact | Marked pair | Core rows |
|---:|---:|---:|---|---|
|9|9|9|0,2|4,1,2,7|
|15|15|17|0,3|10,4,1,6,15|
|23|736|130|4,5|4,5,24,131,107,15,47,119|

In each row the old and flipped profiles have the stated common slope and constants −3 and −1. Both cores can reach their contact set. All other orientations, including the growing block, stay fixed, giving strong tournaments related by this one arc for every positive width. In particular at the base vertex 18, a square crossing used in **both** of the fifty-vertex cycles, the two heights differ by a physical arc reversal for every scale in the theorem. Expressing the same relation at every numerical square edge would require additional marked hosts.

## Joining copies by decompositions of 25

Take blocks indexed by \(j=0,\ldots,r-1\) and replace \(k\) in block \(j\) by \(k+3j\). Consecutive center scales have ratio 64. The blocks are disjoint since each earlier maximum coefficient is 25. In either local cycle the following edges occur with **both** assignments of the two opposite sheet signs:
\[
              16+9=13+12=19+6=25 .
\]
For any consecutive blocks, remove the edge with bases \(16,9\) from the lower cycle and choose either \(13,12\) or \(19,6\) in the upper one. If 16 carries sign \(s\), give 13 or 19 the same sign. In place of the deleted edges use the two cross-block edges between opposite signs. Their base sums are
\[
\begin{aligned}
16+64\cdot12&=28^2, &9+64\cdot13&=29^2,\\
16+64\cdot6&=20^2, &9+64\cdot19&=35^2.
\end{aligned}
\]
The common lower scale is a square. These two cross edges join the two cycles into one. The incoming and outgoing removed edges of any interior block are vertex-disjoint, so the joins can be performed successively. There are two choices of the local cycle at each block and four choices at each junction (the sign and the upper decomposition). The cross edges identify which two edges were removed, so an unoriented output cycle recovers all of these choices. This proves the claimed lower bound on distinct cycles. Additional edges on the full induced graph are allowed.

On the earlier 36-vertex lift the square of an attained multiplier preserves **unshifted** square sums by ordinal multiplication. Here the assembly preserves the actual offset of two at every scale, and the two good cycles retain identical splice edges. The finite initial signature needed by this construction consists of the two path-cover readers for each odd part and the placement of three marked edges of square label 25 in each lifted cycle. It gives an inverse for these exact heights. Exhaustion of all positive odd integers by tournament heights remains open.

## One hundred exact heights cover all odd residues modulo 24 on two separated copies

**PROVED.** For every \(k\ge6\) change the sheet displacement from 1 to 3 on a second copy of the construction. Its shifted labels are \(4^kv\pm3\), hence its required heights are \(4^kv+1\) and \(4^kv-5\). The replay certifies 26 more one-return profiles (same thirteen odd parts and the two new constants), each with \(G/a\le2^{11}\). Hence the second copy also consists entirely of exact strong heights. It again has exactly two Hamilton cycles, since cross-sheet displacements cancel and same-sheet shifted sums are 2 modulo four.

At every such scale the four offsets \(-5,-3,-1,1\) modulo eight give \(3,5,7,1\) respectively. Independently the base label runs through all three residues modulo three, with \(4^k\equiv1\pmod3\). Thus the 100 heights meet **all twelve odd residues modulo 24**, six in each of the two copies. The [eleven-vertex threshold](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md) applies when the precedence tournament is transitive; these 100 heights have explicitly selected strong hosts without a prescribed precedence order.

In fact the two copies have no edges between them under the shifted square-sum test. A mixed pair would give \(4^k(v+w)+c\) square, for \(c\in\{-4,-2,2,4\}\). The two middle values are impossible modulo four. Dividing by four at \(c=-4\) gives 3 modulo four. At \(c=4\) put \(t^2=4^{k-1}(v+w)+1\). Of the consecutive even numbers \(t-1,t+1\), one has exactly one factor of 2, so the other is divisible by \(2^{2k-3}\). Hence \(t\ge2^{2k-3}-1\), whereas \(t^2\le50\cdot4^{k-1}+1\). These inequalities conflict already at \(k=5\) and diverge thereafter. Each copy is connected by its Hamilton cycles, so the induced graph has **exactly two components**. This is a finite residue-complete bank at every scale in which geometric communication between two halves fails for an explicit integer reason.

The separation extends to the \(r\)-block construction whenever \(k\ge\max(6,5+3(r-1))\). In the same argument replace 50 by \(50\cdot64^{r-1}\). The block-splicing argument applies unchanged in each component with displacement 1 or 3. Thus at these scales there are two connected components on \(50r\) heights apiece, each with at least \(2^r4^{r-1}\) Hamilton cycles.
