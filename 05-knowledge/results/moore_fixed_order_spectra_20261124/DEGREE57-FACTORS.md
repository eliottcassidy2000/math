# Many perfect matchings multiply native degree-57 Hamilton paths

**Conditional on the existence of a degree-57 Moore graph.** Use the native 3250-site tournament \(T\) and its marked tile mate \(T^\ast\) as in [the original bounds](DEGREE57-HEIGHT.md). Every complement spine satisfies
\[
 2^{14437}<H(T),H(T^\ast),\qquad H(T)<2^{28860},\quad H(T^\ast)<2^{28863},
\]
and **one can choose a common spine** on which both counts exceed \(2^{18706}\). The first improvement is uniform; the second is not a lower bound for arbitrary spines. All the counts are odd. Thus no native endpoint at any spine hits \(E_p\pm1\) for \(2\le p\le7219\), where \(E_p=2^{p-1}(2^p-1)\). The chosen high-count spine misses these targets through \(p=9353\). These are exclusions, not a claim of graph existence, arithmetic congruences, or which exact heights occur.

## 1. An integer refinement of spectral matching

Put \(n=3250,\ k=1625,\ A=A(G)\). The Moore identity gives
\[
 B=2A+I,\qquad B^2=225I+4J,\qquad B\mathbf1=115\mathbf1.
\]
If the crossing graph for a balanced cut has matching size at most 1577, extend a König cover to size exactly 1577. Its complement leaves cross-anticomplete \(X,Y\), of sizes \(x\le y\le1625\) and \(x+y=1673\). The centered spectrum of \(B\) gives \(520xy\le9\cdot3250\cdot1577\); hence \(48\le x\le54\). For \(g=B\mathbf1_X\), the entries on \(Y\) vanish, those on \(X\) are odd, and the remaining 1577 entries are even. Yet
\[
 \mathbf1^{\mathsf T}g=115x,\qquad \|g\|^2=225x+4x^2.
\]
Start the odd sites at 1 and the even sites at 0. To attain the sum, make \(57x\) increments of 2, greedily at the cheapest squared-norm cost. The marginal costs per even site are \(4,12,20,\ldots\); per odd site they are \(8,16,24,\ldots\). Each site's costs are increasing, so the global list of cheapest increments is feasible and exact. The excesses of the minimum over the required square norm are \(56,124,184,236,280,316,344\) in increasing order of \(x\). Thus every balanced cut has a matching of size at least 1578.

## 2. Count compatible block orders by switching

Use exactly 1578 matched backward pairs as directed blocks in any spine, together with fixed Hamilton prefix/suffix paths on the unmatched sites from the respective halves. For the mate, count in its converse. Put a dummy anchor in the block order: its outgoing and incoming transitions represent joining the prefix and suffix. Of the \(N=1579\) objects, each has at most \(D=59\) forbidden distinct successors and predecessors. Without the anchor this follows from Moore degree 57 plus at most one additional cross-spine pair in the mate converse; the anchor costs at most one more per block on each side. The anchor itself has at most 58 on each side. There are at most \(E=52929\) forbidden ordered pairs: at most \(\lfloor65\cdot3250/4\rfloor=52812\) cross edges by the cut spectrum, one cross-spine pair in the mate converse, and at most 116 anchor exclusions.

Here is an elementary factorial counting lemma. For cyclic orders of \(N\) objects modulo rotation, with \(E\) forbidden oriented adjacencies and maximum forbidden indegree and outdegree \(D\), the number avoiding all forbidden adjacencies is at least
\[
  (N-1)!\left(\frac b{b+1}\right)^E,\qquad b=N-4D-3>0.
\]
Add the forbidden adjacencies one at a time. Among cycles passing all previous restrictions but containing the newly forbidden \(u\to v\), swap \(v\) with an object \(w\) other than \(v\) and its two neighbors. Write \(u,v,z\) and \(a,w,c\) for their two local triples. The four new adjacencies are \(u\to w,\ w\to z,\ a\to v,\ v\to c\). Each rules out at most \(D\) choices of \(w\), even when tested against the full forbidden list: the maps sending \(w\) to its predecessor or successor in the old order are bijections. So there are at least \(b\) switches to passing cycles. Given a passing cycle in the image, the old cycle and the choice of \(w\) are recovered by swapping \(v\) with the current successor of \(u\). Thus the number passing is at least \(b\) times the number failing at this stage, yielding the product bound. Cutting each resulting cycle at its anchor and expanding injectively gives actual Hamilton paths. Exact integer comparison gives the uniform exponent 14437.

## 3. Choose a spine with enough room for a 13-factor

There is a perfect matching in \(G\). Every nontrivial edge cut has at least 57 edges: singleton cuts attain it, and for \(2\le s\le1625\) the centered spectrum gives the strict lower bound \(50s(n-s)/n>57\). After deleting a nonempty set \(S\), each odd component sends at least 57 edges to \(S\); there are at most \(57|S|\) such edges, so Tutte's criterion applies. For \(S=\varnothing\) use connectedness and even order.

Split a fixed perfect matching across a balanced cut with independent fair choices of which end goes to the left. At any vertex the cross-degree is \(1+\mathrm{Bin}(56,\tfrac12)\): two other neighbors cannot be in one matched pair because that would form a triangle. Each pair of matched blocks supports at most one nonmatching edge. Two such edges sharing a vertex would form a triangle; two disjoint ones would form a quadrilateral with the matching. Consequently the nonmatching crossing indicators, each the exclusive-or of two independent cut bits up to a constant, are pairwise independent. The total cross-edge count has mean 47125 and variance 22750. The union bound for minimum degree below 16, plus Chebyshev for more than 47625 edges, is
\[
 \frac{3250}{2^{56}}\sum_{j=0}^{14}\binom{56}{j}+\frac{22750}{500^2}<1.
\]
Select a cut with both desired bounds. Each half of the complement has minimum degree at least 1567 and hence a Hamilton cycle by Dirac. Open one cycle in each half, choosing the first vertex in the right half to avoid a graph edge from the final vertex of the left half. This produces a complement spine.

The crossing graph on this cut has a 13-regular spanning subgraph. The max-flow cut test requires \(e(X,Y)\ge13(s+t-k)\) for all opposite sets of sizes \(s,t\) with \(s+t>k\). It follows by assigning capacities 13 from the source to all left sites and from all right sites to the sink, and capacity 1 to each left-to-right graph edge. Assume \(s\le t\) by swapping sides. Put \(a=s+t-k,\ u=k-t\). Two lower bounds suffice:

- If \(W=16s-13a\) satisfies \(2W\ge u\) and \(W(W-u)\ge us(s-1)\), the \(u\) sites outside \(Y\) receive at most \(W\) edges from \(X\). Indeed, each pair in \(X\) shares at most one neighbor. If \(h\) edges go to those \(u>0\) sites, Cauchy–Schwarz applied to their degrees gives \(h^2-uh\le us(s-1)\). For \(u=0\) the claim is immediate. Minimum cross-degree 16 then gives \(e(X,Y)\ge13a\).
- Otherwise the centered spectrum of \(B\) gives \(2n\,e(X,Y)\ge115st-15\sqrt{st(n-s)(n-t)}\). This supplies the required bound if \(Z=115st-2n\cdot13a\ge0\) and \(Z^2\ge225st(n-s)(n-t)\).

The integer replay verifies the first test for 32955 size pairs and the second for the remaining 628014. No candidate Moore graph is used by the replay.

## 4. Count both the covers and their orders

By [Schrijver's perfect-matching bound](https://ir.cwi.nl/pub/2154/2154D.pdf), the 13-factor has at least \((12^{12}/13^{11})^{1625}\) perfect matchings. For **each** such matching, order the resulting 1625 backward blocks using the switching lemma. This time the anchor is unrestricted; \(D=58\) and \(E\le47626\), allowing for the extra mate-converse cross-spine pair. All these paths alternate right and left, starting on the right. A path uniquely reads back both the matching (its right-to-left steps) and the order of the blocks. Thus the multiplication counts disjoint paths, not multiple descriptions. Both endpoints of this one marked tiling line satisfy
\[
 H\ \ge\
 \left(\frac{12^{12}}{13^{11}}\right)^{1625}
 1625!\left(\frac{1391}{1392}\right)^{47626}>2^{18706}.
\]

The distinction between a single matching and the space of matchings is the source of the larger gain. A prospective search for lower native heights now has two separate obstructions to monitor: near-factorial ordering once a matching is fixed, and exponential multiplicity of matchings when the cut itself is sufficiently regular. Neither controls exact odd residues.

Reproduce with `python3 degree57-factors.py`. The square-moment refinement and factorial-scale idea came from a follow-up research pass; the switch argument and many-matching factor bound give the stronger exponents here.
