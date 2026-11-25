# Two-coordinate rigidity and exact searches beyond congruences

The common idea is to keep two coordinates, rather than a single readout. Reversing a marked edge gives two independent directions of probes for a Moore deck. Reversing the low bits gives a monotone search coordinate for a short-guard tournament even though its ordinary binary value is not monotone.

## 1. A rectangular uniqueness criterion for the normalized degree-57 deck

This is conditional on comparing **two valid completions**, not a claim that one exists. Fix the marked root–hub edge \((r,p_0)\). Write \(p_i\) for the other heads, \(w_x\) for the satellites of \(p_0\), and \((i,x)\) for the interior satellite adjacent to \(p_i,w_x\). Re-rooting at \(p_0\) with hub head \(r\) interchanges \(p_i\) and \(w_x\)'s roles and replaces \((i,x)\) by \((x,i)\). In either chart a grid vertex is adjacent to one other grid vertex in each other branch, and none in its own. Thus the interior adjacency has degree 55 and one neighbor in each other row and column.

Freeze all incidences of any \(f\) interior rows and \(g\) interior columns, including incidences with the unfrozen lines. There is at most one valid completion for the pairs
\[
 (22,19),\quad(23,16),\quad(24,13),\quad(25,9),\quad(26,5),\quad(27,0),
\]
their transposes, and supersets of the frozen lines. The final pair is the earlier one-axis bound (the edge already fixes the hub).

Let \(R\) be the frozen \(f\times g\) rectangle and \(X\) be the active grid vertices, with sizes \(P=fg\) and \(x\). For each newly gained edge \(zv\), probe \(v\)'s own head and hub satellite, and then its one neighbor in each frozen row and column. Each probe is unchanged; as in the [frozen-probe argument](GUARDS-AND-DECK-DISTANCE.md), its unique old common neighbor with \(z\) must be deleted, and each such deletion is distinct by girth five. There are \(2+f+g-c(v)\) distinct probes, where \(c(v)=d_R(v)\). Writing \(F=\max(f,g),\ t=F+2,\ u=2+f+g\), we have
\[
 r_z\ge u-c(v)\ge t,\qquad r_z\le55-F.
\]
The upper bound holds because edited partners must be in different unfrozen rows and different unfrozen columns, excluding the site's own. In the graph of newly gained edges, girth at least five and minimum degree \(t\) force \(x\ge t^2+1\). Moreover each gained neighbor of \(v\) has gained degree at least \(u-c(v)\); their distance-two branches are disjoint. Thus \(x\ge1+t(u-c(v))\). It follows that every active vertex obeys the individual lower bound
\[
 c(v)\ge \ell_x=\max\{f+g+F-53,\ u-\lfloor(x-1)/t\rfloor\}.
\]
The centered spectrum of \(2A+I\) also bounds \(2n e(X,R)\le115xP+15\sqrt{xP(n-x)(n-P)}\), \(n=3250\); use this to discard impossible \(x\).

### Apply the exact square moment across all four rectangles

The remaining three grid regions, of sizes \(P_A=f(56-g), P_B=(56-f)g, P_C=(56-f)(56-g)\), consist respectively of row-only, column-only and entirely unfrozen sites. Let \(m=2e(R)\) be the (even) internal adjacency sum. The sums of \(d_R\) within the four regions \(R,A,B,C\) are exactly
\[
 m,\quad(f-1)P-m,\quad(g-1)P-m,\quad(57-f-g)P+m .
\]
For example each vertex of \(R\) has \(f-1\) neighbors in the other frozen rows, giving the first two sums together. The grid degree 55 then fixes the last one. By \(A^2+A=56I+J\), the global square sum of \(d_R\) is \(56P+P^2-m\). The \(f\) frozen heads each contribute \(g^2\), the \(g\) frozen hub satellites each contribute \(f^2\), and the other nongrid sites contribute zero. So the **grid** square budget is exactly \((56-f-g)P+P^2-m\).

For each integer \(x\in[t^2+1,P_C]\) surviving the first spectrum test and each possible even \(m\), balance the entries of each region as tightly as its required sum allows. Impose the floor \(\ell_x\) on \(x\) entries of region \(C\) (the active ones). If the mean in \(C\) is below \(\ell_x\), the square minimum sets these \(x\) entries to \(\ell_x\) and balances the rest, by moving units from any higher active entry to a lower free entry. If the mean is at least the floor, ordinary balanced entries suffice. Reject if there is not enough sum to meet the floor or the total square minimum exceeds the grid budget. For \(m\), use the obvious row/column degree cap and the centered eigenvalue bounds
\[
 65P^2/n-8P\le m\le50P^2/n+7P.
\]
The exact replay rejects all remaining \(x,m\) pairs for the five nontrivial pairs in the statement (490, 545, 596, 674, 578 respective \(x\)'s survived the first spectrum test). Dropping \(g\) by one in each case leaves some pairs passing these necessary inequalities; this is **not** evidence of a different valid completion. It identifies where this relaxation stops.

## 2. Physical +2 pairs reach every admissible mixed congruence

Use the [six-site fork with a one-site guard](GUARDS-AND-DECK-DISTANCE.md). Write the even payload as \(2v\), with fixed leading prefix (for instance 1010), and put
\[
 q(v)=v-(h(2v)-1)/2,\qquad
 H(v)=1662v-58q(v)+427.
\]
The actual two counts are \(H(v),H(v)+2\). The weighted inversion formula gives \(q(2v)=q(v),\ q(2v+1)=q(v)+\lceil v/2\rceil\): equivalently \(q(v)\) is the sum of same weights \((1,1,2,4,\ldots)\) on pairs of 1 bits by distance. In particular \(0\le q(v)\le v/2\) also follows from \(h(2v)=h(v)+v+(v\bmod2)\) and \(h(v)\ge1\).

For any modulus \(M\), keep the state \((v\bmod2M,q\bmod M)\). From any first coordinate, a block of \(K\) bits with \(2^K\ge2M\) can set that coordinate to any prescribed value. To control the other coordinate, first set the first to zero. Write \(M=2^a m\) with \(m\) odd; choose \(L\ge2\) with \(2^L=1\bmod m\). From first coordinate zero, append \(L\) ones then \(a+1\) zeros. The first coordinate returns to zero, while the increment in \(q\) is \(2^{L-1}-1\), a unit modulo \(M\) (it is odd and is \(-1/2\bmod m\)). Repeating this loop sets the second coordinate, with a final known correction for the block prescribing the first. This reaches every state from any prefix. With \(M=2^a\), the dyadic case can be routed more directly in \(O(a^2)\) bits: each block \(1\,0^{i+1}\,1\,0^{a+1}\), for \(0\le i<a\), preserves first coordinate zero and adds \(2^i\), so concatenate the blocks of the desired binary increment.

Each state has exactly two predecessors (last bit fixed by parity; then two old first coordinates and exactly one old counter for each). The uniform distribution is stationary, and the zero-bit self-loop gives aperiodicity. Therefore any fixed prefix has uniform limiting state law under long uniform extensions. Since \(\gcd(1662,58)=2\), this gives every residue modulo odd \(M\), and every odd one modulo even \(M\), with equal limiting frequencies **among descriptions**. This is joint mixed-modulus reach, not a claim that all sufficiently large integers in those residues occur.

## 3. A sublinear exact integer inverse for that short-guard family

One can test a prescribed actual lower count \(H\), not just a congruence. First \(1604v+427\le H(v)\le1662v+427\) puts \(v\) in an interval of ratio below 2, meeting at most two binary widths. For a width \(L\), take \(t=\max(0,\lfloor(L-7)/2\rfloor)\). Split \(v=2^t u+w\), and let \(R(w)\) be the \(t\)-bit reversal of \(w\) (including leading zeros in the suffix); \(b_1=R(w)\bmod2\) is the first suffix bit. For \(t\ge1\), the concatenation formula is
\[
 q(v)=q(u)+q(w)+\lfloor u/2\rfloor R(w)
              +(u\bmod2)\bigl(R(w)+b_1\bigr)/2.
\]
The cross term follows from summing \(\lceil\text{prefix}/2\rceil\) at the 1 bits in the suffix. For fixed prefix \(u\), the final term is nondecreasing in \(R(w)\), while \(\lfloor u/2\rfloor R(w)\) increases by at least \(\lfloor u/2\rfloor\) per step. The possible increase in \(1662w-58q(w)\) is at most \(1691(2^t-1)\) even across the entire suffix interval. In our split \(u\ge64\cdot2^t\), so \(58\lfloor u/2\rfloor\) exceeds that maximum; hence \(H(v)\) **strictly decreases in reversed-suffix order** at fixed \(u\). When \(t=0\) there is just one candidate.

Enumerate the \(O(\sqrt H)\) possible high prefixes, and binary search the low bits in reversed order in each. Check actual equality and the sentinel/payload restrictions. The cost is \(O(\sqrt H\,\log^2 H)\) small arithmetic/bit scans, allowing extra factors for large-integer arithmetic. It is complete in this chart and can return multiple preimages: for 546598033 it gives both even payloads 659056 and 666074. The remaining barrier is proving a suitable chart always returns a nonempty answer, not finding a representation if one is hidden among the short-guard candidates.

Two scope checks for further work:

- Gap control alone is cheaper: [the one-apex inverse](FIBRE-TRANSFER.md) with empty flank and other flank \(2d-1\) gives actual counts \(h(d-1)\) and \(h(d-1)+2d\le3d\), at logarithmic order. The [six-site determinant](ENDPOINT-AND-SIGNED-DECK.md) adds fixed transfer guards, not the first small realization of bare gaps.
- A fixed strong one-apex flank has at most **eight** contact choices for any selected start or end sum, despite realizing every possible sum. For starts, read the tail weights \(p_i\), \(i\ge2\), backwards as in the determinant proof. Each is the old total or the old total plus one. In the equality case the only overlap of the subset sums is at the touching endpoints, both uniquely represented (all earlier weights are positive), so their multiplicity is at most two. Adding the inner apex and site 1 can each at most double this. Sharpness: the word 10110 has weights \((5,5,2,1,1,1)\) as a multiset and eight representations of 7. In particular for fixed coprime strong flank heights and fixed nonzero determinant half-gap, there are at most 64 contact choices for those two particular flank tournaments. More absolute-height freedom requires varying the cores or allowing a genuinely wider interface.

Run `node two-axis-replay.js` for the rectangle frontier, the explicit absolute inverse on the known collision and small positive cases, and finite state checks for moduli through 128.
