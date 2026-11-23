# The second base-16 Kaprekar orbit: six reversals on eight vertices

The second width-two base-16 gap orbit is \(3,11,5,7,3\) under \(g\mapsto|2g-17|\). Its required ordinary Hamiltonian-path counts are \(45,165,75,105,45\). As in the [first orbit](BASE16-CYCLE.md), prescribed heights are checkpoints in a physical arc walk rather than operations on counts alone.

## A sharp labelled witness

Start from the eight outgoing rows \((250,252,249,112,160,192,144,8)\), with vertices numbered from zero. Reverse these arcs in order:
\[
(0,3),(1,7),(1,2),(0,3),(1,2),(1,7).
\]
The count trace is
\[
45\longrightarrow87\longrightarrow165\longrightarrow75
  \longrightarrow33\longrightarrow105\longrightarrow45.
\]
Every flipped arc is reversed twice, and the actual tournament returns unchanged. The detour through \(33=3(12-1)\) is another literal Kaprekar phase: its gap 3 is on the base-12 orbit \(1,11,9,5,3,7,1\).

Two independent obstructions make six optimal, even when an endpoint isomorphism is allowed:

1. No single arc leads from height 75 to 105.
2. No height-45 tournament has both a single-arc neighbor of height 165 and one of height 105.

The first forces an extra edit between 75 and 105; the second forces a different extra edit at one of the two transitions incident to 45. In a return up to relabelling, transport the last incident arc by the endpoint isomorphism to compare it with the first arc at the same host. Four prescribed transitions plus two distinct extras give six.

At length six, those two forced extras leave at most two edits for \(45\to165\). At order seven the complete height-45 source deck has nine shapes and gives no height-165 arrival among its 189 single edits or 1890 unordered distinct edit pairs. Smaller witnesses would pad to order seven. Therefore **eight vertices is minimal at the optimum**. The jump \(45\to165\) alone does admit one arc, but not until order ten: the complete order-nine deck has 53 shapes and no one-arc arrival. At order ten, rows \((1006,1020,1016,1008,993,960,384,768,512,64)\) and arc \((1,5)\) witness it. This distinguishes the two separate quotient edges from a compatible turn through a common state.

## Reducing the obstructions to finite exact checks

The strong-component product law and [Busch's exact strong minima](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) give minimal strong heights \(9,15,25,45,75,125\) at orders \(5,\ldots,10\). For height 45 the possible nontrivial SCC height lists are \((45),(3,15),(5,9),(3,3,5)\), and for 75 they are \((75),(3,25),(5,15),(3,5,5)\), each up to ordering. The sum of the nontrivial orders is at most 10 and 11, respectively. Enumerate all ordered lists of actual cores and all positions for ordinal singleton pads.

**Incompatible turn at 45.** The two positive gains 120 and 60 have 2-adic valuations 3 and 2, and neither ratio is exceptional in the [singleton budget](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget). For each prospective edge, keep at most four and three singleton vertices, respectively. Retain the union: the gain test for each edge keeps all its endpoints and all singleton intermediates; keeping extra singleton vertices outside its affected interval does not change either height. An internal edge needs no singleton. Pad the reduced source to order 17. There are 27 ordered core templates and 1069 padded shapes, with 160 arcs to 165 and 859 arcs to 105, but **zero sources admitting both**.

**No \(75\to105\) single edge.** Since the gain 30 has valuation one and the ratio is not three, keep at most two singleton vertices. Except for a strong order-nine height-75 core, the 26 ordered core templates give 271 order-13 padded shapes with no arc to 105. For the exceptional strong order-nine core, a cross edge that misses or does not span it leaves a factor 75. If the reversal spans the core or touches it at an endpoint, it joins the whole affected SCC interval into a strong tournament with at least ten vertices: there is an alternate forward route inside the interval plus the new backward arc. Its height is at least 125. Thus an edge out of this source can reach 105 only by an internal reversal. The complete anchored order-nine deck gives 225 strong height-75 orientations, keeping duplicates under relabelling; none of their 8100 internal arcs leads to 105.

**Replay.** Run [`node kaprekar-base16-second.js`](kaprekar-base16-second.js). The core generator fixes a Hamilton path \(0\to1\to\cdots\to n-1\) and enumerates the free arcs. Insertion of a vertex extends each induced Hamilton-path order without conflating distinct orders, so it is safe to discard a frame if a path-framed one-vertex deletion is already above the height threshold. For threshold 75 the order-seven and order-eight banks have 8033 and 18497 anchored prefixes; 13100 strong anchored order-nine candidates pass the deletion filter before the final height check. The small cores are deduplicated by relabelling for the padded banks. The padded scans use the exact rooted gain law with subset tables for each core, not an order-17 tournament table; the order-nine internal checks and literal walks use ordinary subset/end-vertex DP.

## A general way to solve joint padding constraints

For a fixed ordered nontrivial core list, prescribe a cross-SCC arc by its concrete rooted core endpoints or by a singleton endpoint in a named gap. For a separated arc its gain has the form \(A\,2^m\), where \(A\) comes from the cores, rooted endpoint terms, and untouched factors, and \(m\) is the number of intervening singletons. Thus an exact target gain first fixes \(m\) or rules the skeleton out. Let \(x_j\) be the available singleton count in gap \(j\). Fully interior gaps contribute their sum \(L\); a gap used for a singleton endpoint contributes any value from 0 to \(x_j-1\). For different endpoint gaps the necessary and sufficient padding inequalities are
\[
x_j\ge1\ \text{at endpoint gaps},\qquad
L\le m\le L+\sum_{\text{endpoint gaps}}(x_j-1).
\]
An adjacency branch has \(m=0\), no intermediate core, and uses the adjacent gain formula with its loss term. If both singleton endpoints lie in the same gap, the separated branch instead requires \(x_j\ge m+2\). All other gaps have \(x_j\ge0\).

The sums here are consecutive gap increments. With prefix variables \(P_j=\sum_{i<j}x_i\), the inequalities become difference constraints on \(P_j\), including an optional total size cap. For simultaneous independent outgoing arcs, enumerate the rooted skeletons and use one joint feasibility test: negative cycles obstruct feasibility, and feasible integer potentials specify padding. The arcs can choose positions independently in the same original pads; this does not propose executing them consecutively.

For bounded-edit *walks* rather than simultaneous outgoing arcs, another finite cap preserves return up to relabelling. Keep every edited vertex and the nonsingleton support at every named checkpoint, including both endpoint occurrences. At most \(2r+\sum_i\lfloor3\log_3 H_i\rfloor\) vertices remain for \(r\) edits. An omitted vertex is untouched and is singleton at all named checkpoints, so deletion preserves those heights. Moreover a return isomorphism fixes each omitted vertex: it has the same outdegree at both endpoints, and the outdegrees of singleton SCC vertices are pairwise distinct. Thus restriction preserves the return isomorphism. This gives an all-order finite test without imposing any count bound on the unnamed detours.

## Both base-16 orbits lift together

The [exact transport rule](BASE16-CYCLE.md#transporting-a-literal-orbit-to-larger-bases) for an already-scaled base-\(b\) walk asks for an attained ordinal multiplier \(d(B-1)/(b-1)\), where \(B=d(b+1)-1\). A general sufficient form, for nonnegative even \(x\), is
\[
d=(b-1)x+1,\qquad e=(b+1)x+1,\qquad B=(b^2-1)x+b.
\]
Now \(B-1=(b-1)e\), so the exact multiplier is \(de\), not merely a residue or a rational scaling. Any two attained factors \(d,e\) suffice.

For \(b=16\), take \(x=2^k+2\), \(k\ge6\). The previously derived [signed-shell rule](BASE16-CYCLE.md#transporting-a-literal-orbit-to-larger-bases) realizes both factors:

| Factor | Core outgoing rows | Tag mask | \((D,G,L)\) | Transitive filler order |
|---|---|---:|---|---:|
| \(d=15\cdot2^k+31\) | \((70,20,8,19,37,79,30)\) | 6 | \((81,480,50)\) | \(k-5\) |
| \(e=17\cdot2^k+35\) | \((6,12,8,17,39,79,31)\) | 2 | \((35,136,0)\) | \(k-3\) |

Ordinally joining these two shells to either labelled walk carries both base-16 gap orbits to each base \(B=255(2^k+2)+16\), starting at 16846. For the second walk this gives \(2k+16\) vertices. Both walks still have six reversals; this is an upper bound in the larger bases, not a proof that six is optimal there.
