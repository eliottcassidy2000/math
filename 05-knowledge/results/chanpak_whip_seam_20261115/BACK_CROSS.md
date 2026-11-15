# One backward seam edge: exact gain, and a compensated closure

The whip analogy here is **a literal edge that reaches back across an ordinal seam**. A long forward tail by itself cannot hide its chronological marks in the right-state start coordinate, by the previous [fold-closure obstruction](../chanpak_tail_return_20261115/FOLD_CLOSURE.md). The key distinction is not whether one edge changed, but whether the edge genuinely breaks the ordinal decomposition.

## Rank-two gain of one backward edge

Start from a nonempty ordinal sum (A\to B), choose (a\in A,b\in B), and reverse just (a\to b). Denote by (E_A,O_A) the even/odd **vertex-order** weighted starts at (a) in (A), each weighted by the Hamilton count of its complement. Denote by (E_B,O_B) the analogous weighted ends at (b) in (B). Put (S=E+O,d=O-E) on both sides. Then the exact increase in unmarked spine count is
\[
h_{\rm new}-h_Ah_B=2(E_AO_B+O_AE_B)=S_AS_B-d_Ad_B .
\]
Use the [odd-cycle packing expansion](https://arxiv.org/abs/2307.05569): the old noncrossing packs all survive; any new mixed directed cycle uses the distinguished backward edge. Removing it leaves one nonempty path starting at (a) in (A), then one nonempty path ending at (b) in (B). Their orders have opposite parity for an odd cycle. Conversely any such pair gives exactly one crossing odd cycle, and its complement splits ordinally because both endpoints of the sole backward edge have been removed. The only path pair whose putative forward seam is the removed edge consists of two singleton paths and has even parity, so it does not enter the formula.

Singleton paths give (O_A,O_B>0). An even start exists iff (a) has an outneighbor inside (A): a two-vertex path suffices, and any longer start must first leave (a). Similarly an even end exists iff (b) has an inneighbor in (B). Therefore this increase is strictly positive except when (a) is the sink of (A) and (b) the source of (B). In that exceptional case the reversal just swaps the two adjacent singleton SCCs. In particular a genuine one-edge recontact from the ordinal spine cannot keep the *original* (h).

For an ordinal right tail (B=TT_L) and the back-edge origin at position (j) measured from zero, the exact price is
\[
\Delta h = \begin{cases}
S_A-d_A,&j=0,\\
2^j S_A,&1\le j<L .
\end{cases}
\]
Indeed the weighted ending paths choose any subset of the (j) preceding tail vertices, in order. For positive (j), even and odd ending weights are both (2^{j-1}). This is exact geometric focusing by the depth of the return, not a claim that an arbitrarily large return fixes the state.

A broader necessary guard holds for arbitrary recontact **only between** original SCCs, leaving the interior of every SCC fixed. The Hamilton count cannot fall, and it stays equal iff the result just reorders the old SCCs. Old packings still survive, whereas a mixed directed triangle strictly increases the count. If there are no new mixed triangles, each exterior vertex has uniform orientation to each old SCC: a mixed cut of a strongly connected component must have an edge from a vertex the exterior vertex beats to one that beats it. This would form a mixed triangle. Thus every inter-SCC pair is uniform and the reduced quotient is triangle-free, hence transitive. This does not rule out compensating by a simultaneous edge change *inside* a component.

This distinction was reinforced by two seeded random old-note draws: [the rhombus carrier](../rhombus_overlap_elliptic_20261011.md) distinguishes adding a genuine cross-edge from obtaining a genuinely new geometric carrier, while [adaptive Collatz collision DP](../collatz_collision_dp_20261007.md) retains exact collision coordinates and available partner lengths instead of testing against one fixed partner. Here the real guard is the integer profile, not the presence of a nominal seam reversal.

## Both arms can pay the same price: a literal return with opposite first marks

Take the six-vertex body with outgoing rows (38,16,2,7,13,30), with (h=23). Its even-start weights at vertices 1 and 2 are both 9; reverse the first tail vertex's contact into either vertex. Relative to the plain ordinal insertion, both arms therefore increase (h) by (18), to (41). The single-back-edge extended spine rows and complete apex masks are:

| First tail vertex beats | Extended spine rows | Apex mask | Right state and count ((h,N,S,d)) |
|---|---|---:|---|
| Vertex 2 | (102,80,2,71,77,94,4) | 13, tail mark **0** | (41,221,295,17) |
| Vertex 1 | (102,16,66,71,77,94,2) | 112, tail mark **1** | (41,219,295,23) |

The original body is identical; the two branches differ in which old vertex the tail end reaches back to and in the apex-mask contacts. Unlike repainting a still-ordinal tail, this really closes arms having opposite first chronological bits. Their only unmatched right datum is the signed bias (d). A common next singleton painted 1 does not change (N), sends (S\mapsto 2S+h=631), and resets (d=h=41) in both. A common 0 then gives **891 and 893**, with exactly the same ((h,S,d)=(41,1262,0)); both graphs are strong. The +2 persists for every common transitive suffix, with exact outputs
\[
891+1262t+41p(t)+\{0,2\},\qquad t\ge 0 ,
\]
where (p) is the [binary suffix recursion](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md). It is a literal compensated seam return at a *new* (h), not a closure of the inherited eight-count window.

Replay: `python replay_back_cross.py` tests the exact weighted endpoint formula on seeded literal ordinal sums, the depth law for transitive tails, and the above populated states with direct subset DP before and after alignment.
