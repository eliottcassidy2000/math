# Exact carrying across growing signed seams

Archive anchors, sampled from older seam work: [September 21 arithmetic seams](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/arithmetic_seams_20260921_dynamics.md) distinguishes a genuine scalar trajectory from a projected lift that forgot a sheet; [September 25 Mills shells](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/seam_mills_20260925.md) retains the remainder required to decode an output from a shell. No number-theoretic theorem from either file is needed here. The analogous discarded coordinate is the carry between the large signed term and the two-contact remainder.

## A seam test for growing banks and finite walks

At each vertex of a proposed path of repairs, write the physical count as \(H=M B+C\), with integer \(B,C\), one common \(M>0\), and \(C=h(D)-J(D,V)\). In the sentinel-head architecture, \(M=2^k\) is available whenever both guards have length at least \(k\). Allow the cores, masks and guard lengths to vary along the walk. If exact edge differences all have magnitude at most \(d\), and the residuals range in an interval of width \(w<M-d\), then \(B\) is constant along the walk: on an edge,
\[
 M\,\Delta B=\Delta H-\Delta C,\qquad
 | \Delta H-\Delta C|\le d+w<M.
\]
Consequently a +2 walk has at most \(\lfloor w/2\rfloor\) edges. A core that grows with the height need not have small \(w\); the bound identifies the exact threshold it must cross. This argument does not need a fixed finite bank or asymptotic comparisons.

A modular +2 cycle can nevertheless look complete. On a transitive four-vertex core with both guards of length at least 3, take masks whose ordered two-contact counts \(J\) are \(6,4,2,0\), respectively \(\{0,1,2\},\{0,3\},\{1,2\},\varnothing\). The states give the cyclic residues \(3,5,7,1\pmod 8\). The closing edge at guard heights at least 4 cannot be an exact +2: its residual difference is \(-6\), so it would require \(\Delta B\ne0\) with signed high difference \(8\), ruled out by divisibility by 16. This is a literal example of a full modular loop without the corresponding absolute return. The earlier one-core four-count family in the moving-coefficient packet gives an even closer near-miss: reversing its order gives three literal +2 edges and then a closing edge that reads +2 mod 8 but is actually -6. The cyclic quotient alone cannot distinguish the three edges from the spurious fourth.

## A growing +2 bridge that pays exactly for the carry

Put an unmarked guard of length \(m+1\) on each side of a transitive core \(TT_{2m}\), for \(m\ge1\). Use the sentinel \(t\) and head \(x\) from the old tape: all tape vertices beat \(t\), \(t\) beats \(x\), and \(x\) beats precisely the unmarked tape vertices. The result is strong, since the unmarked first guard is reached from \(x\) and the whole tape reaches \(x\) through \(t\). Compare:

| mask on \(TT_{2m}\) | \(E\) | \(L\) | \(J\) |
| --- | --- | --- | --- |
| the singleton at local index \(m-1\) | \(2^{m-1}\) | \(2^m\) | 0 |
| all \(2m\) positions | \(2^{2m}-1\) | \(2^{2m}-1\) | \(2^{2m}-2\) |

The full-mask row counts all optional ordered two-slot covers with first, second, or both slots nonempty. On a transitive core such covers are simply assignments of the vertices to one of the two slots. The singleton row works by splitting before or after the marked position, and its simultaneous end/start count is zero. The two full readings are, in that order,
\[
 2^{4m+2}-2^{2m}+1,\qquad 2^{4m+2}-2^{2m}+3.
\]
To get the +2, the high signed term **increases by \(2^{2m}\)** and the residual **decreases by \(2^{2m}-2\)**. Thus the original address shell really moves, rather than being frozen. Modulo \(2^{m+1}\) this crosses exactly \(2^{m-1}\) units of the high quotient. The case \(m=2\) gives 1009 and 1011. The masks differ on \(2m-1\) head arcs; no uniform edit-cost claim is made.

## A second bridge pays an odd-times-power-of-two carry

Work directly on the whole transitive tape of order \(N\), indexed from zero. Whenever its first and last vertices are unmarked, the same formula reads
\[
 H=1+2^N + \sum_{v\in V}(2^v-2^{N-1-v})-J(V),
 \quad
 J(V)=\sum_{u<v\in V}\left(2^{v-u-1}+\mathbf 1_{\{v=u+1\}}\right).
\]
For \(J\), one order of the prescribed endpoints permits every intermediate vertex to choose its slot; the opposite order requires adjacency. All vertices exterior to the two prescribed endpoints have forced slots.

Take an even \(r\ge2\) and \(N=3r+3\). One mask is the single position \(r+1\), giving signed contribution \( -2^{r+1}(2^r-1)\), with \(J=0\). For the other mask take two centered consecutive runs, each of length \(r\), separated by one unmarked position. Reflection symmetry makes the signed contribution zero. Within each marked run the interaction is \(2^r-2\); across the runs it is the dyadic geometric product \(2(2^r-1)^2\). Thus the whole \(J\) is
\[
 2(2^r-2)+2(2^r-1)^2=2^{r+1}(2^r-1)-2.
\]
Hence the two outputs differ by +2, now carrying the nontrivial odd multiplier \(2^r-1\) in the high coordinate. At \(r=2\) the outputs are 489,491. The minimal unmarked guard length among these configurations grows at least as \(r/2\).

## A six-edge exact carry walk inside one transitive core (finite witness)

For an old transitive tape of order 20, put three unmarked positions on each side of the central length-14 core. Interpret the following length-14 binary masks least-significant-bit-first in tape index (the rightmost printed bit is position 0 of the core):

1. 11101111010011
2. 11100010000011
3. 11000100010101
4. 01100111101010
5. 01011001010010
6. 11001000110101
7. 10010011111110

The seven physical counts, in this order, are \(1054527,1054529,\ldots,1054539\). The respective J values are \(12194,10912,8270,2268,1658,8648,4702\). Every edge pays a nonzero true carry: since its physical difference is +2 and its J difference is never \(-2\), its signed high difference is never zero. Consecutive masks reverse 4,6,10,8,7,8 literal incident arcs. This is a single finite witness, not a uniform seven-count theorem.

The two infinite constructions prove exact carry cancellation at increasing guard height but only one +2 edge at each height. The finite walk proves several nonzero carries can coexist even without a nontransitive tape block. It remains open here whether an arbitrarily long such walk can be synthesized at a chosen absolute target. The replay checks the mask formula against a separate direct enumeration of all old-slot assignments.

## Exactly three arc flips suffice, and two cannot pay a growing carry

There is a bounded-edit version of the carry. Take \(i\ge2\) and old transitive tape order \(N=3i+3\). The two masks, in increasing count order, are
\[
 U=\{i-1,i+1,2i+2\},\qquad V=\{i-1,2i+1\}.
\]
Each mask leaves at least \(i-1\) old vertices unmarked on both ends. They differ at exactly three vertices, so the lifted tournaments differ by exactly three head-arc reversals. If \(S(V)=\sum_{v\in V}(2^v-2^{N-1-v})\), then
\[
\begin{aligned}
 S(V)-S(U)&=-3\cdot2^i,\\
 J(V)&=2^{i+1},\qquad J(U)=5\cdot2^i+2,\\
 H(V)-H(U)&=(-3\cdot2^i)-(-3\cdot2^i-2)=2.
\end{aligned}
\]
The \(J(U)\) terms come from distances \(2,i+1,i+3\). The first pair is \(401,403\) at \(i=2\). In ascending order the pair is
\[
 2^{3i+3}-6\cdot2^{2i}-7\cdot2^{i-1}-1,\quad
 2^{3i+3}-6\cdot2^{2i}-7\cdot2^{i-1}+1.
\]
Thus a constant-size genuine carry works at arbitrarily large guard height. The actual residual, rather than a residue of it, supplies the tiny imbalance.

This has a sharp edit threshold **for transitive tapes in this head architecture**, if genuine carry means a nonzero change in the high term \(S\). Fix guard length at least \(k\ge3\) at both ends, the same old tape, and consider two masks differing at at most two positions. For any position \(v\) in the allowed region put \(d_v=2^v-2^{N-1-v}\), \(\lambda_v=2^v+2^{N-1-v}\). A pair-energy change is bounded in absolute value by the sum of all possible interactions touching a changed position (overcounting pairs touched twice is harmless). For one side of \(v\), the sum of \(w(j)=2^{j-1}+[j=1]\) through distance \(b\ge1\) is \(2^b\), and for \(b=0\) it is zero. Consequently
\[
 |\Delta J|\le 2^{-k}\sum_{v\text{ changed}}\lambda_v.
\]
If at most two bits change and \(\Delta S\ne0\), then \(\sum\lambda_v\le3|\Delta S|\). For one noncentral distance level, \(\lambda_v\le3|d_v|\); for two positions at the same level their contributions either cancel or add. For two distinct levels, reflect positions to the left half and suppose their separation there is \(r\ge1\). Name the farther-from-center level \(u\) and the nearer \(v=u+r\). Direct subtraction gives
\[
 \frac{\lambda_u+\lambda_v}{|d_u|-|d_v|}=\frac{2^r+1}{2^r-1}\le3.
\]
This also handles the exact center, if present. The signed change is at least the displayed difference in absolute value, since the only coefficients are signs. Therefore if \(\Delta S\ne0\),
\[
 |\Delta H|\ge(1-3\cdot2^{-k})|\Delta S|\ge2^k-3>2,
\]
using divisibility of \(\Delta S\) by \(2^k\). The three-flip construction attains the edit threshold at every \(i\ge4\). This does not exclude two-flip repairs with \(\Delta S=0\), nor apply to different tape orders or changing interior arcs.
