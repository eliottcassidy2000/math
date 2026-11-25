# Seven Gewirtz dualities need Fano sources; four +2 readings remove bad primes

A finite restriction on one surviving Moore transport branch and a complete multichannel decoder for the existing tournament constructions. Neither gives a full Moore deck or an absolute selector for every odd integer.

## 1. The seven-duality case cannot use only pencils, or only one Fano source

Continue the [full Gewirtz transport test](GEWIRTZ-AND-DUAL-GUARDS.md), and consider a sharply transitive set with exactly seven type-swapping elements. At a source point coclique their seven line images have multiplicity vector \(b\) with \(Nb\) odd at every projective point, from \(6a-3Nb=(16-7)\mathbf1\). Let \(E\) be the odd-multiplicity line support. Each point is on at least one member, so \(|E|\in\{5,7\}\). For \(t=N\mathbf1_E\), the projective Gram equation \(N^\mathsf TN=4I+J\) gives \(\sum t_i=5k,\ \sum t_i^2=4k+k^2\), \(k=|E|\). Thus \(u=(t-\mathbf1)/2\) is nonnegative integral with \((\sum u_i,\sum u_i^2)=(2,4)\) or \((7,7)\), respectively.

- **Pencil source:** for \(k=5\) one entry is 2 and the rest zero, giving all five lines through a point, plus two copies of one line.
- **Fano source:** for \(k=7\) there are seven triple points and fourteen single points. Every pair of selected lines meets at a triple point, and each line has three such points. The seven distinct lines form a dual Fano subplane.

The binary line-incidence kernel has weight enumerator
\[
 1+168z^6+210z^8+1008z^{10}+280z^{12}+360z^{14}+21z^{16}.
\]
The replay checks this over all \(2^{21}\) subsets. Its complements are the odd covers, including 21 pencils and 360 dual Fano subplanes. This gives a short local test even beyond seven dualities.

There must in fact be **at least two Fano sources among the 21 points**, and separately among the 21 lines. To prove it, apply the exact small group search below. As the earlier 40320-element census showed, any two correlations in the seven-element set agree on at most one projective source point: their quotient is a point-preserving derangement on hyperovals, with at most one projective fixed point. Thus at most \(\binom72=21\) pair coincidences occur in total. A pencil source with its doubled line outside the pencil uses one coincidence; with the doubled line inside it uses three. A Fano source uses none. If there are zero or one Fano sources, every doubled line must be outside. In the all-pencil case every pair coincides at exactly one point with a different color (the point). With one Fano source exactly one pair does not coincide, and the remaining 20 colors are distinct.

### Exact colored-graph check

Translate the seven correlations on the left by the inverse of one of them, so they become seven point-preserving permutations including the identity. They must still be mutually deranging on the 56 hyperovals. Translation preserves colors and converts the line pencils into five collinear points while taking the dual Fano subplane to an ordinary subplane.

1. With no uncolored pair, the six other elements are among **7560** derangements with one projective fixed point. Adjacency of candidates requires hyperoval derangement relative to each other and one projective coincidence, giving that coincidence's color. Delete conflicts with either candidate's own fixed-point color. The resulting graph has **4,999,680** edges. The candidates fall into two point-preserving conjugacy classes of sizes 2520 and 5040. Choose one representative as the first nonidentity point of the search. Depth-first intersect candidate sets and accumulate all edge colors; reject any repetition. No six candidates survive.
2. With one uncolored pair, choose one endpoint as the identity and the other among the four conjugacy classes of zero-projective-fixed-point derangements, each of size 2880. Every other element is a fixed-one candidate and all its further pair edges are colored. Use the same search, and add the geometric pruning rule at each source. Once a duplicate occurs, the remaining images must extend to a line excluding the duplicated point. With no duplicate yet, the current distinct images must lie in a Fano subplane or be collinear after deleting at most one point. These tests are necessary for the pencil or Fano completions; a branch that fails cannot be repaired by adding more images. The complete search rejects all four representatives.

The receipt for case 1 visits 2, 2744, 253139, 972536, 15024 search nodes at selected-candidate depths 1–5, then none at 6. For case 2 it visits 4, 5528, 230185 at additional candidate depths 0–2, then none at 3. Ignoring the geometry in case 2 **does** give a colored-pattern witness; the incidence condition is necessary for that exclusion. Conjugating by a duality gives the same result for source lines. Translating a 49-duality sharp set by a duality gives the complementary seven-duality condition. A full seven-duality example with at least two Fano sources is not excluded; the full Gewirtz group is not ruled out by this search.

Replay (inside this directory):

    python3 gewirtz-fano-export.py rainbow-data.bin
    c++ -O3 -std=c++17 gewirtz-fano.cpp -o gewirtz-fano
    ./gewirtz-fano rainbow-data.bin

## 2. Four actual +2 channels, and why each is needed for this universal decoder

Use even payload \(B_0=2v\) with digit counter \(q(v)\) from the [short-guard construction](TWO-COORDINATE-RIGIDITY.md). The first three lower readings are from the established charts A at guard width 1, B at width 1, and A at width 2. Normalize by subtracting their odd constants and dividing by 2:
\[
 X=(H_{A,1}-427)/2=831v-29q,\quad
 Y=(H_{B,1}-493)/2=798v-29q,\quad
 Z=(H_{A,2}-743)/2=1471v-29q .
\]
The last lower reading is from a one-apex contact seam with flank integers \(2v,2v+1\), whose ceil-half difference is one. By the [two-flank formula](FIBRE-TRANSFER.md) and the binary height recursion,
\[
 W=(H_{\rm flank}-1)/2=4v^2+\lfloor9v/2\rfloor-2q .
\]
Each channel is a separately constructed pair of actual tournaments one physical reversal and two Hamilton paths apart.

The inverse is exact:
\[
 v=28(Z-X)-543(X-Y),\qquad
 q=831v-X-14(4v^2+\lfloor9v/2\rfloor-W).
\]
For measurements of the four physical lower heights modulo \(4M\), first normalize modulo \(2M\); the formulas recover \(v\bmod2M\) and \(q\bmod M\). The doubled modulus for \(v\) is sufficient for the floor. Every such state is reachable after any leading prefix by the previous routing theorem. For exact heights, also require positive integral \(v\), the prefix when specified, and \(q=q(v)\); recompute all four readings before accepting them. These are joint inverse measurements, not arbitrary independent heights.

Moreover **no triple of these four specified channels** has the same universal modular recovery property. Without \(W\), all three remaining channels are blind to a change of \(q\) by 2 modulo \(29\), while agreeing modulo \(58\) after normalization. Otherwise the two retained guard slopes coincide at an odd prime \(p\ne29\) dividing their difference: \(p=3,5,673\) for the three choices. With \(v\) held on one parity and \(q\) chosen to maintain the common guard reading modulo \(p\), \(W\) reduces to a quadratic in \(v\) with nonzero leading coefficient, which cannot be injective over \(\mathbf F_p\). Choose the representatives of \(v\) divisible by 4 and \(q\) even to match the extra dyadic bit of the normalized readings. The replay exhibits actual modular collision states for each omitted channel, in omission order \(X,Y,Z,W\):
\[
\begin{array}{c|c|c}
 M &(v,q)&(v',q')\\ \hline
673&(116,500)&(120,378)\\
5&(0,0)&(4,6)\\
3&(4,0)&(8,0)\\
29&(4,0)&(4,2)
\end{array}
\]
All are realizable as congruence states by routing at modulus \(2M\), which ensures a valid lift for the floor. This minimality is only among subsets of **these channels**, not among all tournament constructions.

For \(v=5\) the actual four pairs are **8679→8681**, **8415→8417**, **15395→15397**, and **241→243**; independently recounted by endpoint-subset DP. The first three use word 1010. The last uses words 101000010 and 101010010 with only the central contact changed. Run `node four-seam-decoder.js` to check all four witnesses, 4096 digit values, modulus lifting, and the four indispensable-channel fixtures.
