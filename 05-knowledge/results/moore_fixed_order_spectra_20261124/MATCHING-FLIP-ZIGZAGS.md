# The first nonzero matching-flip interaction is an Euler zigzag

This complements the [three-vertex deletion residual and 295-site Moore code](FOURIER-OVERLAPS-AND-THREE-DELETIONS.md). Marked *edge flips* are different from marked *vertex deletions*.

## Exact threshold and an all-orders example

Let E be k pairwise disjoint marked edges of a tournament, k≥1. For S⊆E, let T^S reverse exactly S. Define the mixed path difference
$$
\Delta_E H=\sum_{S\subset E}(-1)^{k-|S|}H(T^S),
$$
where the sum includes S=E. Then:

- At order 2k, the mixed difference is zero.
- At order 2k+1, it is twice the mixed number of directed Hamilton cycles, each cycle counted up to rotation.
- For every k, the bound 2k+1 is sharp. On 0,…,2k start with increasing orientation except for 2k→0. Mark the k edges {i,i+k}, 0≤i<k. Then
$$
                    \Delta_E H=-2E_{k-1},
$$
where E_m is the m-th Euler zigzag number, indexed E_0=1, E_1=1, E_2=1, E_3=2, E_4=5, E_5=16. The residuals start −2,−2,−2,−4,−10,−32,−122,… for k=1,2,…

For the threshold, apply the [Grinberg–Stanley odd-cycle collection formula](https://arxiv.org/abs/2307.05569). Index by potential directed cycles before evaluating edge orientations; summing with alternating signs over S cancels every collection omitting any marked edge. If q odd cycles carry the matching, they need at least 2k+q vertices: an odd cycle of length ℓ supports at most (ℓ−1)/2 disjoint edges. At order 2k+1 only one spanning cycle survives, with weight two. In particular the seven-site example uses marks (0,3),(1,4),(2,5). The path cube is [33,33,69,65,69,71,117,113] and the cycle cube is [1,1,5,4,5,5,13,11], with mask bit i the flip of {i,i+3}.

To evaluate the whole family, contract each marked edge in a surviving directed spanning cycle to an oriented block. Starting at 2k, the first arc is forced to 0, so block 0 points low-to-high and has mixed sign −1. All later blocks have state F (low-to-high, sign −1) or R (high-to-low, sign +1). The next block after 0 must be R. Either state of the last block can return to 2k. For the list of the remaining m=k−1 block indices, ascent and descent transitions, with the sign of the new block included, are respectively
$$
A=\begin{pmatrix}0&1\\-1&1\end{pmatrix},\qquad
D=\begin{pmatrix}0&0\\0&1\end{pmatrix}.
$$
Let e=(0,1) and let u be the column of ones. After removing the leading minus sign, the signed cycle sum is the sum of e·(transition product)·u over permutations of m indices (and is 1 for m=0). Write A=D+C, so C²=−I and D=eᵀe. Expanding by the positions where C is chosen forces runs of increasing indices. For run lengths ℓ_1,…,ℓ_j, there are m!/∏ℓ_i! permutations with those chosen ascents; other comparisons are unrestricted. Each nonfinal run contributes e C^(ℓ−1)eᵀ, and the last contributes e C^(ℓ−1)u. Their exponential series are sin z and sin z+cos z−1. Thus the full exponential generating series is
$$
1+\frac{\sin z+\cos z-1}{1-\sin z}
 =\frac{\cos z}{1-\sin z}=\sec z+\tan z.
$$
By [André's formula for alternating permutations](https://arxiv.org/abs/0912.4240), this gives E_m. Restoring the sign and the cycle weight proves the result. This measures a mixed difference across the cube, not the difference of the two extreme tournament heights and not a general single-arc target selector.

## Keep a marked Hamilton spine in the compact Moore code

The earlier 295-site Fourier chart uses at most 43,206 bits when phase origins are normalized. Reserving all consecutive forward arcs in a linear order on **296** sites still leaves 43,365 free arc bits. Fix 157 for eight disjoint contiguous gadgets and use the others for precisely the same Fourier payload. The decoder still knows the marked vertex order (the display itself need not recover that order from an unmarked tournament).

Take seven eight-site gadgets and one six-site gadget. In each, start with increasing orientation and reverse the nonspine bits in lexicographic pair order selected by mask 1172854 (eight sites) or 163 (six sites). The first-to-last path counts are 24 and 3. After reversing just the local spine they have respectively 23 and 3 last-to-first paths. Concatenating at the prescribed endpoints along the global spine gives at least 24^7·3 and 23^7·3 paths. The second bound counts the converse of the mate that flips every nonspine arc. Thus both endpoints are **above 2^33**, and all 43,206 data bits remain available.

For an upper bound add one apex with both directions to each site. Hamilton paths inject into permutations counted by the cycle-cover permanent. The augmented row degrees are between 1 and 296, and the mean is 296·299/(2·297). By Brégman and AM–GM,
$$
H\le \left[\frac{296\cdot299/(2\cdot297)+5}{e}\right]^{297}<2^{1730}.
$$
The replay certifies (d!)^(1/d)≤(d+5)/e for every integer d in this finite range by exact powers, using rational upper and lower bounds on e from its factorial series. These uniform count bounds do not evaluate the hypothetical graph. Using an unmarked count by itself would not decode this tournament.

Run the accompanying exact replay for the cube, all-orders recurrence controls, code capacity and path bounds.
