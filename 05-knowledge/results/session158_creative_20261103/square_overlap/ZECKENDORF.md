# Square-sum Zeckendorf tournament code

Status: **proved construction**. Inspired by the nonconsecutive Zeckendorf normal form, binary 1-2-4 weights, and square-sum paths. This proves an injective family of path counts, not complete odd-count coverage; novelty is not claimed. Except for the core and the first attachment, these tournaments are generally not strong.

Let F₂=1, F₃=2, F_(i+2)=F_(i+1)+F_i. Write N=∑_(i in I) F_i in its Zeckendorf form, where the finite set I⊂{2,3,...} has no consecutive indices; take I empty for N=0. Put T_n=n(n+1)/2 and S={i+2:i∈I}.

On the positive labels through the largest T_(n+1) with n∈S (or through 10 when S is empty), orient every edge from smaller to larger, except reverse exactly these:

* 8→1 and 10→6;
* T_(n+1)→T_n for every n∈S.

All reversed edges have square-sum endpoints: 1+8=9, 6+10=16, and T_n+T_(n+1)=(n+1)². Their tournament Hamilton-path count is

    H_N = 265 ∏_(i∈I) (1 + 2^(i+2)).                         (1)

The map N↦H_N is injective. Indeed divide H_N by 265. For any nonempty product of distinct factors (1+2^k), the smallest exponent is ν₂(product−1). Divide by that factor and repeat to recover all exponents, I, and hence N. Thus N is carried by a canonical nonconsecutive digit pattern in an actual tournament invariant. This code has O((log(N+2))²) vertices, because the maximum Zeckendorf index is O(log(N+2)) and the largest triangular label is quadratic in that index.

Examples for N=0,1,2,3,4: I is {}, {2}, {3}, {4}, {2,4} and H_N is 265, 4505, 8745, 17225, 292825, respectively. In particular, N=1 uses labels 1,...,15 and the square-sum backward edges 8→1, 10→6, 15→10; its tournament is strong. No general strong-connectivity claim is made for the whole family.

## Proof of the count

The first two reversed edges give the strong crossing core [1,10] with exactly 265 paths, by the separate crossing formula in NOTE.md. A single interval of distance d with only its extreme edge reversed has 1+2^(d−1) paths: either no reversed edge, or choose independently which intermediate vertices come before its reversed edge. Such intervals that are strictly separated concatenate, so their path counts multiply.

The indices in S are nonconsecutive. Consequently the added intervals [T_n,T_(n+1)] are pairwise disjoint. The only potential contact with the core is n=4, where T₄=10 and T₅=15. We verify that this shared-endpoint attachment still multiplies by 1+2⁴=17.

If the final reversed edge 15→10 is absent from a Hamilton path of the 15-vertex tournament, every core vertex occurs before 11,...,15, giving 265 paths. Otherwise use 15→10. On deleting 11,...,15, the remaining order is a core Hamilton path except potentially at the predecessor of 10. The only lower vertex to which 10 points is 6, so an invalid projected adjacency must be 6,10. Furthermore, 10 cannot then have a later core vertex, since its only outgoing edge to the core is 10→6 and 6 is already used. This exceptional case would therefore be a path on [1,9] ending at 6.

No such path exists: after 10 is removed, only the edge 8→1 is reversed. To place vertex 9 before terminal vertex 6, a path must descend using 8→1, but 9 cannot occur before 8 in a path with no other descent. Therefore all projections are valid. Conversely, each of the 265 core paths expands uniquely by placing each of 11,12,13,14 either immediately before 15 (after the core prefix) or after the core suffix, in increasing order. Thus there are 2⁴·265 paths using 15→10, and the total is 17·265.

If n=4 is selected, all later intervals are separated from the enlarged core [1,15]. If n=4 is absent, every selected interval is separated from the original core [1,10]. In either case the selected intervals are mutually separated, so their Hamilton-path counts multiply. For interval [T_n,T_(n+1)] the difference is n+1, so its factor is 1+2^n. Substituting n=i+2 proves (1).

## Checks

`python3 zeckendorf_check.py` independently computes the 10-vertex core count, checks that the core-minus-10 endpoint count at 6 is zero, and computes the full 15-vertex attachment count by subset DP. It verifies Zeckendorf encoding/2-adic decoding for N=0,...,9999. These finite checks do not substitute for the proof.
