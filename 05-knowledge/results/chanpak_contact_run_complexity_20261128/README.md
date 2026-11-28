# Contact-run complexity and a dyadic inverse for Chan–Pak

**Results.** Fixed feedback-vertex budget and fixed contact-run budget allow only finitely many heights supported on a fixed prime alphabet. For one apex, the obstruction is effective: if \(H\) is supported on a fixed prime set and sufficiently large, every witness requires
\[
r\gg_P\sqrt{\frac{\log\log H}{\log\log\log H}}
\]
active runs. In the other direction, every even run budget \(r\) gives a directly invertible family with the sharp growth order \(\Theta_r((\log X)^r)\). An exact census finds the one-apex witnesses at \(7\cdot3^t\) through \(t=13\).

## 1. Bounded contact complexity, not just bounded reversal cost

Let \(F\) be a set of \(q\) vertices whose deletion leaves a transitive tournament. In its linear order, partition the remaining vertices into maximal runs of equal \(q\)-bit incidences with \(F\). Each run is a transitive module. Fix the orientations on the named vertices of \(F\) and the ordered run types; these are a finite set of templates at bounded \(q,r\).

In a Hamiltonian path, compress each maximal uninterrupted piece of a run to a symbol. A skeleton uses each exceptional vertex once and visits each run \(i\) some \(k_i\ge1\) times. Between consecutive visits to the same run there must be an exceptional vertex: successive distinct ordinary runs alone strictly increase. Thus \(k_i\le q+1\). With \(a_{\boldsymbol k}\) skeletons of a given visit vector, the exact count for run lengths \(\ell_i\ge1\) is
\[
H(\boldsymbol\ell)=\sum_{\boldsymbol k}a_{\boldsymbol k}
\prod_i k_i!\,S(\ell_i,k_i).
\]
The factor counts surjections from the ordered vertices of the run to its distinguishable visits, since within a visit their order is forced.

When some \(k_i\ge3\), use the first such \(i\) and interchange the first two excursions between visits to \(i\). Their endpoints are the same symbol, so the resulting skeleton is valid and the visit vector is unchanged. They contain disjoint nonempty sets of named exceptional vertices, so this is a free involution. Hence \(a_{\boldsymbol k}\) is even. Expanding the Stirling factors by inclusion–exclusion gives
\[
H=\sum_{\boldsymbol b} c_{\boldsymbol b}\prod_i b_i^{\ell_i},
\quad 1\le b_i\le q+1.
\]
For all-odd base vectors the only odd coefficient is at \((1,\ldots,1)\): the \(k_i\ge3\) vectors contribute evenly, and any \(k_i=2\) contributes an even coefficient at its only odd base, by \(2^{\ell_i}-2\). The all-single-visit term is the odd Rédei count of the quotient tournament.

Suppose infinitely many different heights are supported on a finite prime set \(P\). Pass to one template, a prime \(p\in P\) with \(v_p(H)\to\infty\), and a subsequence on which each length is either constant or tends to infinity. In the equation obtained by adjoining the term \(-H\), partition the nonzero terms into minimal vanishing blocks, passing to a subsequence with fixed partition. All terms are signed units in a fixed finitely generated group (adjoin primes of the bases and coefficients to \(P\)). In each block, divide by one of its terms and apply [Evertse–Schlickewei–Schmidt](https://annals.math.princeton.edu/2002/155-3/p04) to the resulting nondegenerate equation. There are finitely many ratios, so each block lies entirely among terms with bounded \(p\)-valuation or entirely among terms whose valuations tend to infinity. The bounded terms therefore sum to zero. But they include the all-ones term, with odd coefficient; every other term in their sum is even, and \(-H\) is not in their sum. Contradiction.

This strictly extends the [fixed-reversal result](https://github.com/eliottcassidy2000/math/commit/1a976f9964d6186512c0c4ddadb9739c7da982dc). The endpoints of a fixed set of backward arcs leave a transitive spine, with a bounded number of contact runs since all surviving incidences are forward in the original order. Conversely the one-apex word \(1^m0^m\) has just two runs but needs \(m\) backward arcs in every optimal order. Pair the \(m\) 1-sites with the \(m\) 0-sites to obtain arc-disjoint triangles through the apex; each needs a backward arc. Placing the apex before the spine attains \(m\).

**Bounded-parameter exact inverse.** A run on a triangle of the quotient can have length at most \(\lfloor\log_2(H-1)\rfloor\) for \(H>1\): induce on the run and one site from each of the two other parts of the triangle. The induced count is \(1+2^{\ell_i}\) by the two-run formula below with lengths \(1,\ell_i\), and Hamiltonian path counts are monotone under taking induced subtournaments (insert omitted vertices one at a time). A quotient vertex on no triangle is a singleton strong component, because its in-neighbors dominate its out-neighbors. Its run can therefore be shortened to length one without changing \(H\). This gives \(O_{q,r}((1+\log H)^r)\) tuples for a complete inverse, for given budgets. The case \(H=1\) has the one-vertex witness.

## 2. One apex: signed intervals and a direct inverse

Let 1 mean the apex beats a site. Initial zeros and final ones are ordinal padding, so write the nonempty active word as \(1^{\ell_1}0^{\ell_2}\cdots 1^{\ell_{r-1}}0^{\ell_r}\) with \(r\) even. Then
\[
H=r+1+\sum_{1\le i\le j\le r}(-1)^{j-i+1}2^{\ell_i+\cdots+\ell_j}. \tag{1}
\]
Indeed, the parts of a Hamiltonian path before and after the apex must be increasing. If both are nonempty, let \(i\) be the last site before the apex and \(j\) the first after. For \(i<j\) they must be adjacent, giving one path for a boundary 01. For \(i>j\), an inversion 10 in spine order, the sites between \(j\) and \(i\) can each be assigned independently to either side. There are also two paths with an empty side. For a 1-run \(a\) and later 0-run \(b\), the weighted inversions sum to
\[
(2^{\ell_a}-1)(2^{\ell_b}-1)2^{\ell_{a+1}+\cdots+\ell_{b-1}}.
\]
Expansion contributes each nonempty interval exactly once: choose the left 1-run at its left end if odd or one earlier if even, and the right 0-run at its right end if even or one later if odd. Each adjacent 10 gives one empty-interval term. With the \(r/2-1\) boundaries 01 and two empty-side paths this gives (1).

If \(0<\ell_1<\cdots<\ell_r<\ell_1+\ell_2\), the increasing singleton exponents precede every interval of length at least two. Starting at \(z=H-(r+1)\), read \(v_2(z)\), add back its power of two, and repeat \(r\) times. The decoded lengths give an actual witness exactly when they satisfy the inequalities and re-evaluation of (1) returns \(H\). This is a deterministic integer decoder for this family, not just a residue selector. Take one length from each consecutive disjoint bin of width \(\lfloor M/(2r)\rfloor\) inside \([M,3M/2)\); all decode separately and give heights \(<2^{3rM}\). The first and last sites give an inversion of weight \(2^{n-2}\) for width \(n\), so \(n\le2+\lfloor\log_2(H-1)\rfloor\). Counting length tuples proves the matching upper bound. Taking \(M=8r\) also gives \(\exp(\Omega(\sqrt{\log X}))\) distinct heights up to \(X\) with no bound on \(r\).

## 3. Quantitative prime-support obstruction for one apex

The interval formula gives at most \(K=r(r+1)/2+O(\log(r+1))\) signed powers of two for \(H\), expanding the constant \(r+1\) in binary. Cancel opposite signs at equal exponents and carry equal signs to the next exponent. This terminates since the term count drops, yielding \(H=\sum_{j=1}^k\varepsilon_j2^{e_j}\), \(k\le K\), \(e_1>\cdots>e_k=0\), \(\varepsilon_j\in\{\pm1\}\). The last exponent is zero by oddness. The leading sign is positive since lower distinct powers have sum \(<2^{e_1}\). The original exponents are bounded by \(n\) (including the binary expansion of \(r+1\)); the sum of absolute powers never increases on normalization. Thus \(e_1\le n+\log_2 K=\log_2 H+O(\log\log H)\). Also \(e_1\ge\log_2 H-1\).

Put \(g_i=e_1-e_i\) and \(P_i=\sum_{j\le i}\varepsilon_j2^{e_j}=2^{e_i}u_i\); then \(1\le u_i<2^{g_i+1}\). For \(i<k\) the omitted tail is nonzero, again by distinct powers, and \(|P_i/H-1|<2^{e_{i+1}+1}/H\). Write \(H=\prod_{p\in P}p^{v_p(H)}\). [Matveev's logarithm bound](https://doi.org/10.1070/IM2000v064n06ABEH000314), in its positive-rational product form, applied to
\[
P_i/H=u_i\,2^{e_i}\prod_{p\in P}p^{-v_p(H)}
\]
gives \(\log |P_i/H-1|>-C_P(1+g_i)(1+\log\log H)\) for large \(H\). The fixed bases have constant logarithmic heights; \(u_i\) has height \(O(1+g_i)\), and all exponents are \(O(\log H)\). Omit \(u_i\) when it is one, or any zero exponent. Comparing the two bounds and using the estimate for \(e_1\) gives
\[
1+g_{i+1}\le C'_P(1+\log\log H)(1+g_i).
\]
Iterate from \(g_1=0\) to \(g_k=e_1\gg\log H\). Thus \(k\gg_P\log\log H/\log\log\log H\), and \(K=O(r^2)\) proves the claimed lower bound. In particular along \(H=7\cdot3^t\) the active run requirement is \(\gg\sqrt{\log t/\log\log t}\) for large \(t\). The constants are effective; no numerically useful threshold is claimed. This is the signed-power version of the gap argument in [Bugeaud–Kaneko](https://irma.math.unistra.fr/~bugeaud/travaux/NonzeroSunitsRev.pdf); ordinary binary Hamming weight cannot substitute for signed weight, since subtraction creates arbitrarily long borrow strings.

## 4. Exact one-apex census on the third ternary chain

For a prefix starting at 1, keep its binary value \(x\), last bit \(b\), and path count \(H\). Appending \(e\) increases \(H\) by \((1-e)(x+b)\); then update \(x\gets2x+e,\ b\gets e\). Indeed each new terminal zero adds the weighted inversions whose total is \(x\); the net change from the empty-side and 01-boundary counts is \(b\) when \(e=0\), zero otherwise. This is the prefix recurrence from the [third-chain construction](../session302_ternary_ray_birth_20261113/README.md).

The replay examines all 7,735,713 active words with count \(\le7\cdot3^{13}\), using the width bound 25 and pruning only when the monotone prefix count exceeds the target. These are all the hits; each displayed word has exactly one other active word at that height, its reverse-complement.

| \(t\) | \(H\) | Active word | Runs |
|---:|---:|:---|---:|
| 2 | 63 | 1001100 | 4 |
| 8 | 45,927 | 10010010011000110 | 10 |
| 9 | 137,781 | 101010000111001100 | 10 |
| 12 | 3,720,087 | 10001111110011110101100 | 10 |
| 13 | 11,160,261 | 110000101011100000000100 | 10 |

Run with `node replay.js`. The replay also compares (1) and the prefix recurrence on all 32,767 active words of widths 2–16, with direct subset-DP comparisons through width 8, and performs 1,500 direct decodes. The three higher one-apex hits continue the ten-run pattern from exponent 8. The finiteness theorem rules out any fixed run count for an infinite subsequence of this chain, even when one allows an arbitrary bounded exceptional set. None of the results settles whether one apex meets the chain infinitely often, or supplies a universal inverse for Chan–Pak.
