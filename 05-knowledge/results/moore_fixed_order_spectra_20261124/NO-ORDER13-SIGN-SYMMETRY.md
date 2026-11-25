# The Gewirtz branch-sign model has no order-13 symmetry

This strengthens [the earlier exclusions](NO-LARGE-PRIME-SIGN-SYMMETRY.md). It applies to the hypothetical full Gewirtz transport **branch-type matrix**, not to arbitrary Moore-graph automorphisms. The group acting on branch signs up to gauge now has no prime divisor \(p\ge13\). The smaller-prime cases still have open orbit counts.

Throughout, \(C\) has diagonal \(-1\) and entries \(\pm1\), is symmetric, and obeys \(C^2=0\bmod3\). Admissible swap counts lie in \(10,13,\ldots,46\). For distinct rows at Hamming distance \(w\), their route count is \(w+1\) for an off-diagonal negative sign and \(w-2\) for a positive one. In particular all row distances are multiples of three. Gauge away any order-13 twist as in the preceding proof.

## 1. One, two and four moving cycles

The earlier one-frequency scalar argument excludes one moving cycle at every prime other than 3. With exactly two moving cycles of length \(p\), trace zero at nonconstant frequencies forces both internal sign blocks entirely negative. Indeed the sum of their two length-\(p\) polynomials is constant at every lag modulo 3, with lag-zero value \(-2\), and the three possible integer sums are distinct modulo 3. Within either orbit the row distance comes entirely from shifting the cross block; it is even as well as divisible by three. The negative direct sign then requires every nonzero cyclic shift distance to be at least 12. But a length-\(p\) binary word of weight \(j\) has averaged nonzero shift distance \(2j(p-j)/(p-1)\le(p+1)/2\), impossible for \(p=5,7,11,13,17,19\).

With four moving 13-cycles there are five fixed sites. Gauge each fixed row to agree positively on the first cycle. The five profiles must be distinct (equal profiles would give distance at most five). Choose five of the eight normalized profiles, and a symmetric five-by-five fixed sign submatrix with prescribed negative diagonal. For each fixed pair its Hamming distance is \(13\) times the profile distance plus the five-coordinate row distance. Exactly **57,344** possibilities (\(\binom85 2^{10}\)) are checked, with none having all pairwise distances divisible by three. This is a necessary-condition test independent of unknown cross blocks.

## 2. Three moving cycles: reduce the fixed block to two profile patterns

With three moving cycles, 18 sites are fixed. On the three moving cycles, the nonconstant Fourier part of \(C\) is square-zero and has rank at most one at each root. Trace zero forces the three internal circulant sign blocks to be the same \(g\), as for order 19. Name the three upper cross blocks \(f,h,e\), in positions 12,13,23. Work directly modulo \(\Phi_{13}\); it need not be irreducible. At every primitive root the two-by-two minors imply equal cyclic autocorrelations for \(g,f,h,e\) modulo 3. Equality modulo the cyclotomic polynomial implies equality of the full cyclic autocorrelation because lag zero is identical. Together with these norm constraints, test the three off-diagonal entries of the nonconstant square-zero equation:
\[
 \bar f h=ge,\qquad fe=gh,\qquad h\bar e=gf\pmod{\Phi_{13}} .
\]
No division by \(g\) is used, so the zero-frequency-component cases of \(g\) in a factor of \(\Phi_{13}\) are retained. Enumerate the 64 symmetric diagonal blocks, independently rotate/switch \(f,h\) by choosing representatives (also \(f\le h\) by swapping the last two orbits), and allow all \(e\) in the same norm bucket. This gives 8324 tuples; 353 pass these equations, and 96 meet all within-orbit route counts (there is no fixed-coordinate difference within an orbit).

For each pair of moving orbits, the Hamming contribution from the 18 fixed sites is constant across relative cyclic shifts. Test each possible contribution and retain it only if all 13 shifts give an admissible route count. Normalize the first moving-orbit sign of each fixed row to \(+1\), giving four profiles \(p_r=(1,\pm1,\pm1)\). The three possible fixed-coordinate distances are determined by the four profile multiplicities \(n_r\) (sum 18). This leaves 3286 profile/tuple descriptions.

For any profile class of size \(n\), let \(k_i\) be the number of positives in row \(i\) of the fixed block at targets in that class. By symmetry, the total of fixed-coordinate distances over pairs in that profile is \(\sum_i k_i(n-k_i)\). Yet a negative pair needs at least 9 differences, and a positive pair at least 12. Since the diagonal is negative, each positive pair is counted twice in the internal \(k_i\)'s. Therefore
\[
 9\binom n2
 \le \sum_i k_i(n-k_i)-\tfrac32\sum_{i\text{ in class}} k_i
 \le(18-n)n^2/4+n(n-3/2)^2/4 .
\]
This forces \(n\le6\), leaving 552 descriptions.

Now let \(A\) be the unknown 18-by-18 fixed sign block and \(P\) its profile rows. Let \(M\) be the symmetric three-by-three matrix of moving-block sign sums modulo 3. The fixed-to-moving equation in \(C^2=0\) reads \(A P=-P M\) since 13 is 1 modulo 3. For a row of profile \(p_r\), the allowed per-profile positive counts \(k_s\) must obey
\[
 \sum_s(2k_s-n_s)p_s=-p_r M,\qquad 0\le k_s\le n_s-[r=s].
\]
Enumerate these four small integer counts. For any nonempty subset of profiles sum their distance inequalities, then upper bound the right side by maximizing each fixed row's joint contribution over the allowed \(k_s\). This discards all but **17** descriptions, falling into just these two kinds after switching moving signs and relabeling profiles:

1. Three classes of size 6, with \(M=0\).
2. Profile sizes \(3,5,5,5\), where \(n_{+++}=3\) and \(M\) has diagonal 2 and every off-diagonal 1 over \(\mathbf F_3\).

## 3. Close both patterns by parity, not by searching all fixed adjacency matrices

For the first kind the three occupied profiles are independent over \(\mathbf F_3\). The fixed-to-moving equation forces every row sum within each class to be zero modulo 3, hence \(k_s\in\{0,3,6\}\). For a class of six the distance bound attains its required value only if every column into that class has \(k=3\): the twelve external columns can contribute at most \(12\cdot9\); the six internal columns, with the penalty, at most \(6\cdot4.5\); their sum 135 exactly matches \(9\binom62\). Thus every fixed row has 9 positives across the three classes. All fixed-row distances are even. But there are negative pairs within a class (each row has only three positive neighbors among five), and equality requires their distance be exactly 9, impossible.

For the other kind, record, for each profile, the split into even and odd total fixed-row positive counts. Same-parity rows of one profile are at least distance 12, since their distance is both even and divisible by three. In a class of size 5, all five cannot have one parity: the distance sum would be at least \(12\binom52=120\), but at most \(18\lfloor25/4\rfloor=108\). Thus there are just \(4\cdot4^3=256\) splits to check. For opposite-parity rows of the same profile the distance is at least 9 for a negative pair and 15 for a positive pair. For one profile with parity counts \(n_0,n_1\), this requires total distance at least
\[
 12\left(\binom{n_0}2+\binom{n_1}2\right)+9\,n_0 n_1
 +6\,(\text{positive opposite-parity pairs}).
\]
For a row of that class and parity group size \(m\), at least \(\max(0,k+1-m)\) of its \(k\) internal positives have opposite parity, because the row cannot count itself. Charge three per directed incidence to pay the final six-per-pair term. As above, for each subset of profile classes, upper bound the resulting joint distance score by its allowed row-wise maximum, now with prescribed total parity \(\sum_s k_s\bmod2\). The **256** splits all fail. This last inequality uses only the integer profile degrees and parity, so no numerical or heuristic solver is involved.

This closes all order-13 orbit patterns. In the smaller primes 5, 7 and 11 only the one- and two-cycle cases are closed by this proof. The main Gewirtz transport question remains open with no symmetry imposed.

Replay:

    c++ -O3 -std=c++20 order13-signs.cpp -o order13-signs
    ./order13-signs
