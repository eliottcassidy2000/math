# First absolute noncube overlaps in the short-guard chart

This continues the [literal parity-cube bridge](../pollock_chanpak_literal_parity_cube_20261126/README.md) and the [small-edit obstructions](../pollock_chanpak_small_edit_rigidity_20261127/README.md). The Pollock target is \(N=\sum_i F(x_i)\), \(F(x)=(x^3-x)/6\), with all \(x_i\ge 2\). The common readout with Chan–Pak is \(H=2N+1\). Here the equality is over the integers, not just a residue.

## 1. Four exact matches and a first-in-chart certificate

These are **all** the matches through upper-endpoint height 17,023,575 between a noncube four-active unit step preserving the first two moments and the [one-site-guard six-vertex fork](../moore_fixed_order_spectra_20261124/TWO-COORDINATE-RIGIDITY.md) with a 1010-leading payload \(v\). No other tournament charts are excluded.

| Lower and upper tournament counts | Payload \(v\) | Five lower indices | Five upper indices |
|---|---:|---|---|
| 8,541,511 → 8,541,513 | 5185 | 281, 25, 31, 145, 70 | 282, 12, 48, 140, 70 |
| 16,818,421 → 16,818,423 | 10259 | 336, 80, 86, 200, 150 | 337, 67, 103, 195, 150 |
| 17,021,085 → 17,021,087 | 10415 | 20, 290, 254, 162, 182 | 21, 277, 271, 157, 182 |
| 17,023,573 → 17,023,575 | 10337 | 327, 55, 111, 237, 108 | 338, 38, 170, 184, 108 |

In particular the first common noncube target is \(N=4,270,755\) within this chart. The third row uses the negative-product orientation of the same triple shape as the first two.

**Actual tournament.** Give a transitive prefix the apex bits of the binary word of \(2v\), using 1 for an arc out of the apex. Follow it by a six-vertex block \(Y\), a guard vertex, and a last vertex, in that order; all cross-block core arrows go forward. Internal backwards edges of \(Y\) have mask 177 or 181 in lexicographic order of pairs, least significant bit first. The apex beats \(Y\)-slots 0, 1, 4 and the guard, and loses to the last vertex. The only changed arc is \(Y\)-pair (0,3). The apex reaches all earlier-to-later blocks via the first prefix vertex; the final vertex can be reached from every core vertex and links back to the apex. Both tournaments are strong. The replay independently enumerates paths by splitting at the apex and enumerating Hamiltonian orders on the two induced subsets of \(Y\). It also checks the exact inherited index formula
\[
H(v)=1662v-58q(v)+427,\qquad
q(2v)=q(v),\quad q(2v+1)=q(v)+\lceil v/2\rceil,\quad q(0)=0.
\]

## 2. The Walsh reduction and why the search is complete

Ignore the retained fifth index. For lower and upper four-index lists \(A,B\) of equal first and second moments, form \(P(Z)=\sum_{b\in B}Z^b-\sum_{a\in A}Z^a\). If the \(F\)-gap is one, then \(P(Z)=(Z-1)^3R(Z)\) with \(R(1)=1\): the third derivative at one is \(6\). Consequently \(P(-1)\equiv8\bmod16\), and \(|P(-1)|\le8\). Equality forces all four lower indices to have one parity, and all four upper indices the opposite parity.

Write \(S\) for their common sum, \(U,V,W\) for the three pair differences of the four-coordinate Walsh transform. Each difference is \(S\bmod4\), since any pair of indices has even sum. The \(F\)-gap equals the Walsh product difference divided by 16. Thus \(S\) cannot be divisible by four, which would make both products divisible by 64. Setting \(s=S/2\) and \((u,v,w)=(U,V,W)/2\) gives two odd triples of equal squared norm with signed products differing by two. Conversely, invert by
\[
\frac12(s+u+v+w,\ s+u-v-w,\ s-u+v-w,\ s-u-v+w).
\]
The squared norm determines the second moment at a given \(s\); the \(F\)-gap is half the product difference. Permuting the triple or changing two signs only permutes the four indices. Opposite product signs force both magnitudes to be one, the original parity cube. For equal signs, sort positive odd magnitudes and test both sign orientations. In fact noncube triples cannot share a magnitude: a common divisor of two odd products differing by two must be 1; removing that 1 gives two pairs of equal squared norm whose products differ by two, hence whose even coordinate sums have squares differing by four. That is impossible for positive even sums.

For a complete bounded search, put \(B_0=371\), the largest possible index through the stated cutoff. For either four-index list, \(p_3\le3(H_{\max}-1)+4B_0\). The Walsh triple norm \(q_0\) satisfies
\[
q_0=p_2-\tfrac14 p_1^2\le\tfrac34p_2,
\qquad 16q_0^3\le27p_3^2
\]
by positivity and Hölder. Thus \(q_0\le163884\). The exhaustive replay finds only two noncube shapes in the entire norm range, including both minima:
\[
(65,71,185),\ (53,89,181)\quad(q_0=43491);
\qquad
(17,73,199),\ (11,143,157)\quad(q_0=45219).
\]
It tries every odd \(s\le2B_0\), both signs, and every positive common fifth index within the cutoff. Independently, \(H(v)\ge1604v+427\) from the guard formula gives a finite payload bound; the replay scans every 1010-leading \(v\) in it. The intersection gives exactly the four rows.

## 3. Rational-function obstruction to a cubic join lift

An ordered join multiplies the path count, so the tournament side can multiply both endpoints by \(k^3\) whenever that factor is itself a realizable height. A suggested Pollock analogue scales the four centered offsets by \(k\) and sets their center \(c(k)\) and common fifth index \(z(k)\). Even allowing **arbitrary rational functions** \(c(k),z(k)\in\mathbb Q_2(k)\), no such identity can hold at infinitely many \(k\).

Indeed \(2\sum F(x)+1 = 1+\tfrac13\sum(x^3-x)\). For fixed centered offsets \(a_1,\ldots,a_4\), expanding the four active terms gives the rational-function identity that would be required:
\[
4(c^3-c)+(z^3-z)+3k^2c\sum_i a_i^2+k^3\sum_i a_i^3
=3k^3H-3.
\]
If one of \(c,z\) has a pole at \(k=0\), the lowest power must come from a cubic term. Unequal pole orders cannot cancel. With equal pole orders the leading coefficients \(C,Z\in\mathbb Q_2^\times\) would satisfy \(4C^3+Z^3=0\), impossible because \(2+3\nu_2(C)\) and \(3\nu_2(Z)\) cannot agree. Both functions are therefore regular at zero. Evaluating there leaves \(4(c(0)^3-c(0))+z(0)^3-z(0)=-3\). If either value has negative 2-adic valuation, the negative valuations of the two terms again cannot coincide, so the sum cannot have valuation zero. If both values are integral, the sum is even. This proves the obstruction, and it already applies on just one endpoint.

The four exact overlaps survive because they need not scale all offsets together. Neither the first-in-chart search nor the rational-function obstruction settles Pollock or Chan–Pak.

Replay: run the two JavaScript files with Node. The scan uses integer arithmetic in the stated range; the norm ceiling is evaluated with BigInt. The path replay does not use the guard formula.
