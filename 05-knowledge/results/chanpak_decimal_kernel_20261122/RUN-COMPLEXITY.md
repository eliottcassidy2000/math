# How many binary runs does an absolute Farey crossing need?

**PROVED.** Fix any two finite nonempty tagged flanks selected by the unit Farey determinant around the exponent-15 middle. If exact crossings of the Euclid centers have both exterior bit lengths tending to infinity, the product of their numbers of runs of ones tends to infinity. For exterior lengths both proportional to \(p\), that product is \(\Omega(p/\log p)\). Run lengths and positions may be chosen afresh for each target. In each of the eight canonical spines, the hypothesis that both lengths grow can be dropped. For any fixed bound on the runs in each address, every sufficiently large Euclid center lies at positive **relative** distance from the accessible heights of that chart.

The estimates allow arbitrary fixed neighboring endpoint fractions and arbitrary fixed tagged realizations of those fractions in the two-growing-flank result. Their constants depend on the chart. An exponent-dependent choice of unbounded core sizes and the individual exponent-61 crossing remain open.

## The odd factor in the leading coefficient

Write the lower chart as
\[
 H=k+tx+cy+qxy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y)).
\]
All entries are nonnegative, \(h\) is odd, \(c,t>0\), and \(0\le f(z)\le z\). The specified physical middle has \(q=q_A\,320\,q_V\). Thus \(q\) has an odd divisor greater than one for every choice of endpoint tags and cores. This survives selecting both endpoints by their Farey determinant. [Chart and join](E15-FAREY.md) · [Endpoint selector](TWO-BOUNDARIES.md)

For a positive address \(x\) of width \(I\), each run of ones contributes a difference of two powers of two to \(x/2^I\), with exponents nonpositive. For instance \(110011_2/64=(1-2^{-2})+(2^{-4}-2^{-6})\). Every limit of such normalized addresses with boundedly many runs is dyadic. Pass to a subsequence on which each of the finitely many run endpoints either stays constant or tends to infinity. The terms in the latter group vanish. All limits lie in \([1/2,1]\).

If \(I,J\to\infty\), every term in \(H/2^{I+J}\) outside \(q(x/2^I)(y/2^J)\) tends to zero. Suppose along exponents \(p\to\infty\) we had \(H/E_p\to1\). Since \(E_p=2^{2p-1}(1-2^{-p})\), the integers \(2p-1-I-J\) are bounded. On a subsequence they equal a fixed \(j\). Bounded runs would force \(2^j=q\alpha\beta\) with \(\alpha,\beta\) positive dyadic rationals. The odd divisor of \(q\) makes this impossible in \(\mathbb Z[1/2]\). This also applies to any sequence of heights asymptotic to increasing powers of two. The earlier rational-slope obstruction treated a prescribed affine scaling of the addresses. [Earlier unequal-flank inverse](../chanpak_absolute_liouville_leech_20261122/UNEQUAL-FLANK-INVERSE.md)

## An explicit lower bound with arbitrary run positions

**PROVED.** Suppose both addresses are positive and \(H=E_p-1\). Let \(R_x,R_y\) count their runs of ones, \(K=4R_xR_y\), \(m=\min(I,J)\), and \(D=q+c+t+2h+e+r+k+1\). Then
\[
 m\le 1+\lceil\log_2(2D)\rceil+
             K\lceil\log_2(2qK)\rceil .
\]

Set \(s=I+J\), \(j=2p-1-s\), and \(G=\lceil\log_2(2qK)\rceil\). Here \(2^j\ge q/4>1\) because \(H\ge qxy\ge q2^{s-2}\). In the expansion
\[
 S=q(x/2^I)(y/2^J)-2^j
     =-2^j+q\sum_{\nu=1}^{K}\pm2^{-d_\nu}
\]
all \(d_\nu\) are nonnegative integers. Each partial sum retaining the terms with \(d_\nu\le d\) has magnitude at least \(2^{-d}\). After multiplication by \(2^d\), its numerator is nonzero modulo an odd divisor of \(q\). Initially take \(d=0\). If the next occupied exponent exceeds \(d+G\), the remaining tail has absolute value at most \(qK2^{-d-G}\le2^{-d-1}\). Otherwise advance \(d\) to that exponent. There are at most \(K\) advances and \(d\le KG\) throughout. Whether a large gap occurs or the terms run out, \(|S|\ge2^{-KG-1}\).

In the other direction \(E_p=H+1\le D2^s\) and \(2^{p-1}\le\sqrt{E_p}\). All chart terms except \(qxy\), together with the constant 1, have sum at most \(D2^{s-m}\). Hence \(|S|\le D2^{-m}+\sqrt D\,2^{-s/2}\le2D2^{-m}\). Combining the two estimates proves the inequality. For \(m\) proportional to \(p\) it implies \(R_xR_y=\Omega(p/\log p)\). The estimate permits either one of the run counts to stay small.

## When one address stays bounded

**PROVED.** For any chart as above, if both run counts are bounded, \(x\) is fixed, \(y\to\infty\), and \(H/E_p\to1\), at least one of \(c+qx,c+h+qx\) must be a power of two. For fixed \(y\) and growing \(x\), the alternatives are \(t+qy,t+h+qy\).

Along a subsequence \(y/2^J\to\alpha\in[1/2,1]\) and \(f(y)/y\to\beta\in[0,1]\). Both limits are dyadic. For the latter, the reversed-digit formula is \(f(y)/y=1-\sum_i y_i2^{-i-1}+O(J/y)\), where \(y_i\) is the bit of weight \(2^i\). Reversal preserves the run count, so the previous compactness argument applies. [Reversed phase formula](REVERSED-PHASE.md)

After passing to a further subsequence, \(2p-1-J\) is fixed and \(2^{2p-1-J}=\alpha(c+qx+h\beta)\). Each positive factor belongs to \(\mathbb Z[1/2]\). The positive units in this ring are the integer powers of two. Thus the second factor is a power of two, necessarily an integer since \(c\ge1\). Consequently \(h\beta\) is an integer. The oddness of \(h\) forces the dyadic \(\beta\) to be an integer as well, leaving only 0 and 1.

**FINITE-EXACT inputs to the all-exponent corollary.** In all eight canonical mask-12 charts, the largest powers of two dividing \(c,c+h,t,t+h\) are \(2,1,16,1\), respectively. The corresponding powers dividing \(q\), in order of spine \(4,8,\ldots,32\), are \(65536,131072,65536,262144,65536,131072,65536,524288\). Each of the four slopes exceeds its own power-of-two divisor. Adding a nonnegative multiple of \(q\) preserves its two-adic valuation, excluding a power of two. Any sequence with \(H/E_p\to1\) has a subsequence on which each address is fixed or tends to infinity. At least one tends to infinity. The two-growing-flank argument and the two bounded-coordinate tests exclude all cases. Thus for each fixed run budget on these charts some \(\epsilon>0\) gives \(|H/E_p-1|\ge\epsilon\) for all sufficiently large \(p\) and all such address pairs.

## A finite bound for two periodic substitutions

The [mixed-word drift theorem](../chanpak_exponent31_seams_20261122/README.md) already excludes infinitely many exact crossings when at least one growing mixed periodic block occurs on a fixed affine repetition schedule. For \(x=U W^{an+b}V\), \(y=U'(W')^{a'n+b'}V'\) there is a bound on the number of coincidences. Take \(p=p_0+dn\), with \(d>0\), and omit indices before every growing repeated block has appeared at least once. The difference \(H-(E_p-1)\) is a rational linear combination of \(1,n\) and at most five exponentials \(2^{sn}\) with distinct positive \(s\) in
\[
       \{a|W|,\ a'|W'|,\ a|W|+a'|W'|,\ d,\ 2d\}.
\]
The drift theorem supplies a positive linear coefficient if a growing block is mixed. If both addresses grow through pure repeated blocks, the odd factor in \(q\) excludes the identically-zero function. The bounded-coordinate slope test does the same for one growing pure block on a canonical chart. Any nonzero function of this form has at most six real zeros. Differentiating twice removes the affine part, and a sum of \(k\) distinct real exponentials has at most \(k-1\) zeros by Rolle induction after dividing by its smallest exponential. The purely affine case has at most one zero. If both address bit rates equal \(d\), at most three indices remain possible.

From this directory run the following for the eight chart certificates, exhaustive small sparse-product tests and examples of the reversal estimate.

~~~sh
python3 -B run_complexity.py
~~~
