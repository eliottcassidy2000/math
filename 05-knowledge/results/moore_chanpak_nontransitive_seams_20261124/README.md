# Nontransitive fillers in the coupled pentagon

This extends [the coupled pentagon rail](../moore_chanpak_counts_20261124/moore-coupled-rails.md). Equal transitive widths are sufficient and necessary **among pairs of chains**, not among all fillers.

## An exact nontransitive unit seam

Use the prepared seven-site tournament with out-neighbor masks
\((122,100,25,82,34,76,20)\), substitute fillers \(X,Y\) in positions 3 and 4 (zero based), then reverse the physically unchanged exterior contact \((5,6)\). The first certified witness is
\[
 X=(66,284,409,209,33,463,150,275,89),\quad
 Y=(62,124,56,112,32,64,21).
\]
Its 21-vertex counts are **19,116,233,339 → 19,116,233,337**. Both fillers are rigid and nontransitive and have different orders. A full subset Hamilton-path DP on the two 21-vertex tournaments [replays the counts](literal.cpp) independently of the transfer formula.

The complete census through order eight finds 1,923 distinct triples \((h,Q,R)=(Z_1,Z_2,Z_3)\) and no unit pair with a nontransitive filler. Thus **9 is the smallest possible maximum filler order for this particular chart**. The two additional 9+9 pairs in [certificates.py](certificates.py) are exact witnesses but do not constitute an order-nine census or a minimum combined-order proof.

## The arithmetic mechanism

For a filler let \(\rho=(h,Q-2h,R-3Q+3h)\), the ordered counts with exactly one, two or three nonempty path slots. The two endpoint heights are \(\rho_X^{\mathsf T}M\rho_Y\) and \(\rho_X^{\mathsf T}M'\rho_Y\), where direct enumeration of the exterior words gives
\[
 M=\begin{pmatrix}109&117&32\\94&117&40\\24&36&16\end{pmatrix},\qquad
 M'=\begin{pmatrix}111&100&24\\111&117&36\\32&40&16\end{pmatrix}.
\]
This recovers the signed gain
\[
 \Delta=2h_Xh_Y+5(Q_Xh_Y-h_XQ_Y)+4(R_XQ_Y-Q_XR_Y).
\]
In particular the gain compares \(\theta=(4R-5h)/Q\) with an exact rational tolerance. Both fillers must have primitive profile \(\gcd(h,Q,R)=1\): the content of either side's coefficient vector in the gain equals the other side's profile gcd. For example the \(X\)-coefficient vector is \((2h_Y-5Q_Y,\,5h_Y-4R_Y,\,4Q_Y)\); its content is odd (the middle entry is odd), divides \(Q_Y\), then \(h_Y\), then \(R_Y\), and conversely. This has a symmetry interpretation. A filler automorphism extends to both endpoint tournaments, acting freely on each Hamilton word, so its group divides both heights and the gain. **An exterior unit seam cannot leave any nontrivial internal filler symmetry intact**, since tournament automorphism groups have odd order. This argument applies to any exterior flip of an ordinary substitution chart, not just the pentagon.

**Exact finite search if either filler is fixed.** For a tournament \(T\) of order \(n\), \(Q/h\ge n+1\) because there are \(n+1\) distinct prefix cuts of every Hamilton word. Also \(R/Q\ge 1+n/2\): split the first path of an ordered two-cover at any of its cuts, including the ends; the resulting ordered three-cover recovers the two-cover and the cut. The mean length of the first path is \(n/2\), by interchanging the colors. If \(X\) is fixed and \(|\Delta|=2\), it follows that the order \(n\) of \(Y\) satisfies
\[
 2n+4\le
 \frac{4R_X+5h_X}{Q_X}
  +\frac{5+2(h_X+1)/Q_X}{n+1}.
\]
Indeed divide the gain identity by \(Q_XQ_Y\), bound \(h_Y/Q_Y\le1/(n+1)\) and use \(h_Y\ge1\). The same inequality holds interchanging the names, since the estimate uses absolute values. The right side decreases and the left increases, giving an explicit finite search bound independent of the unknown filler's structure. For either filler in the 21-vertex witness the partner order is at most 15. This rules out an infinite unit family even if only **one** filler is kept fixed while the other is allowed to vary arbitrarily. All-orders equal-chain families evade it because both grow.

For a prescribed pair \((N,N+\Delta)\) with \(Y\) fixed, the two matrices also give an integer-lattice test before enumerating fillers. Put \(w=M\rho_Y\), \(d=(M'-M)\rho_Y\), and use columns \(v_i=(w_i,d_i)\); rank two holds here, and the lattice index is the gcd of the three \(2\times2\) minors. Then \((N,\Delta)\) is in the integer image iff adjoining this vector does not lower that gcd. For unit witnesses it is appropriate to also check the finite order bound and demand an actual nonnegative path profile. For the three exhibited fixed right fillers the lattice index is one, so this particular congruence test does not by itself select the witnesses.

## Prime filters

Every exponent of \(Z_k(T)\) has the parity of \(|T|\), by the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569). Thus \(Q\equiv(-1)^{|T|}h\bmod3\) and \(R\equiv0\bmod3\). Modulo three the two endpoint heights are \(h_Y(Q_X+2h_X)\) and \(h_X(Q_Y+h_Y)\). They can both avoid divisibility by three only if \(X\) has odd order, \(Y\) even order, and \(p=h_Xh_Y\) is prime to three. In that case the base height is \(p\bmod3\), the new height is \(-p\bmod3\), and \(\Delta\equiv p\bmod3\). Equal-width transitive fillers for this chart never enter that lane.

The unit condition excludes \(X\) singleton or \(Y\) of order two in the viable lane. For a singleton \(X\), \(\Delta=12h_Y+7Q_Y-8R_Y\le-5Q_Y\) using \(n_Y\ge2\) in the two cut bounds. For \(Y\) of order two, \(\Delta=16R_X-31Q_X-18h_X>2\) using \(n_X\ge3\). Hence the viable sizes are odd \(n_X\ge3\) and even \(n_Y\ge4\). Now \(Q_X\equiv4C(X)\bmod8\), where \(C\) counts directed Hamilton cycles up to rotation, and \(Q_Y\equiv0\bmod8\): in the odd-cycle expansion for \(Z_2\) the only term of valuation below three in these sizes is a spanning odd cycle, contributing four. In \((h,Q,R)\) coordinates the first height matrix is \(\bigl(\begin{smallmatrix}11&-9&0\\-20&33&-8\\0&-12&16\end{smallmatrix}\bigr)\). Thus the first height is \(3p\bmod8\) and \(\Delta/2\equiv p+2C(X)\bmod4\). Together with the mod-three relations and \(\Delta=\pm2\), the **smaller** endpoint is \(11\) or \(23\bmod24\) if \(C(X)\) is even, and \(5\) or \(17\bmod24\) if it is odd. These are filters, not sufficient criteria for a twin-prime pair.

The companion [catalog-pressure result](catalog-pressure.md) shows why collecting fixed-height witnesses, or even allowing some catalogs of unbounded template order, still does not give a Chan–Pak selector.

A genuine prime-to-composite seam in this lane is
\[
 \mathbf{85{,}587{,}188{,}881}\ \longrightarrow\
 85{,}587{,}188{,}879=37\cdot167\cdot2099\cdot6599.
\]
It uses order-nine \(X\) and order-eight \(Y\), both rigid and nontransitive; their exact masks and profiles are in [certificates.py](certificates.py). The prime claim is from complete trial division to the integer square root, not a probable-prime test. Both heights also replay literally on the 22-vertex host.

The exhaustive parity-lane search with odd \(X\) through order 9 and even \(Y\) through order 8 checks 32,820 distinct left profiles and 1,103 right profiles, including the transitive fillers in those bounds. **Exactly one** pair has gain of magnitude two, the prime-to-composite pair above; hence there are no twin-prime seams in that order box for this chart. The parity filter cannot discard a twin prime involving 3 because every endpoint already contains the corresponding seven-site template as an induced subtournament and has at least its Hamilton count (extend paths by insertion). This is not a claim beyond the order bounds or about other coupled charts.

For a stronger formal filter with fixed \(u=h_X,v=h_Y,q=Q_X,r=Q_Y\), put \(D=(\Delta-2uv-5(qv-ur))/4\). In the viable lane the mandatory congruences are
\[
 R_X\equiv3^{|X|}u\pmod{48},\qquad R_Y\equiv3^{|Y|}v\pmod{144}.
\]
Indeed in \(Z_3\) every exponent has order parity; replacing \(3^{a}\) by \(3^{n}\) in a term is valid modulo 16 because their difference is divisible by 8 and nonempty packings have even coefficients. Positive exponents are at least one for odd order and two for even order. Consequently the linear equation \(rR_X-qR_Y=D\) has integral solutions respecting these two congruences exactly if \(D\) is integral and
\[
 D\equiv r3^{|X|}u-q3^{|Y|}v\pmod{\gcd(48r,144q)} .
\]
The modulus in the original gain test is at least \(768\). This solves only the formal congruence; it does not realize the corresponding profiles.

## Reproduce

- Run **python certificates.py** for the endpoint matrices, all four explicit witnesses, rigidity, trial division, and finite partner bounds.
- Compile **literal.cpp** with optimization and run it to verify the smallest exhibited host by literal 21-vertex subset DP (about 350 MiB). Pass any argument to use the 22-vertex prime-to-composite fixture (about 704 MiB).
- Run **python census.py** for the complete order-eight no-hit census. It constructs all ordinary classes from smaller ones by adjoining one vertex in every way and canonicalizing within score classes. Ordinary isomorphism preserves the profile.
- To reproduce the parity-lane order-nine scan, generate **n8.data** with **python census.py --rows > n8.data**. Compile **prime-lane.cpp** with C++17 and optimization and run **./prime-lane 9 0 n8.data**. Starting with all order-eight isoclasses and adjoining a vertex in all 256 patterns covers every order-nine isoclass; deduplication by exact profile loses no possibility for this chart.
