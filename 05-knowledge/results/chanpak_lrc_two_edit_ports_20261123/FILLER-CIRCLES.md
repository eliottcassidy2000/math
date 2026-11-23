# Lonely Runner circles from literal third-slot filler axes

**PROVED** for the marker and covering transfers; **FINITE-EXACT** for the arithmetic certificates below. This uses the [14-pair literal bank](THIRD-SLOT-JACOBIAN.md), the [earlier modular wheel](WHEEL.md), the Lonely Runner [prime-colour decoder, THM-3825](../../../01-canon/theorems/THM-3825-prime-colour-valuation-two-cube-decoder.md), and only the proved max-gap mechanism from the [corrected anchor-floor note S153](../../../07-reflections/the-anchor-floor-turns-the-maxgap-crux-into-finite-order-statistics-klein-S153.md). It does not use the refuted AP-minimality or one-anchor uniform claims.

## The 14 directions are labels, not a proof about 14 runners

Each pair has readings \((D_i,Q_i,R_i^+),(D_i,Q_i,R_i^-)\). Set \(r_i=R_i^+/R_i^-\); the private-marker array is exactly \(v_{p_i}(r_j)=\delta_{ij}\). This is a diagonal version of the prime-colour sidecar: valuations recover which slot choices were made even after forgetting the factors' order. For a fixed quota \(t_i\) of each type, ordinal composition gives
\[
 R=R_0\prod_i r_i^{k_i},\quad 0\le k_i\le t_i,\qquad
 h=C+20R,\quad C=11D+42Q,\quad d=-3D-6Q.
\]
\(D,Q,d\) remain fixed and both edit orders still lead exactly to \(7h\). The exact relation lattice \(\{z:\prod_i r_i^{z_i}=1\}\) is zero by the private valuations. This is **not** the relation lattice \(v^\perp\cap\mathbb Z^{13}\) of integer speeds used by [LRC(14)'s covering transference, THM-4009](../../../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md). A modular projection below produces that kind of resonance, at the cost of losing absolute information.

**Important arithmetic guard.** In the full equal-quota bank (\(t\ge1\)), the common divisor of all \(D,Q,R\) across the grid is exactly \((3^6 5^4 11)^t=5{,}011{,}875^t\). Hence every raw start in that full grid is composite. Seven is the maximum subbank size whose profiles are all primitive, at every positive equal quota and every multiplicity choice. Exactly two seven-subsets work, with 1-based bank numbers \(\{1,2,3,4,7,8,12\}\), \(\{1,2,3,4,8,9,12\}\). To certify this, intersect the union of primes in selected \(D_i\), the union for \(Q_i\), and the union for the two alternatives of \(R_i\). A prime in the intersection can be made to divide all readings by making a suitable \(R\) choice; an empty intersection precludes this at any choice. The [exact replay](filler-circle-check.py) scans all subsets. Primitive input does not guarantee a prime output.

## 1. A real max-gap selector at fixed Jacobian

Take the second primitive seven-subbank. Orient each ratio to exceed one, and use \(b=302043/301803\) from pair 2 as the main clock. For each combination of zero, one or two switches in each of the other six axes, take its product \(A_j\). Find \(c_j\) such that \(\bar A_j=A_jb^{-c_j}\in[1,b)\). The real phase is \(\log_b\bar A_j\) on a circle. Sorting these is an exact finite version of the S153 max-gap move: the largest uncovered interval governs worst-case rounding, rather than average density.

The 729 phases are distinct and include 1 in the multiplicative representation. The maximum \(c_j\) is 83. The largest multiplicative gap, including the wrap to \(b\), is \(G<1+4338/10^9\), certified using rational comparisons only. At each common quota \(t\ge84\), put \(X'= (X-C)/(20R_0)\) where \(R_0\) now uses the smaller member of every pair. For any target \(b^{83}\le X'\le b^t\), set \(q=\lfloor\log_b X'\rfloor\) and take the largest \(\bar A_j\le X'b^{-q}\). The main exponent \(q-c_j\) lies in \([0,t]\) and the other exponents lie in \([0,2]\). Thus the physical start \(\widehat h=C+20R_0 b^{q-c_j}A_j\) satisfies
\[
 0\le X-\widehat h < 4.338\cdot10^{-6}(X-C).
\]
Every such candidate has the same \(D,Q,d\) and order \(63t+10\) at that quota. By using more anchors, the worst relative error can be made arbitrarily small on analogous intervals: any ratio outside the clock axis has irrational logarithm in base \(b\) by the private valuations, so its nonnegative powers have dense phases. For each fixed desired error take a finite prefix before increasing the quota. This is a real *relative* selector on a specified interval, not a shrinking absolute error or an exact prescribed-integer selector.

## 2. A single circular clock: literal runner phases

Use the original ordered ratios for the full bank. Give every type \(N\) copies, let \(h_0\) use all minus choices, and let \(h_i(k)\) switch only axis \(i\) exactly \(k\) times. Then \((h_i(k)-C)/(h_0-C)=r_i^k\). This observable is a quotient of physical heights after removing a known offset, not another count. Clock steps replace factors; they are not the two marked external reversals.

At \(\ell=977,\ m=61,\ \zeta=77\), the character \(u\mapsto u^{16}\) has logs in base \(\zeta\)
\[
 v=(17,15,48,49,1,43,50,33,30,42,40,52,26,10).
\]
Applying this character to the literal quotient realizes \(kv_i\bmod61\). These are fourteen distinct nonzero sign classes. The exact max-min circular distance over all 61 ticks is \(5/61\), attained at tick 9. This exceeds \(1/15\) for these fourteen moving runners plus a stationary one. The repo's LRC(14) convention is thirteen moving runners plus a stationary one; deleting a direction also retains a witness here. More generally all displayed \(v_i\) are units, so using \(k_i(s)=[v_i^{-1}(a_i s+b_i)]_{61}\) realizes any desired finite speeds and offsets at this field with quota 60, though these calibrated exponents need not advance by one.

For a raw common exponent without resets, the markers give a stronger existence theorem: for every \(m\ge2\) coprime to all markers and any \(a\in(\mathbb Z/m)^{14}\), there are infinitely many primes \(\ell\equiv1\bmod m\) and primitive roots \(\zeta\in\mathbb F_\ell\) of order \(m\) with \(r_i^{(\ell-1)/m}=\zeta^{a_i}\) for every \(i\). Over \(K=\mathbb Q(\zeta_m)\) the private markers are unramified, so their valuations prove independence of the 14 classes in \(K^\times/(K^\times)^m\). The [Kummer pairing](https://www.jmilne.org/math/CourseNotes/FT421.pdf) gives all \(m^{14}\) radical automorphisms, with the chosen one sending the \(i\)-th radical to \(\zeta_m^{a_i}\) times itself. The radical field is Galois over \(\mathbb Q\). [Chebotarev](https://doi.org/10.1007/BF01206606) supplies rational Frobenius classes of the specified element; its trivial cyclotomic restriction forces \(\ell\equiv1\bmod m\). Choose a prime above \(\ell\) with this Frobenius representative. Reduction of the radical quotient is \(r_i^{(\ell-1)/m}\), and reduction of \(\zeta_m\) supplies the character generator. Exclude the finitely many nonunit primes for the literal quotient.

For integer speeds \(a_i\), the discrete grid maximum approximates the continuous loneliness within \(\max_i|a_i|/(2m)\), by taking a nearest grid point to a maximizing time. This shows a faithful finite simulation for arbitrarily fine grids. It does not prove a uniform lower bound for the simulated clocks; the hard universal Lonely Runner problem is unchanged.

## 3. Exact modular wheels with physical counts

At the 14 private primes choose all quotas equal to a positive multiple of \(\mathrm{lcm}(p_i-1)\) and switch all copies of a type, or none. Then \(D\equiv Q\equiv1\) at each marker, and the vector of \(R\) residues is any desired vector of zeros and ones. Each switch affects only its own marker. The four actual heights reduce to \((53,329,113,371)\) or \((73,445,157,511)\), respectively; all are units at every private marker. Thus all these choices pass that fixed wheel with the same Jacobian, despite the **forced** small-prime divisor outside the wheel. This is a hostile case for substituting modular coverage for absolute prime coverage.

For one variable exponent at primes where all readings are units, a forbidden affine reading of \(R\) rules out at most one clock residue mod \(\operatorname{ord}_{\ell}(r_i)\) (unless it vanishes identically). Across a common full clock cycle the elementary union bound works without independence: the bad fraction is at most the sum of these reciprocal periods, with multiplicities. This is the same avoidance step as the earlier [physical wheel](WHEEL.md). The more concrete joint **unit** wheels and a prime-producing subbank are in [the CRT companion](FILLER-CRT.md).

The transferable mechanism is finite phase control at fixed \(D,Q,d\), and an actual covering estimate for the absolute height after subtracting the offset. Neither encodes all integer targets or resolves Chan–Pak or Lonely Runner.
