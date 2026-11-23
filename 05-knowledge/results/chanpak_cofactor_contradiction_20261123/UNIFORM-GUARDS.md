# Which modular guards survive arbitrary uniform contexts?

**PROVED scoped / FINITE-EXACT.** Companion to [the cofactor contradiction test](COFACTOR-TEST.md). The factor-seven fork depends on two slot readings, not on universal modular persistence. These are provably different requirements.

## 1. The integer observation lattice

Let \(M\) be a nonempty tournament of order \(n\). Write \(P_k\) for its number of ordered covers by \(k\) nonempty directed paths, and \(Q(t)=\sum_k \binom{t}{k}P_k\) for its ordered covers by \(t\) slots allowing empties. In a uniform exterior, each outside vertex beats all of \(M\) or loses to all of \(M\). Let \(g(M)\) be the gcd of the resulting Hamilton path counts, including the empty exterior. Then
\[
g(M)=\gcd_{1\le k\le n} (k-1)!\,P_k.
\]

Indeed, collapse the \(k\) maximal \(M\)-segments in a Hamilton path into symbols. The remaining data consist of a prefix, \(k-1\) nonempty exterior bridges and a suffix. Any internal ordered \(k\)-cover fits that exterior skeleton. The \(k-1\) bridges can be permuted freely: each starts at a vertex defeated by \(M\) and ends at a vertex beating \(M\). Since the bridges are nonempty and disjoint this action is free. Thus the coefficient of \(P_k\) is divisible by \((k-1)!\).

Conversely, define \(U_r\) by adding the exterior transitive order \(b_1,c_1,\ldots,b_r,c_r\), with \(M\to b_i\) and \(c_i\to M\). Between successive module visits each bridge uses at least one \(b\) and one \(c\). There are at most \(r+1\) module segments. At that maximum each bridge is \(b_i\to c_j\) with \(i\le j\), and all exterior vertices are used. Summing indices forces \(i=j\) on each bridge; the \(r!\) bridge orders all work. This makes the coefficient of \((k-1)!P_k\) a unit triangular matrix in \(H(U_0),\ldots,H(U_{n-1})\). In particular these actual reads generate the full integer observation lattice (integer combinations need not themselves be individual counts). For \(r\ge1\), \(U_r\) is strong: every outside pair \(b_i,c_i\) lies on a directed triangle with each vertex of \(M\).

For every prime power \(p^a\), the test can stop once \(v_p(k!(k-1)!)\ge a\): \(k!\mid P_k\) by free permutation of its ordered nonempty slots, so all such weighted terms vanish automatically. More explicitly, only the ladder reads \(U_{k-1}\) for \(1\le k\le n\) with \(v_p(k!(k-1)!)<a\) are needed. The test certifies the gcd even when the module has no corresponding automorphism.

## 2. The field test and two sharp order thresholds

By [the odd-cycle-packing expansion of Grinberg and Stanley](https://arxiv.org/abs/2307.05569), \(Q\in\mathbb Z[t]\) is monic of degree \(n\), has parity \(n\), and vanishes at zero. A disjoint packing of directed odd cycles contributes \(2^{c}t^{c+u}\), with \(c\) cycles and \(u\) unused singletons. For an odd prime \(p\),
\[
p\mid g(M)\quad\Longleftrightarrow\quad Q\bmod p\text{ is the zero function on }\mathbb F_p.
\]
In the lattice formula only \(k\le p\) can matter, with factorial invertible there. Evaluating \(Q\) at \(0,\ldots,p-1\) is triangular for \(P_k,\ k<p\); \(P_p\equiv0\) automatically by the factor \(p!\). Hence by parity only \(Q(1),\ldots,Q((p-1)/2)\) are needed. If \(p\) is erasable, one of these reads is nonzero, and the corresponding triangular change of basis guarantees an actual eraser among \(U_0,\ldots,U_{\min(n-1,(p-3)/2)}\).

It follows that universal \(p\)-persistence needs \(n\ge p\). Universal \(p^2\)-persistence needs **\(n\ge2p\)**. To see the second bound, suppose \(p\le n<2p\). All \(P_k\) with \(k\le p\) vanish modulo \(p^2\); all later ones are divisible by \(p\). In particular \(Q(t)\equiv0\bmod p^2\) for \(1\le t\le p\). Parity gives \(Q(-t)\equiv0\bmod p^2\) there as well, so division by \(p\) gives
\[
\sum_{j=1}^{n-p}(-1)^{p+j}\frac{P_{p+j}}p
 \binom{t+p+j-1}{p+j}\equiv0\bmod p.
\]
For \(1\le t\le p\), Lucas's digit formula reduces the binomial to \(\binom{t+j-1}{j}\bmod p\): if the low digit wraps, both sides vanish. These are polynomials of distinct degrees \(j\le p-1\) with invertible leading coefficients. Vanishing at all \(p\) residues forces every \(P_{p+j}/p\) to vanish modulo \(p\). If \(n>p\) this contradicts \(P_n=n!\); if \(n=p\), the imposed \(p^2\mid P_p\) already contradicts \(P_p=p!\).

A cyclic tournament on \(p\) vertices achieves order \(p\), by a rotation acting freely on full Hamilton paths in every uniform exterior. Two such modules in ordinal sum achieve order \(2p\) for \(p^2\), by their independent rotations. Placing them and a singleton in a directed triangle achieves the strong upper bound \(2p+1\); for \(p=3\) the strong order-6 example with rows \((48,49,3,7,12,28)\) and \(Q(t)=t^6+16t^4+28t^2\) already achieves the lower bound and has \(g=9\).

For a wider order bound, take \(1\le k<p\). Universal \(p^{2k}\)-persistence requires \(n\ge(k+1)p\). The raw top factorial already forces \(n>kp\). If \(kp<n<(k+1)p\), then \(p^{k+1}\mid P_i\) for \(i\le kp\) by the weighted test, and \(p^k\mid P_i\) for \(i>kp\) by \(i!\mid P_i\). For the base-\(p\) digits \(t_0,t_1\) of \(t\bmod p^2\), Lucas now gives
\[
 (Q(t)/p^k)\bmod p = \binom{t_1}{k}\,A(t_0),\qquad
 A(z)=\sum_{j=1}^{n-kp}(P_{kp+j}/p^k)\binom{z}{j}.
\]
Here \(A(0)=0\). For \(1\le t\le p-1\) the right side is zero; on replacing \(t\) by \(-t\), the upper digit is \(p-1\) and \(\binom{p-1}{k}\ne0\bmod p\). Parity therefore forces \(A\) to vanish at all residues. The top coefficient contradicts \(v_p(n!)=k\). At the boundary \(n=(k+1)p\) the same argument applies to the intermediate indices, since \(P_n=n!\) itself is divisible by \(p^{k+1}\). All \(P_i\) are then divisible by \(p^{k+1}\), so \(Q(t)\equiv0\bmod p^{k+1}\) for every integer \(t\). Taylor expansion at \(a\) on the \(k+1\) inputs \(a+jp\), \(0\le j\le k\), gives a Vandermonde system with determinant prime to \(p\) when the unknowns are the first \(k+1\) Taylor coefficients multiplied by their corresponding powers of \(p\). These coefficients consequently vanish modulo \(p\). Thus every residue is a root of multiplicity at least \(k+1\) and
\[
 Q(t)\equiv (t^p-t)^{k+1}\pmod p
\]
at the order boundary. This identity is only necessary for the original prime-power guard.

## 3. Cost in the actual contradiction statistic

For a strong persistent component, [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives \(H\ge L_p=5^{(p-1)/3}\). If all \(K_p\) components with \(p\mid H\) in a host of height \(x\) were required to be persistent, then
\[
p^{K_p}\ \le\ x^{3\log p/((p-1)\log 5)}.
\]
For \(p\ge 11\) this exponent is below \(1/2\): \(5^{10/3}>11^2\) and \(\log(p)/(p-1)\) decreases. Thus this approach cannot reach the square-root excess even with no cofactor loss. More generally, \(J\) persistent components among \(K_p\) cost \(x\ge p^{K_p-J}L_p^J\). The stronger order bound replaces \(L_p\) by \(5^{(2p-1)/3}\) for a component required to retain \(p^2\) in all contexts. In particular a strong realization of \(p\) or \(7p\) for \(p\ge11\) necessarily has a uniform eraser: \(7p\le p^2<L_p\). These are restrictions on *individual* SCCs, not on the product of slots of an ordinal sum.

The \(7\) case also has an exact cost. A strong module with both \(7\mid H\) and \(7\mid Q(2)\) has height at least **49**, and this is attained even for universal persistence. At order \(n<7\) the field test forbids persistence. Busch's exact minima start at 25 for \(n=7\), 45 for \(n=8\) and increase thereafter; among odd multiples of \(7\) below \(49\) only order \(7\), height \(35\), needs checking. The exhaustive anchored replay scans all \(2^{15}\) tournaments extending the fixed directed path on seven labels. Every tournament has a Hamilton path, so this covers all types. None of the 70 height-35 tilings passes even the \(Q(2)\) test, although all pass \(Q(3)\). At height 49, 98 tilings pass and 245 fail (tilings, not isomorphism classes). In particular \(7^{K_7}\le\sqrt{x}\) under this strategy, even when each counted SCC protects only the first two reads. The minimum 49 also holds without the strong restriction: 7 and 21 are absent, and at height 35 the only nonsingleton SCC has count 35 (otherwise a factor would have count 7). Singleton padding multiplies \(Q(2)\) by a power of 2 and cannot repair the guard.

## 4. Track only the fork's two reads

Vanishing of \(Q(1)\) and \(Q(2)\) is far cheaper than universal persistence. More generally the first \(R\) readings, with \(2R<p\), require at least \(2R+1\) vertices by the distinct roots \(0,\pm1,\ldots,\pm R\); the zero at \(0\) has multiplicity at least two if \(n\) is even. For \(R=2\) the five-vertex bound is attained at \(p=13\).

Here are literal strong modules, specified by outgoing adjacency bitmasks (least significant bit is vertex 0). All data are exact. The first two rows explain the distinction even at a fixed height; the last four rows show two-read persistence on the prime and single-seven lanes, at orders below \(p\) and heights below \(p^2\).

| Outgoing rows | \(H\) | \(Q(2)\) | gcd of all uniform reads |
|---|---:|---:|---:|
| (106,52,121,114,97,64,2) | 49 | 812 | 7 |
| (10,116,121,114,97,65,1) | 49 | 872 | 1 |
| (10,20,25,18,1) | 13 | 104 | 1 |
| (38,124,120,113,97,64,1) | 37 | 740 | 1 |
| (50,116,121,115,96,64,1) | 29 | 580 | 1 |
| (58,76,121,112,98,66,1) | 91 | 1196 | 1 |

The persistent height-49 fixture has no order-seven rotation (its outdegrees differ), so the criterion detects more than symmetry. Its first exterior reads are \(49,861,15729\), versus \(49,921,16929\) for its erasable twin. At height 13 one has \(Q=t^5+8t^3+4t\equiv t(t^2-1)(t^2-4)\bmod13\); \(U_1\) retains 13, while \(U_2\) erases it (1293). This exhibits exactly the phase selectivity the universal strategy lacks.

For the factor-seven fork write \(h=7(40Q(2)+H)\), with offsets \(\pm2H\) at the center \(3h\). For \(p\ge11\), making the center and both exits divisible by \(p\) is equivalent to killing both \(H\) and \(Q(2)\) modulo \(p\). In contrast, a genuine unit fork has \(H=1\); the exits must leave the zero residue of any prime dividing the center. A controller planning a unit move should therefore permit the two prescribed changes of residue, rather than protect the prime on all three states. None of these local tests supplies a new arbitrary-\(p\) inverse or resolves Chan–Pak.

Replays: `python3 uniform-receipt.py` and `python3 bridge-floor-check.py`.
