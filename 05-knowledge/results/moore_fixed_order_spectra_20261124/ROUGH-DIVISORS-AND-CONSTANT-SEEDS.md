# Rough-divisor transfer, the exponent-23 seam, and the order-11 residual

## 1. Sharper distinct-height bound from rough divisors

For either variable-tail guard, write the [spaced quotient](PRODUCT-SPACED-DISTINCT-SEAMS.md) as \(x(au+b)\). With no prescribed height residue, the number of distinct literal \(+2\) pairs through height \(X\) is at least
\[
 \Omega\!\left(\frac{X^{2/3}}{(\log X)^{\eta}(\log\log X)^{3/2}}\right),
 \qquad \eta=1-\frac{1+\log\log 2}{\log 2}.
\]
The two guards separately satisfy this bound, including either scalar ordering. The earlier elementary bound remains valid in each fixed admissible residue class; the improved logarithmic factor below is not claimed with that restriction.

Let \(d=\gcd(a,b),\ a'=a/d,\ b'=b/d\). Fix \(w\ge4\) at least as large as the prime factors of \(a'\). [Ford, Theorem 1(ii) and Remark 1](https://arxiv.org/abs/1901.02548) gives \(\gg N^2/[(\log N)^\eta(\log\log N)^{3/2}]\) distinct \(w\)-rough integers in \((4N^2,8N^2]\) with a divisor \(y\in(N,2N]\). Here \(w\) is fixed and \(N\) tends to infinity. For each integer select one factorization \(n=xy\). Pigeonhole \(y\bmod a'\), retaining a positive constant fraction with the same unit residue \(r\). Choose the least positive \(t\) with \(tr\equiv b'\pmod{a'}\) (if \(a'=1\), take \(t=1\)); this is bounded by \(a'\). Then \(u=(ty-b')/a'\) is an integer, positive for large \(N\), and \(u\asymp N\) uniformly over the finitely many possible \(t\). Both factors \(x,u\) lie in \(O(1)\) dyadic boxes with sides comparable to \(N\); retain the largest box. The known offset and spacing for that box make \(x(au+b)=dt\,n\) the exact height quotient. Distinct integers \(n\) thus give distinct heights. With \(T=O(N)\), the heights are \(O(N^3)\). Take \(N\) a sufficiently small constant times \(X^{1/3}\).

Ford's ordinary multiplication-table bound also gives the matching order for the number of *quotients* in a bounded collection of comparable boxes. This is not an upper bound of the same order for distinct heights: the remainders can distinguish diagrams with the same quotient.

## 2. The exponent-23 target really occurs as a literal seam

For \(E_p=2^{p-1}(2^p-1)\), a 45-vertex tournament and a single reversal give exactly
\[
 E_{23}+1=35{,}184{,}367{,}894{,}529
 \quad\longrightarrow\quad
 E_{23}+3=35{,}184{,}367{,}894{,}531.
\]
Use the six-site A core, masks 177 and 181, with apex out-neighbors 0,1,4 in the core. Prepend the transitive payload `101000000100011011010110010100101110` and append the transitive tail `10`, using the usual common apex contacts. These are \((v,u)=(21511975575,1)\). The composite value \(2^{23}-1=47\cdot178481\) makes this a useful positive control alongside the negative prime-Mersenne result; this witness is in the full two-tail chart and is not a claim about membership in a spaced box.

The replay independently counts the full 45-vertex orders without the height formula or \(2^{45}\) subset DP. Removing the apex gives two lanes, each visiting blocks in ordinal order. For the first lane keep its last contact, and for the second keep its first contact, allowing each lane to be empty. There are nine joint states. A chain uses two choices per site, without enumerating all chain partitions. The six-site core is exhaustive over the split position in each permutation. Concatenation overwrites the first lane's last contact if the next block has any first-lane sites, while the second lane's first contact is unchanged if already occupied. The result is compared with the digit formula. As additional specific checks the two charts give \(E_{19}-39\to E_{19}-37\) at \((2422,44328)\) and \(E_{19}+401\to E_{19}+403\) at \((0,27887181)\).

## 3. Eight constant seeds for the five-cycle order-11 case

For a hypothetical [Gewirtz branch-sign matrix](NO-ORDER13-SIGN-SYMMETRY.md) with five moving 11-cycles and two fixed sites, first gauge away the twist and normalize the first fixed row to all positive on the cycles. Its profile differs from that of the second in \(\delta\) positions. If the fixed-to-fixed sign \(e=+1\), their row distance is \(11\delta+2\) and the route count is \(11\delta\); if \(e=-1\), the distance is \(11\delta\) and the route count is \(11\delta+1\). Membership in \(10,13,\dots,46\) leaves only \((\delta,e)=(2,+1),(3,-1)\).

Set \(A=\left(\begin{smallmatrix}-1&e\\e&-1\end{smallmatrix}\right)\) and let \(P\) be the two profiles as rows. On orbit-constant vectors the operator modulo 3 is
\[
 \begin{pmatrix} A & -P \\ P^{\mathsf T} & S\end{pmatrix},
\]
since \(11=-1\bmod 3\). Thus \(A^2=PP^{\mathsf T},\ PS=-AP,\ S^2=P^{\mathsf T}P\). These equations, plus symmetry of \(S\), suffice for the constant block. Pick one cycle where the two profiles disagree and one where they agree. Choose the symmetric three-by-three complementary principal block in \(3^6\) ways; the linear equation determines the other two rows columnwise, with symmetry determining their entries in the complementary rows first. Checking symmetry and the square equation leaves 12 matrices for each profile, in four classes of sizes 2,2,2,6. The replay prints the eight orbit representatives.

There is a further exact constraint on the missing sign words. At each nontrivial 11th root the five-by-five Fourier matrix is square-zero and hence traceless. The sum of the five diagonal circulant words is therefore a multiple of \(\Phi_{11}\) over \(\mathbf F_3\). Its lag-zero coefficient is \(-5\); hence at every lag either zero or three of the five internal signs are positive. Symmetry identifies opposite lags. Since 3 has order five modulo 11, the nonconstant Fourier problem splits into two degree-five factors exchanged by inversion. Symmetry makes their Fourier matrices transposes of one another under that identification; the remaining square-zero problem can be tested in one \(\mathbf F_{3^5}\) block of rank at most two. Realizing the corresponding \(\pm1\) words and route counts is still open. No order-11 exclusion is claimed here.

Replay: `python3 rough-divisors-and-constant-seeds.py`.
