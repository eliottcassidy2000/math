# Uniform carrier Fourier profile: fast voltage loops and removal of the intermediate range

**Proved, for the central two-tape domain.** The short-loop clock is \(M\log^2(2M)\) for every checksum modulus, not just powers of two. The resulting Fourier profile is uniform on the entire circle: all peaks of order \(T\log T\) occur within polylogarithmic multiples of \(1/T\) of the finitely many rationals with reduced denominator dividing \(2a\), where \(a\) is the carrier's checksum coefficient. This does not control the integral finely enough for a single target.

## 1. Constant phase separation at every modulus

Start from the [cyclic-mixing packet](../chanpak_cyclic_mixing_20261115/README.md). Retain its base address graph modulo \(2M\), voltage \(q(x)=p(x)/2\bmod M\), and root \(o=0\). Its Jensen estimate says that for any padded length-\(L\) walk from \(o\) to \(v\), with voltage \(d\), the squared gap \(|\chi(d) f(v)-(P_\chi^Lf)(o)|^2\) is at most \(B=2(2M)L(\|f\|_2^2-\|P_\chi^Lf\|_2^2)\). One only needs two loops with phases separated by an absolute constant: comparing them bounds \(|f(o)|\) by a constant times \(\sqrt B\). The zero loop then bounds \(|(P_\chi^L f)(o)|\) by the same order, and using a walk to every \(v\) bounds \(\|f\|_2^2=O(B)\).

Here is the short loop construction from the background pass. Write \(M=2^j m\) with \(m\) odd. Take \(r\) minimal with \(2^{r-1}>2M\), and put \(b=2^{r-1}-1\). For reversed \(r\)-bit addresses write
\[
 \Delta_r(c)=\operatorname{rev}_r(c+2M)-\operatorname{rev}_r(c).
\]
If an \(r\)-bit integer \(n\) has \(k<r\) trailing ones, then \(\operatorname{rev}_r(n+1)-\operatorname{rev}_r(n)=3\cdot 2^{r-k-1}-2^r\), by flipping those \(k\) ones to zeros and the next zero to one. The trailing-one counts at \(b,b+2M\) are \(r-1,j+1\), so
\[
 \Delta_r(b+1)-\Delta_r(b)=3(2^{r-j-2}-1).
\]
At width \(r+1\), reselect \(b=2^r-1\). The two displayed differences have gcd 3. Thus for each character order \(s\mid M\), \(s\notin\{1,3\}\), some individual \(\Delta_h(c)\), \(h\in\{r,r+1\}\), \(c\in\{2^{h-1}-1,2^{h-1}\}\), is nonzero modulo \(s\).

Choose such a suffix. The equal-width words for \(2^h(2v)+c\) and \(2^h(2v)+c+2M\), \(0\le v<M\), end at the same address modulo \(2M\). The exact [binary split](../chanpak_true_cut_fourier_20261115/README.md) gives their voltage difference
\[
 q(c+2M)-q(c)-v\Delta_h(c).
\]
The parity boundary vanishes since the upper block \(2v\) is even. As \(v\) varies, the character of the difference averages to zero, so some value has real part at most zero. That pair has phase separation at least \(\sqrt2\). Append the same length-\(\lceil\log_2(2M)\rceil\) reset word to both to return to \(o\), preserving their voltage difference. All lengths are \(O(\log(2M))\). Characters of order 3 instead use a generating loop from the original packet with nontrivial phase, whose distance from 1 is \(\sqrt3\). Pad all words with leading zeros at \(o\). Trivial character contraction on centered functions is as in the original packet.

Consequently the squared centered norm contracts by \(1-c/(ML)\) every \(L=O(\log(2M))\) bits, for an absolute \(c>0\). The total variation error from a point after \(k\) bits is \(O(M\exp[-ck/(M\log^2(2M))])\). The true-cut residue and rational-character error bound in the original packet therefore holds with \(\Lambda(M)=M\log^2(2M)\) for **every** \(M\). For \(M\le(\log T)^\theta\), any fixed \(\theta<1\), its error is \(O_{F,\delta,\theta,A}(T(\log T)^{-A})\) for every fixed \(A\), including after multiplication by any fixed power of \(\log T\). This also gives the original near-target shell corollary on arbitrary moduli in this range. [Replay](verify.js) checks the reversal, drift, and separated-phase identities on 1,125 modulus/order pairs and 160,991 exact drifts.

## 2. Remove the intermediate denominator range without interval completion

Fix a physical form \(F(x,y)=Qxy+f(x)+g(y)+n\), with the documented coefficients, and \(D=\lceil T^\delta\rceil\), \(0<\delta<1/2\). Let \(S_W(\alpha)\) sum \(e(\alpha F)\) over \(x,y\ge D,\ F\le W\), \(0\le W\le T\). Replacing the domain by \(Qxy\le W\) costs \(O_F(T\log(2T)/D)\) uniformly; this changes only the domain, not the phase. Choose \(1\le K\le D^{1/4}\). Partition the \(x\)-rows into consecutive integer intervals starting at \(X\) with full length \(\lceil X/K\rceil\), truncating the last at \(\lfloor W/(QD)\rfloor\). Replace each column-height profile by its height at the right endpoint. For \(x\) in a row interval the height change is \(O(T/(Kx)+1)\), so the discarded area is \(O_F(T\log(2T)/K+T/D)\). Discard the final possibly short row interval at cost \(O_F(T/K+T/D)\), and any rectangles of column height below \(\sqrt D\) at total cost \(O_F(T/\sqrt D)\).

Each remaining rectangle has \(U\gg D/K\) rows, \(V\ge\sqrt D\) columns, and the number of rectangles is \(O(K\log(2T))\). The sum of their \(UV\)'s is \(O_F(T\log(2T))\), so the sum of \(\sqrt{UV}\) is \(O_F(\sqrt{KT}\log(2T))\). If \(|\alpha-u/q|\le1/q^2\) in reduced terms, the approximation to \(Q\alpha\) has reduced denominator between \(q/Q\) and \(q\) and error at most \(Q\) times the inverse square of that reduced denominator. Apply the [weighted bilinear estimate](../chanpak_true_cut_fourier_20261115/README.md) on each rectangle with the original one-variable factors retained. After summing, for \(q\le T\),
\[
 |S_W(\alpha)|/(T\log(2T))
 \ll_F K^{-1}+q^{-1/2}+D^{-1/4}
     +\sqrt{K\log(2T)/D}+\sqrt{Kq\log(2q)/T}.
\]
Unlike Fourier completion of all short intervals, this estimate does not multiply the \(q^{-1/2}\) term by a logarithm. The price for coarse rectangles is the adjustable \(K^{-1}\) term.

## 3. A full uniform profile

For every fixed \(0<\gamma<1/2\), uniformly in \(W\le T\) and \(\alpha\in\mathbb R/\mathbb Z\),
\[
 S_W(\alpha)=
 \sum_{\substack{r\bmod 2a\\ \|\alpha-r/(2a)\|\le 2\log^4 T/T}}
 {\cal A}(r/(2a))\, I_W(T(\alpha-r/(2a)))
 +O_{F,\delta,\gamma}\big(T\log^{1-\gamma}T\big),
\]
for large \(T\). The difference in the integral is taken on the circle with the unique small signed representative, \({\cal A}\) is the arithmetic factor in the sharp-cut packet, and
\[
 I_W(\xi)=Q^{-1}\int_{QD^2}^{W}\log(u/(QD^2))e(\xi u/T)\,du
\]
is zero when \(W<QD^2\). There is at most one summand for large \(T\).

Indeed approximate \(\alpha\) by a reduced \(u/q\) with \(q\le \lfloor T/\log^4 T\rfloor\) and \(|\alpha-u/q|\le 1/(q\lfloor T/\log^4 T\rfloor)\) using Dirichlet approximation. For \(q>(\log T)^{2\gamma}\), the preceding rectangle bound with \(K=\lceil(\log T)^\gamma\rceil\) gives the desired error. For \(q\le(\log T)^{2\gamma}\), apply the uniform rational-character theorem with \(M=q\) and center \(2u/(2q)\); the offset parameter is at most \(2\log^4 T\) and partial summation still yields an error smaller than \(T\) times any fixed negative power of \(\log T\). If the center's reduced denominator does not divide \(2a\), the factor is zero. The claimed expression in the specified neighborhoods follows by applying that same theorem directly at each center there (regardless of which approximant Dirichlet selects). Finally integration by parts on \(I_W\) gives \(O_F(T\log T/(1+|\xi|))\), so already outside radius \((\log T)^\gamma/T\) from each listed rational, the entire sum is \(O(T\log^{1-\gamma}T)\).

This describes the central presentation sum at leading scale on the whole circle. In particular there is no untreated intermediate denominator range for **relative Fourier cancellation**. This is not an estimate of the integrated error on the scale of an individual coefficient, which is the scale needed for Chan–Pak. The triangle inequality with the displayed global error is insufficient; the next question is cancellation of the residual against a prescribed Fourier frequency or a real exact repair in the remaining short interval.
