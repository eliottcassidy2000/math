# Growing-modulus Chan–Pak factors by short voltage loops

**Proved.** For the physical two-tape carriers in the [sharp-cut Fourier packet](../chanpak_true_cut_fourier_20261115/README.md), the arithmetic factors hold uniformly up to \(M\le(\log T)^\theta\) for every fixed \(\theta<1/3\). For powers of two the range includes every \(\theta<1\). The background estimate used an \(M^6\) mixing clock (giving \(\theta<1/6\)); this argument uses checksum translation symmetry to replace it by \(M^3\log^2(2M)\), or \(M\log^2(2M)\) for powers of two. This is not an exact-target theorem.

## Short loops on the address graph

Write \(q=p(x)/2\). The fair-bit transition is \((x,q)\mapsto(2x+b,q+d_b(x))\), modulo \(2M,M\), where \(d_1=0\) and \(d_0=(x+(x\bmod2))/2\). Uniform joint stationarity follows by counting two incoming choices at each destination. Put \(S=2M\) and \(l=\lceil\log_2 S\rceil\). In the first coordinate any \(u\) can be joined to any \(v\) in exactly \(l\) steps by appending the integer \(v-2^l u\bmod S\) in \(l\) bits. The voltage of a walk is its added checksum, independent of starting checksum.

Use \(o=0\) and choose length-\(l\) walks \(A_v:o\to v\), \(B_v:v\to o\), of voltages \(\alpha_v,\beta_v\). The loops \(A_vB_v\) and \(A_v e B_w\) for every base edge \(e:v\to w\) generate the full voltage group. Their voltage differences include \(d(e)+\alpha_v-\alpha_w\), by subtracting \(A_wB_w\). Summing these differences around a closed walk gives its voltage. The [fixed-length suffix programmer](../chanpak_absolute_steering_20261115/PREFIX_PROGRAMMING.md) ensures every voltage occurs on a return from \(o\) to \(o\). Thus every nontrivial checksum character separates one of these short loops from the zero loop. Pad at \(o\) with zeros to common length \(L=2l+1\).

## Fourier contraction

On a checksum character \(\chi\), the transition on address functions is \(P_\chi f(v)=\tfrac12\sum_{e:v\to w}\chi(d(e))f(w)\). In the uniform norm on the base let \(f_j=P_\chi^j f\), \(\epsilon_j=\|f_j\|_2^2-\|f_{j+1}\|_2^2\). Stationarity gives the exact variance identity
\[
 \epsilon_j=\frac1{2S}\sum_{e:v\to w}
 |\chi(d(e))f_j(w)-f_{j+1}(v)|^2 .
\]
Along any length-\(L\) walk \(o\to v\) of voltage \(d\), telescope in reverse time order from \(f_L(o)\) to \(f(v)\), multiplying each difference by the character of the preceding edges. Each time layer is used just once. Cauchy–Schwarz gives
\[
 |\chi(d)f(v)-f_L(o)|^2\le 2SL\sum_{j<L}\epsilon_j=:B .
\]
Use both the zero loop and the separated padded loop of voltage \(d_*\). For \(g=|1-\chi(d_*)|\), this gives \(g|f(o)|\le2\sqrt B\) and \(|f_L(o)|\le(1+2/g)\sqrt B\). Using the padded walks to other \(v\) gives \(|f(v)|\le(2+2/g)\sqrt B\); since \(g\le2\),
\[
 \|P_\chi^L f\|_2^2 \le (1-g^2/(72SL))\|f\|_2^2.
\]
For nontrivial characters \(g\ge2\sin(\pi/M)\ge4/M\). For the trivial character on the mean-zero space, averaging \(|f(v)-f_L(o)|^2\le B\) gives \(\|f\|_2^2\le B\). Iterate and take the maximum norm over the orthogonal Fourier summands. For the full transition and its adjoint the centered operator norm at time \(r\) is at most \(C\exp(-cr/\Lambda(M))\), where \(\Lambda(M)=M^3\log^2(2M)\), for absolute positive \(C,c\). Starting from a point, total variation is at most \(CM\exp(-cr/\Lambda(M))\), since the centered point density has \(L^2\) norm at most \(\sqrt{2}M\).

**Dyadic phase kick.** For \(M=2^a\) use \(o=-1\bmod2M\) and its zero-voltage one loop. The word \(0\,1^{j+1}\,0\,1^{a+1}\), \(0\le j<a\), returns to \(o\) and adds \(-2^j\bmod M\). The two zero positions are above the final \(a+1\) bits; the first zero contributes zero, and before the second the address is \(-1-2^{j+1}\bmod2M\). A character of exact order \(2^s\) sees phase \(-1\) by choosing \(j=s-1\). Use common length \(2a+3\) with padding at \(o\), and set \(g=2\) above. This replaces \(\Lambda(M)\) by \(M\log^2(2M)\). For \(M=1\) only the trivial character occurs.

## True-cut arithmetic

For \(D=\lceil T^\delta\rceil\), \(0<\delta<1/2\), the error in equation (4) of the sharp-cut packet can be taken as
\[
 {\cal E}=O_{F,\delta}\!\left(
 T\log(2T)M\exp[-c_\delta\log T/\Lambda(M)]
 +T^{1-\delta/2}\log(2T)\right).
\]
The bound is uniform in \(M\), the residue, and the sublevel \(0\le W\le T\). In fact this holds for any bounded-by-one test of the pair of modular states, compared to two independent uniform states on the hyperbolic approximation. Decompose any integer interval into aligned dyadic blocks, discarding those shorter than \(\sqrt D\) at cost \(O(\sqrt D)\) points. In the remaining blocks the initial prefix is deterministic, followed by at least \(\frac12\log_2D\) fair bits. Thus an interval of length \(R\) has discrepancy \(O(R\varepsilon+\sqrt D)\), with \(\varepsilon=CM\exp[-c\log D/\Lambda(M)]\). Mix the \(y\)-interval in each row of \(Qxy\le W\): the sum of row lengths is \(O(T\log(2T))\), and there are \(O(T/D)\) rows. Now apply summation by parts on each dyadic \(x\)-segment of scale \(U\) to the conditional mean over \(y\), weighted by the nonincreasing row length. The cost per segment is \(O((T/U)(U\varepsilon+\sqrt D))\). Summation gives \(O(T\log(2T)\varepsilon+T/\sqrt D)\). The hyperbola-to-true-cut cost stays \(O_F(T\log(2T)/D)\). A rational character is one bounded test, so no modulus-sized loss is incurred when deriving the corresponding near-center major arc by partial summation, with error \(O((1+|\xi|){\cal E})\) at \(r/(2M)+\xi/T\).

For the order-five mask-zero carrier, the residue probability satisfies \(\sigma_M(b)\ge10/(11M)\) for every odd \(b\). In the stated polylogarithmic ranges the error is smaller than \(T\) times any fixed negative power of \(\log T\). For example, each such residue has asymptotically \(\sigma_M(b)(1-2\delta)T\log^{1-A}T/Q\) physical presentations in the true shell \(T(1-\log^{-A}T)<F\le T\), for every fixed \(A>0\), by differencing the two sublevels and applying the hyperbolic area derivative. This gives simultaneous absolute-window and growing-congruence control but the window is much wider than the modulus. Neither exact equality at a prescribed integer nor positive density follows.

The background argument also follows from the same telescope with no Fourier split: pad full-state paths to length \(2M^2-1\) at one self-loop and use the variance in the full space. [Small exact replay](verify.py) checks 197,376 generating-loop words through modulus 256 and all 465 dyadic phase kicks through exponent 30. The mathematical proofs are not restricted to those ranges.
