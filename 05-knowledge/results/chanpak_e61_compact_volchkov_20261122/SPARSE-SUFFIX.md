# Two exact sparse-language inverses for the E61 Fibonacci receiver

The previous receiver selector specifies the count but does not solve the unrestricted suffix at \(E_{61}\). These tests include all 610 normalized receivers through order 9, not just the three at total order 123. The six-vertex marked core and the receiver normalization are unchanged. At \(N=E_{61}-1=2^{121}-2^{60}-1\) the equation is
\[
 N=K_{y\bmod2}+cy+h f(y),
 \qquad f(2z+1)=f(z),\quad f(2z)=f(z)+z+[z\text{ odd}],
\]
with \(f(0)=0,\ 0\le f(y)\le y,\ c,h>0\), and \(h\) odd.

## 1. Eliminate the leading run exactly

Appending \(g\) ones to a positive address \(x\) sends \(x\) to \(2^g(x+1)-1\) and leaves \(f(x)\) unchanged. Appending \(g\) zeros sends \(x\) to \(2^gx\) and increases \(f(x)\) by \((2^g-1)x+[x\text{ odd}]\); only the first zero can use the parity indicator. In a run decomposition starting with \(1^a\), let \(A=2^a-1\). Once every later run length is fixed, both \(y\) and \(f(y)\) are affine forms in \(A\). This yields a single exact division to obtain \(A\), followed by the tests \(A>0\) and \(A+1\) a power of two. The bit-length budget is \(L=\operatorname{bitlength}(\lfloor(N-\min K_i)/c\rfloor)\). Consequently addresses with a fixed run cap \(r\) are completely inverted by \(O(L^{r-1})\) suffix tuples, and prescribed total orders can be imposed without changing the recurrence.

For a four- or five-run word \(1^a0^b1^j0^d1^e\) with \(a,b,j,d\ge1,e\ge0\), put \(X=2^{b+j}A+2^j-1\). Then
\[
 y=2^{d+e}X+2^e-1,\qquad
 f(y)=(2^b-1)A+(2^d-1)X+2.
\]
The slope of the equation for \(A\) is strictly positive. Enumerate all \(b+j+d+e<L\) with those lower bounds, solve the affine equation, and replay on the original bits. The C++ replay checked **111,999,172** four-run tuples and **2,949,122,622** tuples for four or five runs, with no hits. Together with the inherited exhaustion of at most three runs, no address with at most five runs at this target works in any of the 610 charts.

## 2. Wildcard suffixes with a bounded number of defects

Let \(B=2^\ell\), \(y=Bz+r\) with \(0\le r<B\), and let \(r^*\) reverse the \(\ell\) bits of \(r\), including leading zeros. The exact binary split gives
\[
 f(y)=f(z)+(B-1-r^*)z+f(r)+[z\text{ odd}][r<B/2].
\]
For each prescribed bit suffix \(r\) and high parity \(i=z\bmod2\), write
\[
 C=cB+h(B-1-r^*),\quad
 M=N-K_{r\bmod2}-cr-h(f(r)+i[r<B/2]).
\]
The exact remaining equation is \(Cz+hf(z)=M\). A negative \(M\) is impossible; otherwise the only candidates are
\[
 \max(1,\lceil M/(C+h)\rceil)\le z\le\lfloor M/C\rfloor,\qquad
 Cz\equiv M\pmod h,\qquad z\equiv i\pmod2.
\]
The lower bound 1 tests actual words of width greater than \(\ell\), not left-padded short words. Here it is automatic: \((c+h)(2^\ell-1)+\max K_i<N\) for every tested chart and \(\ell\le60\). Solve the congruence by dividing by \(\gcd(C,h)\); the resulting modulus is odd, so the parity and residue combine into one progression. Replay \(f(z)\) only on that list. This is necessary and sufficient and works for arbitrary explicitly listed suffixes, not only homogeneous ones. In a Hamming ball of radius \(s\) about either constant word, the suffix-list size is \(2\sum_{j\le s}\binom{\ell}{j}\). The interval width is at most \(1+hM/(C(C+h))\); the residue often eliminates most of it.

For \(\ell=54,s=3\), the exact certificate has **32,073,800** chart/suffix pairs and only **53,049** surviving prefix checks, with no hit. This extends the homogeneous-tail test: the lowest 54 bits may contain three arbitrary exceptions, including interior defects that split a long terminal run into several shorter ones. It excludes this entire finite union of suffix cylinders, not every binary word.

## Replay

Generate the chart stream in the directory of the previous packet's check.py:

    from check import strong_readers, graft_chart
    for v in strong_readers(9):
        print(len(v), *graft_chart(v), ",".join(map(str, v)))

Use one line (order, c, h, K0, K1, receiver) per chart and run the two C++ replays with that stream on stdin. The full integer target fits the signed 128-bit types used here, and all arithmetic in the tests is exact. None of the exclusions applies to a different marked core, receivers beyond the bank, or addresses in the complementary suffix language.

## 3. Fixed-order two-run elimination by a quadratic

If the suffix width \(L\) is fixed (equivalently the host order is fixed once the receiver is fixed), leave the **first and last** run lengths unknown and fix all \(r-2\) middle lengths. Let \(X,Y\) be the powers of two corresponding to the first and last lengths. They satisfy \(XY=2^M\), where \(M\) is the remaining width. Before the last run, write \(y=AX+B\), \(f(y)=CX+D\), starting the affine recursion with \((A,B,C,D)=(1,-1,0,0)\). The two possible final counts are then exactly
\[
\begin{array}{ll}
\text{last bit }0:& (c+h)A\,XY+h(C-A)X+(c+h)B Y+K_0+h(D-B+1),\\
\text{last bit }1:& cA\,XY+hCX+c(B+1)Y+K_1-c+hD.
\end{array}
\]
The sign in the zero case uses that the preceding run is of ones, including the case of only two runs. Writing the applicable expression as \(aXY+bX+gY+d\), the target reduces to
\[
 bX^2-(N-d-a2^M)X+g2^M=0.
\]
A root must be a power of two \(X\ge2\) dividing \(2^M\), with \(Y=2^M/X\ge2\). These conditions and an exact replay are sufficient. This requires only \(O(L^{r-2})\) middle vectors for fixed \(r\), not \(O(L^{r-1})\). With fewer than two runs, replay the sole constant-one word directly.

For an odd \(r\ge3\), this has **at most one** possible first run for each middle vector: the zero middle run makes \(C>0\) and \(B\le-2\), preserved on all extensions, so \(b=hC>0\), \(g=c(B+1)<0\), and the target read as \(bX+g2^M/X\) is strictly increasing. For an even \(r\ge2\), \(C-A<0\) and \(B<0\): appending zeros preserves \(C-A\), and appending ones strictly decreases it; both updates keep \(B\le-1\). Thus \(b,g<0\), and there are at most two solutions. One may solve via the integer discriminant or split the dyadic powers at the unique turning point by the signs of \(b2^j-g2^{M-j-1}\) and search each monotone branch. The fixed-width reduction does not depend on any conjectural properties of prime counts or valuations.

Additional wildcard replays, all with no hits:

| Low-bit width | Maximum defects | Chart/suffix pairs | Exact prefix tests |
|---:|---:|---:|---:|
| 48 | 2 | 1,435,940 | 9,746,835 |
| 50 | 3 | 25,468,720 | 10,779,819 |
| 52 | 3 | 28,644,380 | 757,927 |
| 54 | 3 | 32,073,800 | 53,049 |
| 56 | 3 | 35,766,740 | 3,676 |
| 58 | 3 | 39,732,960 | 265 |
| 60 | 3 | 43,982,220 | 24 |

The quadratic replay tested **346,491,123** middle vectors for six and seven runs at total order 123 (three eligible charts), and **147,925,910** for six runs at total order 122 (**31** eligible charts across the whole 610-receiver bank). Both lists are empty. At order 123 an independent direct run sweep also checked **3,013,417,231** six-/seven-run words after the necessary leading/trailing bit filters, with the same result. Thus any hit in this construction has at least eight runs at order 123 and at least seven runs at order 122. The latter is not merely a test on the top-order three receivers: the inherited receiver bound at total order 122 is nine, so the bank covers every normalized receiver there.
