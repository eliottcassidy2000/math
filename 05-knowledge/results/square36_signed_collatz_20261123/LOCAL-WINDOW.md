# The sharp window for short square-sum detours

**PROVED.** Every positive even integer \(M\) has a simple path from \(M-1\) to \(M+1\) on positive odd labels, with at most five square-sum edges, in the window \(M\pm O(M^{3/4})\). For paths of bounded length the exponent \(3/4\) is optimal. At a Euclid midpoint \(E_p=2^{p-1}(2^p-1)\) of odd exponent \(p\ge5\), three edges already give this order, and one path uniquely attains the smallest possible window for three edges. These are integer-label statements. [Actual tournament-height intermediates](DETOUR.md) and [the marked-arc lift](LIFT.md) require further information.

## An exact remainder in five reflections

An odd-label square-sum step has the form \(x\mapsto4a^2-x\). It negates the residue modulo four, so a path between the two neighbors of \(M\) has odd length. After an odd number of steps it reaches \(M+1\) exactly when the alternating sum of the squared half-roots is \(m=M/2\). This is the reflection mechanism of [THM-4510, sum graphs and Pythagorean zigzags](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md).

Choose the smallest \(n\ge\sqrt m\) of the same parity as \(m\). Put \(D=(n^2-m)/2\), \(b=\max(1,\lfloor\sqrt D\rfloor)\), and write \(D=ab+c\), \(0\le c<b\). Use the five successive half-roots
\[
      n+a+1,\quad n+a+b+1,\quad n+b+1,\quad n+c+1,\quad n+c.
\]
Expansion gives the alternating sum of their squares as \(n^2-2(ab+c)=m\). Put \(x_0=M-1\) and reflect in order. Positivity holds even at the smallest input. The first and fourth interior labels satisfy
\[
 x_1=4(n+a+1)^2-M+1>0,\qquad
 x_4=4(n+c)^2-M-1\ge M-1.
\]
The second and third are \(x_2=M-1+4b(2n+2a+b+2)\) and \(x_3=M+1+4(2n+2c+1)\). All labels are odd. Removing each closed subwalk between repeated labels leaves a simple path with the same endpoints and no larger window. At \(M=2\) the resulting path is \(1,3\); at \(M=4\) it is \(3,61,39,25,11,5\).

Since \(n-\sqrt m<2\), the deficit is \(O(\sqrt m)\) and \(a,b,c=O(m^{1/4})\). For each successive half-root \(k\), one has \(4k^2=2M+O(M^{3/4})\). Starting at \(x_0=M-1\), the five relations \(x_i-M=(4k_i^2-2M)-(x_{i-1}-M)\) give the claimed common window.

There is also a three-edge construction for every \(M\) with a larger uniform window. Choose \(a\equiv m+1\pmod2\) smallest above \(\sqrt m+2\), put \(b=(a^2+1-m)/2\), and use half-roots \(a,b,b-1\). Their alternating square sum is \(m\), with \(a,b-1>\sqrt m\), so both interior labels are positive. All three roots are \(O(\sqrt M)\).

## Why no bounded-length construction improves the exponent

Write \(W=\max_i|x_i-M|\) for a path of odd length \(\ell\). Set \(m=n^2+n\), \(R=\sqrt m\). Any half-root obeys
\[
 |k_i-R|\le\frac{W}{2R},\qquad
 |k_i-n|\le\frac{W}{2R}+\tfrac12 .
\]
The first bound follows by factoring \(k_i^2-R^2=(x_{i-1}+x_i-2M)/4\). Write \(t_i=k_i-n\). Reducing the alternating square sum modulo \(2n\) gives \(\sum_i(-1)^{i-1}t_i^2\equiv n\pmod{2n}\). Its absolute value is at least \(n\), hence
\[
             W\ge 2R\left(\sqrt{n/\ell}-\tfrac12\right).
\]
In particular any fixed length bound costs \(\Omega(M^{3/4})\) along this infinite sequence.

## The unique narrowest three-edge path at a Euclid midpoint

Let \(n=2^{p-1}=r^2\), with \(p\ge5\) odd. Put \(M=2n^2-n\), so \(m=n^2-n/2\). For each sign, half-roots
\[
            n\pm r/2,\qquad n\pm r,\qquad n\pm r/2
\]
give a path. For the negative sign the interior labels are \(C+1,C-1\), where
\[
    C=2(n-r)^2=M-4nr+3n,\qquad W_* =4nr-3n+1.
\]
Both interior labels are positive. There is no direct edge because the odd part of \(2M\) is 3 modulo 4.

To prove uniqueness, first write the three half-roots as \(n+t_i\) with \(\max|t_i|\le r\). Set \(s=t_1-t_2+t_3\), \(Q=t_1^2-t_2^2+t_3^2\). Then \(2ns+Q=-n/2\), \(-n\le Q\le2n\), and \(Q\equiv s\pmod2\). Since \(n\) is divisible by four, the equation makes \(s\) even, and these inequalities force \(s=0\). It follows that \(t_2=t_1+t_3\) and \(t_1t_3=n/4\). Hence \(|t_2|\ge r\), with equality precisely at \(t_1=t_3=\pm r/2\). The two displayed paths are the only minimizers of the maximum half-root displacement.

Now take any path of length three in the proposed vertex window \(W\le W_*\). Every half-root satisfies
\[
             (n-r)^2-\tfrac12\le k_i^2\le m+W_*/2<(n+r)^2 .
\]
As the roots are positive integers, the preceding classification applies and rules out the positive sign. Thus the negative path is the unique optimum in the actual vertex window. For unrestricted odd length \(\ell\) the preceding congruence argument works with residue \(-n/2\) instead, giving \(W\ge2\sqrt m(\sqrt{n/(2\ell)}-1/2)\). Consequently even allowing more fixed steps preserves the sharp exponent along the Euclid sequence itself.

At \(p=61\), \(W_*=4,951,760,153,682,756,585,775,955,969\). Run `python3 -B detour_window.py` from this directory for the exact path, the small cases, and an exhaustive minimizer check at exponents 5, 7, 9 and 11. Primality of \(2^p-1\) plays no role in the proof.
