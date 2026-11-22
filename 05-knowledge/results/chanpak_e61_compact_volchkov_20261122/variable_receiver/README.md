# Exact single-apex receiver inversion and a Mersenne padding barrier

**Proved.** The familiar exponent-31 core with a variable one-apex receiver admits a four-state read and an exact modular suffix decoder. At any exponent \(p\ge11\) for which \(2^p-1\) is prime, an exact perfect seam in this chart must have at most \(p-2\) trailing receiver sinks. No exponent-61 hit was found in the explicitly bounded scan below.

## Four-state receiver

Let \(V(v)\) be a transitive spine with binary word \(v\), with bit 1 meaning its apex beats the corresponding spine vertex. Its apex is the sole receiver tag for the next extension. At \(v=0\), use just this apex. Let \(Q_v\) count ordered two-path covers allowing empty slots, and \(m_v\) count those with designated first path ending at the receiver apex. Then \(h_v=H(V(v))=1+p(v)\) for the binary recurrence \(p(2s+1)=p(s)\), \(p(2s)=p(s)+s+[s\text{ odd}]\) for \(s>0\), with \(p(0)=0\).

Fix the slot containing the apex, color each spine vertex \(B,A,O\) (before it, after it, other slot), and put each color subsequence in spine order. The only requirements are that the first \(A\), if present, is bit 1 and the last \(B\), if present, is bit 0. The two flags “\(A\) has started” and “last \(B\) is legal” give a four-state counter: \(O\) leaves both fixed, \(A\) starts only on 1, and \(B\) resets legality. Sum the two final legal states and double. For \(m_v\), force \(A\) absent and fix the apex slot. The last \(B\) position, or lack of one, gives
\[
m_v=1+\sum_{i:\text{bit }i=0}2^i,
\]
where \(i\) is indexed from the left starting at 0. The binary recursion's contribution \(s+[s\text{ odd}]\) of each zero at position \(i\) is at most \(2^i\), so \(h_v\le m_v\).

Write \(v=2^t(w+1)-1\) by stripping trailing ones, so \(w=0\) or \(w>0\) is even. Those \(t\) vertices are ordinal sinks in \(V(v)\); consequently \(p(v)=p(w)\), \(m_v=m_w\), and \(Q_v=2^tQ_w\).

## Physical formula and exact suffix inverse

For the [published exponent-31 gate](../chanpak_e31_physical_20261122/README.md), use source rows \((4,5,24,19,35,15)\), tag 1, and mark \((1,3)\). Put this source before \(V(v)\) in an ordinal join, then attach the outer marked apex and reversed-complement right word \(y\). The balanced-core transfer gives
\[
L=A2^t(y+1)+F(1+p(y))-D-25[t=0,\ y\text{ even}],
\quad (A,F,D)=(25(w+1)+28Q_w,\;25h_w,\;10+28m_w).
\]
The marked reversal gives exactly \(L+2\). Small physical DP independently checks the formulas and literal edge and recovers the known exponent-31 witness \((w,t,y)=(6,0,1973896132951836)\). That certificate is not new.

For target \(N\), put \(B=A2^t\), \(u=y+1\), and \(V=N+D\), adjusted by \(25\) when \(t=0\) and \(y\) is even. The equation \(Bu+F(1+p(u-1))=V\) implies
\[
 \left\lceil V/(B+F)\right\rceil\le u\le\left\lfloor(V-F)/B\right\rfloor,
 \qquad Bu\equiv V\pmod F .
\]
This interval follows from \(0\le p(u-1)\le u-1\). Enumerate the intersection with the arithmetic progression and the prescribed parity, if any, and test the full carry read. Its step is \(F/\gcd(A,F)\) since \(F\) is odd. For fixed \(w\), the residue halves modulo that odd step as \(t\) increases. This is complete for each chart with finite candidate count, without case guesses about the binary suffix.

## Mersenne obstruction

Let \(E_p=2^{p-1}(2^p-1)\), with \(2^p-1\) prime. The constants satisfy \(A\) odd, \(A>1\), \(F<D-1<A\), and \(F/A\le25/81<1/2\). Indeed \(h_w\le m_w\), \(h_w\le w+1\), \(Q_w\ge2(w+1)\) by taking a Hamilton path on each side of every partition of the \(1+\operatorname{bitlength}(w)\) receiver vertices, and \(Q_w\ge m_w\).

For \(1\le t\le p-1\), if \(D-1\) and \(F(y+1)\) are both less than \(2^t\), the required congruence
\[
F(1+p(y))\equiv D-1\pmod{2^t}
\]
at an exact \(L=E_p-1\) would force equality. Then \(A2^t(y+1)=E_p\), and primality forces \(A=2^p-1\), so \(y=2^{p-1-t}-1\). This suffix is all ones or empty; \(p(y)=0\), contradicting \(F<D-1\). Thus a necessary condition in this range is
\[
2^t\le D-1\quad\text{or}\quad 2^{2t}\le (F/A)(E_p+D-1).
\]

If \(D-1<2^{p-1}\), no \(t\ge p-1\) works. At \(t=p-1\), use \(F(y+1)\le (F/A)(E_p+D-1)2^{-t}<2^{p-1}\) and apply the contradiction above. For \(t\ge p\), \(L\ge A(2^t-1)\) rules out \(t\ge2p-1\). In the remaining range the required residue of the carry is \(2^t-2^{p-1}+D-1\ge2^{p-1}\), whereas the same bound makes the actual read smaller.

For \(p\ge11\), the other case \(D-1\ge2^{p-1}\) implies active width \(b\ge p-5\) since \(m_w\le2^b\). In a nonempty active tape the first bit is 1 and the last is 0. Forcing them respectively to \(A,B\) and choosing any color at each inner vertex shows \(Q_w\ge2\cdot3^{b-2}\). Hence \(A\ge56\cdot3^{p-7}\ge2^{p+1}\). But \(L\ge A(2^t-1)\) forces \(A<2^{p+1}\) at \(t\ge p-1\). Thus the same exclusion holds and any exact solution satisfies \(t\le p-2\).

## Exponent-61 scan

The C++ replay considers every active even word of width at most 19 (including the empty active word) and every remaining active even word of width at most 60 with Hamming weight at most 3. For each prefix it considers every nonnegative padding width before the obvious height cutoff, with \(t=0\) parities separate. Every chart with at most 64 arithmetic candidates is tested in full; the larger progressions are **not** exhausted.

- 295,436 active prefixes;
- 12,241,836 empty charts or failed divisibility;
- 693,842 nonempty charts scanned;
- 4,761,061 larger charts omitted, all at \(0\le t\le51\);
- 8,630,343 full recurrence tests; **no exact hit**.

The best error among the arithmetic-grid candidates actually tested was \(-840141300527699250\), at \((w,t,y)=(22,45,11408404609076599471)\). It is not an optimum over all seams or even over every address in the variable receiver chart. Replays: run python check.py; compile scan.cpp with C++17 and run the executable with arguments 61 19 64. The exact recorded output is in scan.out.
