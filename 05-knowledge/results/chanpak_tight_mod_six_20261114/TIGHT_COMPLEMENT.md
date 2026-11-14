# A modulus-six selector inside absolute height windows

Let \(g(0)=1\), \(g(2x+1)=g(x)\), \(g(2x)=g(x)+x+(x\bmod2)\) be the known
one-apex transitive-spine recurrence. The tournament on the contact word
\(x\) followed by the reverse complement of \(y\), using their natural
binary widths and with \(x,y>0\), is strong and has
\[
H(x,y)=xy+g(x)+g(y)-1+p(x)p(y),\qquad p(x)=x\bmod2.
\]
This follows by splitting weighted inversions at the word boundary: the
factor \(xy\) counts weights \(2^{j-i-1}\) across blocks and the extra
term is the unit correction at the adjacent positions.

For \(u\ge0\), \(M=2^k\), \(k\ge4\), consider
\[
x=Mu+a,\quad y=Mu+b,\quad a=M/2-1-q,\quad b=M/2+q,\qquad 0\le q<M/2.
\]
All choices except \(u=0,q=M/2-1\) are strong. That boundary has
\(x=0\) and merely gives the transitive count one.
The low words \(a,b\) are complementary at width \(k\).
Their common high prefix is \(u\); the first-order terms in \(xy\)
are independent of \(q\). For a width-\(k\) suffix \(v\), splitting
inversions gives
\[
g(Mu+v)=g(u)+g(v)-1+u Q_v+p(u)(1-v_{\rm first}),
\quad Q_v=M-1-\operatorname{rev}_k(v).
\]
Thus the two suffix coefficients add to \(M-1\), the boundary indicators
add to one, and \(p(x)p(y)=0\). Consequently
\[
H=M^2 u(u+1)-u+2g(u)+p(u)-3+W(a),\qquad W(a)=ab+g(a)+g(b).
\]

Write \(h=q_{\rm first}\) at declared width \(k-1\).
The word \(a\) is the complement of \(q\) at that width, and
\(b\) is a leading 1 followed by \(q\). The reflection identity
\[
g(q)-g(\bar q)=q-\operatorname{rev}_{k-1}(q)+h-p(q)
\]
and the one-bit-prefix concatenation give
\(g(a)+g(b)=2g(q)-q+M/2-2h+p(q)\). Since
\(ab=M^2/4-M/2-q-q^2\), the exact result is
\[
H=B(u,k)+E(q)-2h,\quad E(q)=2g(q)-q^2-2q+p(q),\quad
B(u,k)=M^2u(u+1)+M^2/4-u+2g(u)+p(u)-3.
\]
This cancels the high-prefix terms **over the integers**, not
only modulo a chosen modulus.

For \(q=0,1,2\) the three actual outputs are \(B+2,B,B-2\), since
\(g(0)=g(1)=1,g(2)=3\) and \(h=0\). The center \(B\) is odd. Given any
odd target \(N\in[B-2,B+2]\), its phase modulo six selects
\(q\equiv (B+2-N)/2\pmod3\) in \(\{0,1,2\}\). The same tournament
satisfies \(H\equiv N\pmod6\) and \(|H-N|\le4<6\); hence \(H=N\).
At \(u=0,k=4\) this gives 61, 63 and 65 in one window;
\(u=1,k=4\) gives 573, 575 and 577 in the next example.
There are infinitely many such windows, with two parameters,
but only \(O(\sqrt X+\log X)\) below \(X\): at each \(k\),
\(B\ge M^2u^2+M^2/4-1\), and only \(k=O(\log X)\) contribute.

## Lower windows and small residual phase repair

For every \(1\le j\le k-3\), the choices \(q=2^j-1,2^j\)
occur with \(h=0\). The former has \(g(q)=1\), the latter
\(g(q)=q+1\), so the outputs are \(B+4-4^j,B+2-4^j\).
They cover the two odd residue classes modulo four at their
corresponding two-point height window (error at most 2).

On an aligned four-block \(q=4v+i\), \(i=0,1,2,3\),
the differences \(E(q)-E(4v)\) modulo three are respectively
\[
0,\quad 1,\quad 2-2p(v),\quad 1-2p(v).
\]
This follows directly from the four instances of the \(g\)
recurrence. Every triple of odd residues modulo six occurs,
since \(h\) is constant on an aligned four-block and every output
is odd. If \(D=B-N\ge0\) and the four-block containing
\(\lfloor\sqrt D\rfloor\) is entirely inside the chart, pick a
phase-matching \(q\) there. From \(1\le g(q)\le q+1\) we get
\(H\equiv N\pmod6\) and \(|H-N|\le10\sqrt D+26\). This controls
error by the square root of the residual after choosing the center,
rather than the square root of the original target; it still does
not make that error less than six for arbitrary residuals.

## Exact target inverse for the whole complementary chart

The defect \(q^2+2q-2g(q)-p(q)\) is strictly increasing. For even
\(q\), the next \(g\) weakly decreases; for odd \(q\), the bound
\(g(q+1)-g(q)\le q+1\) gives net increase at least two.
The additional \(2h\) is nondecreasing. Thus \(H\) strictly decreases
in \(q\) at fixed \(u,k\). Also \(W\ge2\), \(W\le M^2/4+2\), and
\(-u+2g(u)+p(u)-3\) is between \(-u-1\) and \(u\).
These give
\[
M^2u^2 < H < M^2(u+1)^2.
\]
The only possible high prefix at prescribed \(N\) and chosen \(M\)
is \(u=\lfloor\sqrt N/M\rfloor\). Binary-search \(q\) and accept
only equality. If \(u=0\) the smallest strong output occurs at
\(a=1,b=M-2\), so \(H\ge M-1\) from the product term. If \(u>0\),
\(M<\sqrt H\). Thus for \(N>1\) it suffices to try powers of two
\(M\le N+1\). This is a polynomial-in-bit-length exact inverse
for this chart, not all tournaments. For \(M\le\sqrt X\)
there are \(O(\sqrt X)\) parameter pairs per width; for larger \(M\)
we must have \(u=0\), and \(ab\ge aM/2\) leaves at most \(2X/M\)
choices. The whole chart has only \(O(\sqrt X\log X)\) outputs up to
\(X\).

## Same-phase approximation at arbitrary odd targets

A separate exact low-bit calculation controls the mod-six phase
uniformly, but not the absolute error to within six. For
\(y=8u+v\), \(p=p(u)\), across \(v=0,\ldots,7\) the values of
\(g(8u+v)-g(u)\) are
\[
(7u+p,3u+p,5u+p+2,u+p,6u+4,2u+2,4u+4,0).
\]
Add \(xv+p(x)p(v)\) and reduce modulo three. The only failing
\(u\bmod6\) state in each \(x\bmod6=0,\ldots,5\) row is
respectively \(5,2,1,4,3,0\). With four low bits, the two three-bit
subfamilies with high prefixes \(2u,2u+1\) cannot both fail.
Every odd residue modulo six occurs among \(H(x,16u+v)\),
\(0\le v<16\).

Given odd \(N\ge289\), set \(x=\lceil\sqrt N\rceil\) and
\(u=\lfloor N/(16x)\rfloor\ge1\). Choose a phase-matching \(v\).
Then \(|xy-N|<16x\), \(y\le x+15\), and \(1\le g(z)\le z+1\)
gives \(|H-N|<18x+17\). Every such target therefore has a
same-phase strong witness at sublinear absolute distance.
Bridging this coarse approximation to the exact windows remains open.

Run **python3 replay_balanced_windows.py** for the complete low-bit phase table, the
complementary formula and modulus-four/modulus-six windows, and
literal subset-path DP on small cases.

Run `python3 replay_balanced_windows.py TARGET` to select a witness from
this complementary chart and print literal outgoing masks, or reject
just this chart. The modulus-six height window is an analytic guarantee
for the stated targets; the program also checks other values by its
exact inverse.
