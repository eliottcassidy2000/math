# The mask-zero carrier's exact parallel-cylinder quotient

Let \(F\) and \(p\) be as in [the local collision theorem](COLLISIONS.md). The previously published cyclic 5760/45 carrier has a different transport system. Here the parameters 100/11 make every persistent **parallel** rectangular alias a transfer of dyadic valuation, with no possible nonzero constant displacement.

First, direct substitution in the recurrence gives
$$
  F(2x,y)=F(x,2y+1)
$$
for all nonnegative \(x,y\). In coordinates \(Y=y+1\) set
$$
 q(x)=\frac{x-[x\text{ odd}]-p(x)}2,\qquad g(Y)=\frac{p(Y-1)}2 .
$$
The same recurrence gives \(q(2x)=q(x)\), \(g(2Y)=g(Y)\). For \(x,Y>0\) with odd parts \(a,b\) and total valuation \(j=v_2(x)+v_2(Y)\) we obtain
$$
F(x,Y-1)=100\cdot 2^j ab+11+22\{g(b)-q(a)\}.
$$
Thus ordered triple \((a,b,j)\) identifies an exact transport orbit. Different triples may nevertheless have equal outputs.

## Classification of all parallel rectangular cylinders

Let \(A,B,C,D\) be powers of two (possibly 1), and \(U,V,W,Z\) suffixes in the corresponding ranges: \(0\le U<A\), \(0\le V<B\), \(0\le W<C\), \(0\le Z<D\). Suppose
$$
F(Ax+U,By+V)-F(Cx+W,Dy+Z)=\text{constant}
$$
for all sufficiently large \(x,y\). Then its constant is zero, and the presentations are related by an integral dyadic valuation transfer.

Indeed, compare coefficients after canceling \(11p(x)+11p(y)\) with the binary split, using even outer variables to suppress any bounded parity terms. The mixed coefficient gives \(AB=CD\). After exchanging sides if necessary, write \(A=kC,D=kB\), \(k\) a power of two. Write \(r_M(s)\) for bit reversal on \(\log_2 M\) bits and \(R_M(s)=M-1-r_M(s)\). The remaining two coefficient equations are
$$
 100C[k(V+1)-(Z+1)]=11[r_A(U)-r_C(W)],\quad
 100B(U-kW)=11[R_D(Z)-R_B(V)].
$$
Decompose \(U=kU_0+u\), \(Z=kZ_0+z\), \(0\le u,z<k\). The reversal identities \(r_A(U)=C r_k(u)+r_C(U_0)\) and \(R_D(Z)=B R_k(z)+R_B(Z_0)\), followed by reduction modulo \(C\) and \(B\), give \(U_0=W,Z_0=V\): these radices are coprime to 11, and the reduced differences have absolute value smaller than the radix. The equations reduce to
$$
 100u=11R_k(z),\qquad 11r_k(u)=100r_k(R_k(z)).
$$
If \(u>0\), the first implies \(R_k(z)=(100/11)u>8u\), so their most significant bit indices differ by at least three. The second implies the opposite: since \(v_2(100)=2\), and the valuation of a nonzero reversed word equals its number of leading zeros, those indices differ by precisely two. Hence \(u=0\), \(R_k(z)=0\), and \(z=k-1\). Now \(Ax+U=k(Cx+W)\), \(Dy+Z=k(By+V+1)-1\), so iteration of the fundamental identity proves the claim, constants and parities included. The proof does not classify identities that exchange the two formal outer variables.

## Exactly how much mass this orbit quotient removes

In the central domain \(x,y\ge D_0=\lceil T^\delta\rceil\), \(T/2<F\le T\), each transport orbit meets the domain in a chain (possibly empty): in \(Y=y+1\) coordinates doubling \(x\) simultaneously halves \(Y\). Move toward smaller \(x\) while \(x\) is even and \(x/2\ge D_0\). The chain's terminal representative is unique. Thus deleting precisely the even-\(x\) presentations with \(x\ge2D_0\) preserves each supported count.

Let \(\ell_T=(1-2\delta)\log T/100\). For every fixed \(r\ge1\), writing the presentation of larger \(x\) as \((2^r u,y)\), the width-\(r\) pair condition is just \(u,y\ge D_0\) and \(T/2<F(2^r u,y)\le T\). Its count is
$$
 (T/2)\ell_T\,2^{-r}+o(T\log T)
$$
by the hyperbola sum: \(F(2^r u,y)=100\cdot2^r u y(1+O(D_0^{-1}))\), and the leading term for \(100\cdot2^r u y\le X\) for \(X\asymp T\) is
\((X/(100\cdot2^r))(\log(X/(100\cdot2^r D_0^2))+O(1))\), summing the floor counts over \(u\). The same calculation with \(r=0\) gives raw mass \((T/2)\ell_T(1+o(1))\). Consequently the canonical chain rule removes asymptotically one half the presentation mass.

The width-\(r\) count is also bounded uniformly by \(O(T\log T\,2^{-r})\), by summing \(T/(100\cdot2^r u)\) over \(u\). The pair has unique width (ratio of its \(x\) addresses), so summing over all widths and both orderings gives total ordered within-orbit off-diagonal energy \(T\ell_T(1+o(1))\). This lower-order contribution cannot supply the order \(T\log^2 T\) second-moment target. Counting each residue modulo 11 separately gives the same asymptotic ratios: the relation \(F(x,y)\equiv x(y+1)\pmod{11}\) uses fixed residue classes of \((u,y)\), and multiplying \(u\) by \(2^r\) just permutes these classes. Every odd residue class modulo 22 is identified by its residue modulo 11.

The quotient keeps one presentation per dyadic transfer orbit, not one per integer output. Any collision between distinct triples \((a,b,j)\) remains.

