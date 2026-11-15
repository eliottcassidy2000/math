# Two-ended twisted collisions and exact inverse carrier search

For the digital checksum use
$$
p(0)=0,\qquad p(2u+1)=p(u),\qquad p(2u)=p(u)+u+(u\bmod2),
$$
and let $G(u)=cu+kp(u)+\phi(u\bmod2)$, with fixed integers $c,k>0$. An additive constant does not affect collisions. The [row hypothesis](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/hypotheses/HYP-9752-row-collision-lemma-a-carrier-row-cu-plus-kp-u-has-linear-energy-which-would-give-the-chan-pak-spectrum-positive-lower-density.md) remains open. The result below improves control of the sector where all differing bits are confined to one short interior span. Separately, there is now a deterministic square-root-scale inverse algorithm for the two-tape carrier, at exact integer resolution.

## 1. Peel both ends of a collision

Let $u>u'$ be distinct $L$-bit inputs, allowing leading zeros. Write
$$
u=2^j(2^s w+x)+t,\qquad
u'=2^j(2^s w+y)+t,\qquad
x\ge 2^{s-1}>y,\qquad 0\le t<2^j.
$$
Here $j$ is the lowest differing bit position, and $j+s-1$ is the highest. Thus $a=x-y>0$ and $d=\operatorname{rev}_s(x)-\operatorname{rev}_s(y)$ are both odd. Use $M=2^s$, $W=2^{L-s-j}$, $B=2^j$, and $N=2^L$, so $WB=N/M$. For $j\ge1$ put $r=\operatorname{rev}_j(t)$ and $\varepsilon=(x\bmod2)-(y\bmod2)$.

The inherited binary split identity is, for $V=2^m\ge2$ and $0\le b<V$,
$$
p(Vz+b)=p(z)+p(b)+(V-1-\operatorname{rev}_m(b))z+[b<V/2](z\bmod2).
$$
Apply it first at the common suffix $t$ and then at the high prefix $w$. In the checksum difference this gives
$$
p(u)-p(u')=p(x)-p(y)-dw-[w\text{ odd}]
 +(2^j-1-r)a+\varepsilon[r\text{ even}].
$$
Here the last indicator is correct also at $j=1$: the MSB of $t$ is the parity of its reverse. The explicit parity readout in $G$ cancels. Hence the *exact* collision equation is
$$
d w+a r+[w\text{ odd}]-\varepsilon[r\text{ even}]
 =\big((c/k)2^j+2^j-1\big)a+p(x)-p(y).          \tag{1}
$$
On each of four parity choices it is a lattice line. If $g=\gcd(a,|d|)$, consecutive integral solutions move by $a/g$ in $w$ and by $|d|/g$ in $r$. Consequently each core gives at most
$$
O\left(1+\min(Wg/a,Bg/|d|)\right)
$$
candidates. We need not assume the intercept is an integer, nor assume generic coprimality. For $j=0$ the original one-ended split leaves at most one $w$ per parity for each pair $x,y$, so $O(M^2)$ suffices for that entire sector.

## 2. Moment of the primitive directions

The pair $(a,d)$ determines the signed binary difference vector $(x_i-y_i)_{i=0}^{s-1}$. Indeed, subtraction of two possible vectors would give a nonzero polynomial with coefficients of absolute value at most two vanishing at $2$ and $1/2$. Remove the largest common power of the indeterminate. The remaining constant term is even by evaluation at $2$ and nonzero, hence has absolute value two. The sum of absolute values of the other terms evaluated at $1/2$ is strictly smaller than two, a contradiction. This is the signature rigidity from the [two-tape displacement packet](../chanpak_dyadic_collision_ports_20261115/COLLISIONS.md), applied here to the two-ended row template. Each zero in the vector allows at most two choices of the shared digit; the other digits are fixed. Thus for $P>1$, with $Q=P/(P-1)$ and $\gamma=\log_2(2+2^P)/P$, the multiplicity function on $(a,d)$ has $\ell^P$ norm at most $M^\gamma$.

Write positive absolute directions $(a,|d|)=g(A,D)$ with $A,D$ coprime. The number of allowed $g$ is at most $M/\max(A,D)$. The normalized geometric-mean factor has bounded moment in the primitive slopes:
$$
\sum_{A,D\ge1}\frac{(AD)^{-Q/2}}{\max(A,D)}<\infty.
$$
For example in $A\ge D$, the inner sum of $D^{-Q/2}$ is $O(A^{1-Q/2})$ for $Q<2$, $O(\log(2A))$ at two, and $O(1)$ above two; every resulting outer sum converges. Thus the $Q$-th moment of $\sqrt{WB}g/\sqrt{a|d|}$ over all $a,|d|\le M$ is $O_Q(M(WB)^{Q/2})$. For an unbalanced rectangle we also have $O_Q(M\min(W,B)^Q\log(2M))$ for the actual minimum factor: when $B\le W$, use $B/D$, sum first $1/\max(A,D)$ over $A\le M$, and then sum $D^{-Q}$. Signs of $d$ cost just a factor two.

Hölder on the multiplicity and direction factors proves that the number at each specified $(s,j)$ is at most
$$
O_P\left(
 M^2+
 M^{1+\gamma-1/P}\min\left(\sqrt{WB},\min(W,B)\log(2M)^{1/Q}\right)
\right).                                                     \tag{2}
$$
In particular this is $O_P(M^2+\sqrt N M^{1/2+\gamma-1/P})$. With $P=14/5$, the reciprocal of twice the exponent of $M$ in the second term is $0.392794163029\ldots$. For any fixed $\epsilon>0$, all pairs with differing span at most $(0.392794163029-\epsilon)L$ contribute $o(N)$ in total (at most $L^2$ sectors, with a fixed power margin). This goes beyond the previous $N^{1/3-o(1)}$ reverse-displacement regime for this short-span class. E.g. if many low bits agree, the ordinary displacement can still be macroscopic while the span estimate applies. Highly unbalanced $(W,B)$ can be treated out toward the $M^2=N$ boundary using the other factor in (2). The next argument strengthens this balanced threshold. No statement here covers the full long-span collision sector, or proves positive density.


### Better balanced bound by slope-adapted splitting

For any fixed primitive absolute slope $D/A$, the number of core presentations with $(x-y,|\operatorname{rev}_s(x)-\operatorname{rev}_s(y)|)=g(A,D)$ for some $g$ is $O(6^{s/2})$, uniformly. This estimate is only for one primitive slope; it does not assume that the signed-digit weights are equidistributed over slopes.

For the proof choose a digit split at $m$, radices $B_*=2^m,C_*=2^{s-m}$ with $C_*/B_*$ within a factor four of $D/A$. This is possible because the possible radial ratios range from $1/M$ to $M$ in steps of four, and $A,D\le M$. Let $(u,u^*)$ be the ordinary and reversed difference signature of the two low blocks of length $m$, and let $(v,v^*)$ be the analogous signature of the high blocks. For each sign of the reversed difference, the slope condition is
$$
D(u+B_*v)=\pm A(v^*+C_*u^*).
$$
Given $(u,u^*)$, this confines the integer $v$ to an interval of bounded length: the uncertainty from $v^*$ is at most $A C_*/(D B_*)=O(1)$. Then $v^*$ is fixed by the equation. Reversing the argument, given $(v,v^*)$ the integer $u^*$ has $O(1)$ choices because $D B_*/(A C_*)=O(1)$, and then $u$ is fixed. Thus the valid signature pairs form a bipartite graph of bounded degree. The squared multiplicity sum for the low signatures is $6^m$: each signed position has two nonzero possibilities of multiplicity one and one zero possibility of multiplicity two. This uses signature rigidity to identify distinct signed vectors. Similarly the squared multiplicity sum for the high signatures is $6^{s-m}$. Cauchy--Schwarz on the signature edges gives the claimed $O(6^{s/2})$.

For primitive directions with $\max(A,D)\le R$, the geometric-mean contribution to the lattice count is at most
$$
O\bigl(\sqrt{WB}\,6^{s/2}\sum_{1\le A,D\le R}(AD)^{-1/2}\bigr)
=O\bigl(\sqrt N\,M^{(\log_2 6-1)/2}R\bigr).
$$
For directions with $\max(A,D)>R$, the primitive second-moment tail satisfies
$$
\sum_{\max(A,D)>R}\frac{1}{AD\max(A,D)}
\ll \frac{\log(2R)}{R}.
$$
To see this, sum first over $D\le A$ to get an $O(\log(2A))$ numerator for each denominator $A^2$, and double for the other side. Restoring $g$ gives an extra factor of at most $M$ before taking the square root. With $P=2$, the full-core signature weight norm is at most $6^{s/2}$, so Hölder now estimates this tail by $O(\sqrt N\,6^{s/2}\sqrt{\log(2R)/R})$. Choose $R=M^{1/3}$, and include the unrestricted $M^2$ bound for the additive 1 in the lattice count. The final balanced bound at fixed $(s,j)$ is
$$
O\left(M^2+\sqrt N\,M^{\log_2(6)/2-1/6}\sqrt{\log(2M)}\right).       \tag{3}
$$
Thus all differing spans at most $(0.44412286644879795-\epsilon)L$ contribute $o(N)$ for any fixed positive $\epsilon$: the exact threshold is $1/(\log_2 6-1/3)$. Bound (2) may still be better for an unbalanced common prefix and suffix. Neither bound controls the full long-span sector.

## 3. Collision orientation and persistent primitive directions

**The physical rows have no negative reversed displacements.** For $u>u'$ with decreasing full reversal, the first differing low bit is zero in $u$ and one in $u'$. This sign of the reversal also equals that of the two-ended core $d$: the higher and lower shared digits cancel upon subtraction. Write $u=2^j X+t$, $u'=2^j Y+t$ with $X>Y$, $X$ even and $Y$ odd. From the checksum recursion and $0\le p(n)\le n$,
$$
p(X)-p(Y)\ge (X-Y+1)/2\ge 1.
$$
If $j\ge1$, the shared-low split adds a nonnegative multiple of $X-Y$ and at worst a correction of $-1$. Then $p(u)\ge p(u')$, the two parity readouts agree, and $G(u)>G(u')$. If $j=0$, use that $p$ is always even by its recursion: its positive integer difference is at least two, whereas the parity term is $-\phi$. Therefore $G(u)>G(u')$ whenever $\phi<c+2k$, which includes both physical rows. This threshold is sharp uniformly in $c,k$: when $\phi=c+2k$ the inputs 2 and 1 collide and their reversal decreases.

For arbitrary fixed $\phi$, negative collisions are confined to an explicitly bounded initial range. Only $j=0$ can occur; write $u=2z$, $u'=2z'+1$ with $D_0=z-z'\ge1$. The exact equality is
$$
\phi=c(2D_0-1)+k[p(z)+D_0+(z\bmod2)+z'-p(z')].
$$
Thus if $C_0=\lfloor\phi/k\rfloor<1$ none occur; otherwise $D_0,p(z),z'-p(z')\le C_0$. For $z'>0$ write $z'=2^a m$ with $m$ odd. The deficit $n-p(n)$ for $n=2v$ equals the deficit at $v$ minus $v\bmod2$; at odd $m$ it is at least $(m+1)/2$. Hence $m\le2C_0+1$. Choose the smallest dyadic $K>C_0+1$. If $2^a\ge K$, then $z\bmod K=D_0$ is not the all-ones word, and a zero digit below that cutoff contributes at least $\lfloor z/K\rfloor$ by the recursion. It follows that $z<K(C_0+1)$. If $2^a<K$, then $z'\le(K/2)(2C_0+1)$ and $z\le z'+C_0$. The case $z'=0$ already has $z=D_0\le C_0$. This proves a finite bound $O((1+\max(\phi,0)/k)^2)$, independent of the dyadic length $L$.

**Positive-slope collisions persist.** For a fixed core with $d>0$ and $\alpha=(1+c/k)a$ integral, fix $w\bmod2=e_w$ and $r\bmod2=e_r$. Extend by a common low suffix of length $j\ge1$. Write $\beta=-a+p(x)-p(y)-e_w+\varepsilon(1-e_r)$. Substituting $w=e_w+2v$, $r=e_r+2z$ into (1) gives an integral lattice iff
$$
 \alpha 2^j+\beta-d e_w-a e_r\equiv 0\pmod{2g}.
$$
On each admissible $j$, the valid $r$ are one residue modulo $2d/g$ inside $[0,2^j)$; every such solution has $w>0$ once $j$ is large, since $\alpha>a$. This gives exactly $(g/(2d))2^j+O(1)$ collisions for these parities at length $j$, on addresses of order $4^j$. For the 3-cycle row take $x=1,y=0$, $w,r$ both even. The identity is $w=7\cdot2^j-r$ for every even $r$ in the interval. The mask-zero example has $x=13,y=2,s=4$, $a=11,d=7$, and its even-even identity is $7w+11r=200\cdot2^j-8$. The congruence fixes one class of $r$ modulo fourteen. These actual infinite collision families saturate the square-root scale of a fixed core but remain far below linear energy.

## 4. Exact inversion for a prescribed physical carrier output

For the mask-zero carrier
$$
F(x,y)=100xy+89x+11\big(p(x)+p(y)+1+(x\bmod2)\big),\qquad x,y\ge1,
$$
there is a deterministic search that returns every exact pair with $F(x,y)=T$, or certifies that this carrier has none. It takes $O(\sqrt T\log^2(2T))$ elementary digit and arithmetic steps on $O(\log T)$-bit words and $O(\log T)$ working bits when solutions are streamed.

1. Enumerate the smaller coordinate $a\le\sqrt{T/100}$, with both orientations, counting the equal pair only once. If $x=a$, the remaining variable $n=y$ satisfies $c=100a,\phi=0,Q=T-89a-11(p(a)+1+(a\bmod2))$ in $cn+11p(n)+\phi(n\bmod2)=Q$. If $y=a$, then $n=x,c=100a+89,\phi=11,Q=T-11(p(a)+1)$.
2. The elementary bound $0\le p(n)\le n$ gives $n\in[(Q-\phi)/(c+11),Q/c]$, and require $n\ge a$. This interval has length $O(T/a^2+1)$ and its endpoints are on scale $T/a$ for sufficiently large $T$, since $a=O(\sqrt T)$ and $Q=T-O(a)$.
3. Choose a power of two $K\asymp\sqrt T/a$ with $K\ge2$. For each possible prefix $w=\lfloor n/K\rfloor$ write the inherited split as
$$
Q=L(w)-11w\operatorname{rev}(b)+R(b),\quad
L(w)=(cK+11(K-1))w+11p(w),
$$
where $R(b)=cb+11p(b)+\phi(b\bmod2)+11[b<K/2](w\bmod2)$ and $0\le R(b)\le(c+11)(K-1)+\phi+11$. The interval for $n$ leaves only $O(\sqrt T/a+1)$ prefixes. Because $w\asymp\sqrt T$ and the range of $R$ is $O(aK)=O(\sqrt T)$, at each prefix the possible reversals lie in an integer interval of bounded length. Reverse each candidate, then check the literal equality. If $w=0$, check its entire suffix interval directly; it can occur only in the bounded small-$T$ regime.

Summing $O(\sqrt T/a+1)$ over all $a$ gives $O(\sqrt T\log(2T))$ prefixes. Each individual evaluation takes $O(\log(2T))$ elementary digit steps. Physical realizability for each returned pair is the inherited carrier construction; here the target itself restricts both the absolute prefix and the reverse suffix before the physical host is selected. This is exact membership for this carrier, not a proof that every admissible odd target belongs to it. The same algorithm works with any fixed positive leading product coefficient and nonnegative bounded-coefficient linear/checksum/bit terms of the two coordinates.

The script `verify_twisted_inverse.py` independently compares (1) with direct digit iteration and compares the decoder's full output fibers against direct enumeration on bounded controls.
