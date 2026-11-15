# Short-displacement collisions cannot account for the carrier energy

Consider the one physical mask-zero form
$$
 F(x,y)=100xy+89x+11\bigl(p(x)+p(y)+1+[x\text{ odd}]\bigr),
quad p(0)=0,\quad p(2w+1)=p(w),\quad p(2w)=p(w)+w+(w\bmod 2).
$$
All addresses below are positive integers. There is a new power-saving collision estimate in a range of displacement, not a total second-moment upper bound.

**Theorem.** Put
$$
 P=14/5,\quad Q=P/(P-1),\quad
 \gamma=\frac{\log_2(2+2^P)}{P},\quad \tau=1+2\gamma-2/P=2.5458626785\ldots.
$$
Fix $0<\delta<1/2$, and let $W\ge 1$, $W\le T^{1/(2\tau)}$. Within the central shell $T/2<F\le T$, $x,y\ge T^\delta$, the number of ordered distinct-address collisions whose addresses are within $W$ in both coordinates, **either in the given order or after transposing one presentation**, is
$$
  O_\delta\!\left(\bigl(\sqrt T\,W^\tau+W^4\bigr)\log^3(2T)\right).
$$
In particular this is $o(T)$ if $W\le T^{1/(2\tau)-\epsilon}$ for any fixed $\epsilon>0$; the threshold exponent is about 0.1964. The balance and the long thin parts of the central shell are both included. Thus this whole displacement range is below the $T\log^2T$ energy scale needed for density one; it cannot determine the remaining, longer-displacement contribution.

## 1. Splitting a close pair, including all high-bit carries

For a power of two $A=2^a$, put $R_A(b)=A-1-\operatorname{rev}_a(b)$, with $0\le b<A$. The exact split is
$$
 p(Au+b)=p(u)+p(b)+R_A(b)u+[b<A/2](u\bmod2)
 \qquad (A\ge2).
$$
For $A=1,b=0$ the identity is simply $p(u)=p(u)$, and we take $R_1(0)=0$ and the extra parity term zero. The split follows by $a$ applications of the defining recurrence.

Take a base radix $V=2^r\ge1$, and consider two integers differing by at most $V$. Their quotients on division by $V$ either agree or differ by exactly one. In the agreeing case take $A=V$, with arbitrary base suffixes $b,b'\in[0,V)$. If the quotients increase by one, let $j\ge0$ be the number of trailing 1s of the smaller quotient. Its low $j+1$ bits form the integer $2^j-1$, while those of the larger quotient form $2^j$. The remaining common prefix is unrestricted. Thus $A=V2^{j+1}$ and
$$
 (x,x')=(Au+B,Au+B'),\quad
 B=(2^j-1)V+b,\quad B'=2^jV+b',
$$
where the remaining inequality $|B'-B|\le V$ is retained whenever needed. The decreasing case exchanges these formulas. Every close pair has one of these descriptions, with at most two orientation types per larger radix and one type at the base radix. There are $O(\log(2T))$ types when the coordinates are $O(T)$. Only the $r$ base suffix bits can vary within each type.

Apply this to both coordinates with base radices $V_x,V_y$. Let
$$
 (x,x')=(Au+B,Au+B'),\quad (y,y')=(Cv+D,Cv+D'),\quad
 h=B'-B,\ k=D'-D,
$$
and put $d_x=\operatorname{rev}_a(B')-\operatorname{rev}_a(B)$, $d_y=\operatorname{rev}_c(D')-\operatorname{rev}_c(D)$.
On each of the four parity classes of the common prefixes $(u,v)$, the exact equation $F(x',y')=F(x,y)$ is an affine line with coefficients
$$
 \alpha=100Ak-11d_x,\qquad \beta=100Ch-11d_y.
$$
Both are nonzero away from the address diagonal: if a shift multiplying 100 is nonzero it dominates the reversal term, and if it is zero the reversal term is nonzero because the other shift is nonzero. The terms $89x$ and the explicit parity of $x$ contribute only to the constant, as do the two binary-split parity terms after restriction to a parity class. No parity constraint is used to infer a nonzero slope.

For a near transpose, relabel the close comparison as ((x,y)) and ((x',y')) but use the equation (F(y',x')=F(x,y)). The slopes are (alpha=100Ak-89A-11d_x), (\beta=100Ch+89C-11d_y) instead. They never vanish: (min_{kinmathbb Z}|100k-89|=11) but (|d_x|<A), and the other slope is analogous. The two offsets (-89A,89C) are fixed in each carry type, so the injective coefficient decoding below is unchanged. This allows pairs that are close only after swapping, which can be arbitrarily far in their original coordinates in the long thin region.

## 2. Signed-digit rigidity and its exact moment

For two $a$-bit suffixes $(B,B')$ write the signed difference of each bit as $s_i\in\{-1,0,1\}$. The pair of ordinary and reversed differences
$$
 h=\sum_{i=0}^{a-1}s_i2^i,\qquad d_x=\sum_{i=0}^{a-1}s_i2^{a-1-i}
$$
uniquely determines the full vector $s$. Indeed, subtract two candidate vectors to obtain a polynomial with integer coefficients of absolute value at most 2 and roots 2 and 1/2. After dividing by its initial power of the indeterminate, a nonzero such polynomial would have even constant coefficient by its root 2, and therefore constant coefficient $2$ or $-2$. At 1/2 the absolute sum of the remaining terms is strictly less than 2, a contradiction.

Inside a quotient/carry type the extra suffix bits are fixed, so the multiplicity of a derivative signature among the $r$ free low bits is bounded by $2^{\#\{i:s_i=0\}}$. Summing its $P$-th power over signatures gives at most
$$
  \sum_{s\in\{-1,0,1\}^r}2^{P\#\{i:s_i=0\}}
    =(2+2^P)^r.
$$
This bound remains valid after imposing a maximum displacement. Each coefficient pair $(\alpha,\beta)$ determines the four integers $h,k,d_x,d_y$: multiples of $100A$ are separated by more than the range of $11d_x$, and likewise with $100C$. Thus in fixed carry types the $\ell^P$ norm of the number of suffix choices over coefficient pairs is at most $(V_xV_y)^\gamma$. This bound does not count collisions; it is a moment of the possible line weights before imposing their intercepts.

## 3. Counting lattice points on the coefficient lines

Let both $x,x'\asymp X$ and $y,y'\asymp Y$, with bounded implicit comparison constants and $X\le Y$. Write $\rho=Y/X$. If $x'=x+h,y'=y+k$, from $x'y'-xy=h y'+kx$ and $0\le p(z)\le z$, the collision identity forces
$$
  |h|\ll 1+W/\rho.
$$
Indeed, the product change in absolute value dominates $100(|h|y'-|k|x)$; the absolute change from the other terms is $O(|h|+X+Y+1)$. Use base radices $V_x\asymp 1+W/\rho$ and $V_y\asymp W$, large enough for the displacement bounds. In the transposed comparison the asymmetric linear term is $O(Y)$ instead, so the same short-coordinate displacement bound holds. If the short coordinate is $y$, exchange the labels $(x,y)$ in this argument; the extra asymmetric $89x$ then has size $O(W)$, harmless since $W\ll Y$ in the shell ranges below.

Assume first $1\le\rho\le W$, so $V_x\asymp W/\rho$. For fixed carry radices $(A,C)$, $|\alpha|\ll AW$ and $|\beta|\ll CW/\rho$. For a fixed suffix choice an affine line in the prefix rectangle has at most
$$
 O\!\left(1+\frac{Y}{AC}\cdot
       \frac{g}{\max(|\alpha|/A,\rho|\beta|/C)}\right),
 \quad g=\gcd(|\alpha|,|\beta|),
$$
integer points, since consecutive integer points change $u$ by $|\beta|/g$ and $v$ by $|\alpha|/g$. The bound applies without any congruence restriction and hence also on each parity class.

For the $Q$-th moment of the second term over nonzero coefficient pairs, write their absolute values as $(g a,g b)$, with coprime positive $(a,b)$. Put $M=\max(a/A,\rho b/C)$; the number of $g$'s is $O(W/M)$. Also $M\ge m_0=1/\min(A,C/\rho)$, and at most $(AC/\rho)t^2$ positive pairs have $M\le t$, even dropping coprimality. A dyadic sum of $M^{-Q-1}$ bounds the moment by
$$
 \ll_Q \left(\frac{Y}{AC}\right)^Q W\,
         \frac{AC}{\rho}\,m_0^{1-Q}
   = \frac{X^Q W}{\max(A,C/\rho)^{Q-1}}.
$$
This includes all four sign choices at a constant multiplicative cost and converges since $Q>1$. Hölder against the suffix multiplicity moment bounds the nonconstant portion of the collision sum by
$$
 \ll_P (V_xV_y)^\gamma X W^{1/Q}\max(A,C/\rho)^{-1/P}.
$$
Sum over $A=V_x2^i,C=V_y2^j$ and their finitely many orientations: $\sum_{i,j\ge0}2^{-\max(i,j)/P}<\infty$. Thus this portion is
$$
  O\bigl(X W^{\tau}\rho^{1/P-\gamma}\bigr)
  \le O\bigl(\sqrt{XY}\, W^\tau\bigr).
$$
The additive 1 in the line count contributes $O((V_xV_y)^2\log^2(2T))$, by the number of possible base suffix pairs and carry classes.

If $\rho>W$ then the short-coordinate radix is $O(1)$, while the other radix is $O(W)$. There are $O(W^2)$ suffix choices per pair of carry types. Since both coefficients are nonzero, fixing the shorter-coordinate common prefix leaves at most one opposite prefix, and there are $O(1+X/A)\ll 1+X$ such prefixes. Hence the total is
$$
  O(W^2(1+X)\log^2(2T))
       \le O((W^2+\sqrt{XY}\,W^{3/2})\log^2(2T)).
$$

## 4. The dyadic boxes really cover all the close collisions

The recurrence gives $0\le p(z)\le z$ by induction. For $x,y\ge T^\delta$,
$F(x,y)=100xy+O(x+y+1)=100xy(1+O(T^{-\delta}))$.
Thus both $xy$ and $x'y'$ are in a fixed annulus comparable to $T$. Assign each collision using the dyadic box of $(x,y)$; there are $O(\log(2T))$ such boxes. Relabel the primed coordinates first if the pair is close only after transpose. Its longer coordinate is $\gg\sqrt T$, so the corresponding primed coordinate is within a constant factor (in fact $1+o(1)$) because $W\ll\sqrt T$. Comparing the products in the annulus then puts the other primed coordinate within a fixed factor as well. The argument of §3 was for arbitrary fixed-factor enlarged rectangles and thus applies to these primed coordinates, regardless of whether one crosses the original dyadic endpoint. Summing its bounds gives the asserted $O((\sqrt T W^\tau+W^4)\log^3(2T))$. Since $4/(2\tau)<1$, both terms are $o(T)$ under a fixed exponent margin.

The estimates concern **off-diagonal** pairs of presentations, independently of any orbit quotient or selection of allowed rows. The obstruction to completing the density-one proof is now demonstrably nonlocal in the underlying two address variables.

[Finite replay](verify_collisions.py): literal binary checks of the signed derivative weights and of the affine equation (80,000 generated quadruples plus exhaustive small carry classes). It also verifies the exact near-transpose coefficients. This is a control and is not used for any asymptotic estimate.
