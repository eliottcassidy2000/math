# Guarded apex words: product defects, a ternary obstruction, and exact renewal

A guarded apex word is a tournament with distinguished apex $z$ whose deletion is given as an ordinal sequence of nonempty blocks. Its first block is a singleton $O$ beaten by $z$, and its last block is a different singleton $I$ beating $z$. Concatenation $WV$ keeps one common apex and adds every arc from the deletion of $W$ into the deletion of $V$; it remains strong. Internally a block is arbitrary.

Write $D$ for the Hamilton count of the apex deletion, $Q$ for its ordered two-slot cover count including empty slots, and $H$ for the Hamilton count including the apex. For the normalized positive transfer of the apex word, let $a,b$ be the one-ended start and end reads and write $A=Da$, $B=Db$. Concretely, $A$ sums $H((T-z)-V(P))$ over all nonempty directed $P$ in $T-z$ starting in $N^+(z)$; $B$ does so over paths ending in $N^-(z)$. All the following quantities are integers.

## Exact composition of the five reads

Expand Hamilton counts by disjoint directed odd-cycle packings, with weight two per cycle. For a cycle using the apex, the rest is an even-order path in its deletion, with initial vertex in $N^+(z)$ and last vertex in $N^-(z)$. The parity projection for even-order paths is the sum of the unsigned read and the signed read at sign $(-1)^{|V(P)|}$. The sum of signed weights over *all* paths, including the one empty path, is zero on any nonempty tournament: expand in cycle packings with all other vertices as singletons, and move any one item between the distinguished path slot and complement slot to cancel. All items are odd-order. When forming ordinal words, this is the zero middle diagonal of the signed transfer. Moreover the signed start read is zero for a word ending in guard $I$; any path starting earlier acquires the cancellative full signed slot sum on the terminal block, and none can start in $I$. The signed end read similarly vanishes for a word beginning with $O$.

The composition is
$$
\begin{aligned}
D_{WV}&=D_WD_V, & Q_{WV}&=Q_WQ_V,\\
A_{WV}&=D_WA_V+A_WQ_V, &B_{WV}&=B_WD_V+Q_WB_V,\\
H_{WV}&=D_VH_W+D_WH_V-D_WD_V+A_WB_V.
\end{aligned}
$$
The last term joins a path starting on the $W$ side to a path ending on the $V$ side. The signed counterpart vanishes by the guards. The preceding terms account for apex cycles confined to one side or no apex cycle. The first four formulas follow by independent path covers in the two modules (the path starting on the left may use an arbitrary, possibly empty, segment on the right).

## A positive mod-four defect and the exact ternary rule

Define $E_O$ as the sum of $H((T-z)-V(P))$ over **even-order nonempty** paths $P$ with both endpoints in $N^+(z)$, and $E_I$ the analogous quantity for two endpoints in $N^-(z)$. The vanishing signed one-ended reads imply
$$
A=H-D+2E_O,\qquad B=H-D+2E_I.
$$
For the first, the signed start read is zero, so $A$ is twice the even-order start read, and the paths ending in the other contact class account for $H-D$ by the odd-cycle expansion. Similarly for the second. Hence
$$
H(WV)-H(W)H(V)=
2E_O(W)(H(V)-D_V)+2E_I(V)(H(W)-D_W)+4E_O(W)E_I(V).
$$
Each $H-D$ is positive (guards give an apex triangle) and even; all $E$ are nonnegative. Thus guarded concatenation is supermultiplicative and preserves the product modulo four. Equality holds precisely if the apex has unique outneighbor on the left and unique inneighbor on the right: two same-type neighbors always give an even two-vertex path with positive complement count.

For the guarded triangle $C=OI$, $D_C=1$, $H(C)=3$, and $A_C=B_C=2$. Thus, exactly,
$$
H(WC)=3H(W)+4E_O(W),\qquad H(CW)=3H(W)+4E_I(W).
$$
Both $3H-2$ and $3H+2$ are excluded from either operation, with no restriction on internal blocks. This is not a no-go for arbitrary two-port substitution, an internal edit, or reorienting arcs across the ordinal cut. Repeated guarded products remain congruent to the product of their inputs modulo four. For singleton blocks, $2E_O=\#\text{ adjacent }OO+\sum_{i<j:W_i=W_j=O}2^{j-i-1}$: this is the sum over even selected subsets between each pair. A nonzero right ternary excess can be exactly four, attained precisely by two O's separated by one or two positions. The left statement is the reverse-converse version.

## Exact scalar recurrence when a guarded block repeats

Put $a=A/D$, $b=B/D$, $d=Q/D$ and $c=H/D-1$. Here $d>1$ since empty and full path slots both contribute. In $W^k$, an apex cycle spanning two copies picks an $a$ path at its first copy and a $b$ path at the last, with factor $d$ per intervening copy. Same-copy contributions add $kc$; there are no signed inter-copy terms. Summing over the distance between first and last copy,
$$
H(W^k)=D^k\left[1+kc+ab\,\frac{d^k-1-k(d-1)}{(d-1)^2}\right].
$$
This includes $k=0$ with $H=1$ from the apex alone. In particular its constant-coefficient recurrence has characteristic dividing $(x-Q)(x-D)^2$.

If $W$ consists entirely of singleton blocks, all three terms $Q^k,D^k,kD^k$ have strictly positive coefficients, so the characteristic is exactly $(x-Q)(x-D)^2$. Indeed $D=1$, $d=2^{|W|}$ and $0<a,b<d-1$ by summing singletons' place values and using both guard types. This proves positivity of the constant and exponential coefficients. For the linear coefficient, write $c=c_++c_-$ from the signed and unsigned cycle reads. Then $c_-$ counts adjacent $OI$; $ab/d=c_++R$ where $R$ sums $2^{-(j-i+1)}$ over reverse pairs $I_i,O_j$ with $i<j$. Each such pair spans an adjacent $IO$ boundary, and the sum of weights of all pairs spanning any fixed boundary is less than one, by two finite geometric sums. Consequently $R\le\#IO=\#OI-1=c_--1$. Also $c\le d-2$, since $c/2$ counts even nonempty selected subsets satisfying the endpoint constraints and excludes the empty subset. It follows that $c(d-1)-ab=d(c_--R)-c>0$, proving strict positivity of the linear coefficient.

## Removing one guard: exact endpoint sign and infinitely many strong $3H+2$ moves

The signed root read is actually a count with a universal negative sign:
$$
\sum_{P:\,\operatorname{first}(P)=v}(-1)^{|V(P)|}H(U-V(P))
=-\,\#\{\text{Hamilton paths in }U\text{ ending at }v\}.
$$
The reversed statement follows by taking converse tournaments.

Here is an algebraic proof. In the commutative ring over $\mathbb Q$ with $x_i^2=0$, let $X=\operatorname{diag}(x_i)$, let $A$ be the zero-diagonal adjacency matrix, let $\mathbf 1$ be the ones column, and put
$$
\mathcal H(x)=\sum_{S\subseteq V(U)} H(U[S])\prod_{i\in S}x_i
=\prod_i(1+x_i)\frac{\det(I+AX)}{\det(I-AX)}.
$$
The determinant identity follows by taking formal logarithms in the nilpotent ideal: even directed cycles cancel in the ratio, each odd directed cycle has log coefficient two, and the prefactor inserts singletons. Exactness follows by exponentiating, since only disjoint squarefree supports survive. The formal path series rooted at $v$ is
$$
\mathcal P_v(x)=x_v[(I-AX)^{-1}\mathbf 1]_v.
$$
Any walk revisiting a vertex vanishes by squarefreeness. Set $B=I-(A^\top+I)X=(I-A^\top X)(I-X)$. The tournament identity $A+A^\top=\mathbf 1\mathbf 1^\top-I$ gives $I+AX=B+\mathbf 1 x^\top$. A rank-one update gives $\operatorname{adj}(I+AX)\mathbf 1=\operatorname{adj}(B)\mathbf 1$ (apply the determinant lemma and inverse formula to $\mathbf 1$; all matrices are units at $x=0$). Also $\det(I-A^\top X)=\det(I-AX)$ by transpose and $\det(I-XA)=\det(I-AX)$. The determinants and $\prod_i(1+x_i)$ cancel in $\mathcal H(x)\mathcal P_v(-x)$ against $\det B$, because $(1+x_i)(1-x_i)=1$. Multiplication by $x_v$ ignores the extra factor $(I-X)^{-1}$ in the inverse of $B$. Consequently
$$
\mathcal H(x)\mathcal P_v(-x)
=-x_v[(I-A^\top X)^{-1}\mathbf 1]_v.
$$
Taking the coefficient of the full vertex set gives the identity.

Now allow the left apex tournament $W$ to be entirely arbitrary (its deletion can be viewed as one ordinal block), and still append the ordinal guarded triangle $C=OI$ with the apex shared. Define $R_O$ to be the number of Hamilton paths of $W-z$ ending anywhere in $N^+(z)$. The signed start read of $W$ is $-R_O$ by the identity. The signed end read of $C$ is zero, so the product formula for the $H$ coordinate remains $H(WC)=H(W)+2D+2A$. Accounting for even-order paths now gives $A=H-D+2E_O+R_O$, and hence
$$
H(WC)=3H(W)+4E_O(W)+2R_O(W).
$$
Thus $3H-2$ is impossible even without the guard on $W$, while $3H+2$ holds exactly when $z$ has a unique outneighbor $u$ and there is a unique Hamilton path in $W-z$ ending at $u$. (If no outneighbor exists then $R_O=0$; if at least two exist then $E_O\ge1$.) The reverse-converse operation gives the corresponding left rule.

This produces infinitely many distinct strong parent counts with actual strong $3H+2$ descendants. Let $A_n$ on $1,\dots,n$ have $i\to i+1$, and $j\to i$ whenever $j-i\ge2$. Its Hamilton path ending at $n$ is unique: the last predecessor is forced to be $n-1$, then the same argument applies to the remaining induced initial interval. It is strong for $n\ge3$ by its spanning directed cycle. Let $W_n$ add the apex beating $n$ only. Then $W_n$ is strong. Appending $C$ keeps it strong, because the apex reaches the suffix through its $O$, the suffix reaches the apex through its $I$, and the old tournament was already strong. The parent counts are unbounded: $H(W_n)\ge H(A_n)\ge 3^{\lfloor n/3\rfloor}$ from disjoint consecutive triangles.

The deletion family also has an exact independent renewal count. Every directed cycle in $A_n$ is a consecutive interval: from its minimum, the directed chain increases in increments one until the first backward arc, which must return to the minimum because all other smaller labels on the interval have already been visited. Hence disjoint odd cycles are the same as disjoint odd intervals of length at least three. Tiling by those intervals with weight two, and by unused singletons with weight one, yields
$$
\sum_{n\ge0}H(A_n)x^n
=\frac1{1-x-2x^3/(1-x^2)}
=\frac{1-x^2}{1-x-x^2-x^3}.
$$

More generally the shared-apex join is supermultiplicative even when neither word is guarded. Define $U_0=H(W)-D_W$, $V_0=H(V)-D_V$, $x=2E_O(W)+R_O(W)$, and $y=2E_I(V)+R_I(V)$, where $R_I$ counts old Hamilton paths starting in the apex's inneighborhood. The signed cross term in the composition is the product of the two nonpositive read values, hence it is $+R_O(W)R_I(V)$. Expanding the positive cross term $(U_0+x)(V_0+y)$ and subtracting $H(W)H(V)$ gives
$$
H(WV)-H(W)H(V)=xV_0+yU_0+xy+R_O(W)R_I(V)\ge0.
$$
If both constituents are strong of order at least three, $U_0,V_0$ are positive even: a strong tournament has each vertex in a directed triangle (otherwise the inneighbors all beat the outneighbors, contradicting strong connectivity). The product rises by exactly two precisely when $(x,y,V_0)=(1,0,2)$ or $(y,x,U_0)=(1,0,2)$. Indeed two nonzero defects already contribute at least four, and with a single nonzero defect the corresponding multiplication must be exactly two.

In fact, in the strong $+2$ classification the zero-defect, two-excess receiver must be the triangle. Suppose $y=0$ and $V_0=2$. Then $z$ has one inneighbor $t$. By the odd-cycle expansion there is exactly one odd cycle using $z$ (its complement has count one), necessarily the unique triangle $z\to u\to t\to z$. Let $X=N^+(z)$. The vertex $t$ beats every vertex of $X-\{u\}$, and $X-\{u\}$ is transitive because it is the complement of the unique cycle and has Hamilton count one. There is no directed three-vertex path within $X$ ending at $u$, since it would give a directed five-cycle through the apex. Thus either $u$ is first in the transitive ordering of $X$, or it has one predecessor which beats all of $X-\{u\}$ and precedes $u$ immediately. To see the two alternatives directly, two inneighbors would already produce the three-vertex path, and a unique inneighbor cannot be beaten by any other vertex of $X-\{u\}$. In the second case there is a deletion Hamilton path starting at $t$ followed by the transitive ordering of $X$, violating $R_I=0$. In the first case, every vertex outside the triangle belongs to a terminal ordinal suffix, violating strength unless there are none. This proves $V=C_3$. The opposite side follows by converse and interchanging factor order.
