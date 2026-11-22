# Complementary-digit cancellation and the exact scalar inverse of a two-flank slice

**Status: PROVED.** One absolute Hamilton count now inverts each prescribed signed-imbalance slice of the physical independent-flank construction (both source and destination), with no scale or displacement bound. The previously unresolved zero-imbalance chamber is included. Each unrestricted fixed slice gives exactly $\Theta_z(\sqrt X\log X)$ distinct physical counts up to $X$. No claim is made for a slice with unrestricted, unknown imbalance or for all of Chan–Pak.

## Inheritance and the exact target

Use the physical hosts in the [three-contact repair](../chanpak_free_seam_second_move_20261117/three_contacts/PROOF.md) and the notation of the [earlier inverse](../chanpak_joint_placement_20261117/inverse/PROOF.md). That inverse addresses a prescribed absolute count and shift, and, under a size and shift bound, handles an unknown shift except at zero imbalance. We address an absolute count directly on each whole imbalance slice. No extra surgery is being asserted in this note.

Let $G=2^k\ge4$, $x>0$ be even, $y>0$, $1\le r<G$ be odd, and $s=G-1-r$. Put $z=x-y$. On the source use
$
a=16x+5,\quad b=16y+4;
$
on the physical destination use $a=16x+9,\ b=16y$. In either case $u=Ga+r,\ w=Gb+s$. These addresses give the lower ordinary Hamilton count
$
N=81+142(u+w)+256uw+5(F(u)+F(w)),
$
where $F(0)=1,\ F(2v+1)=F(v),\ F(2v)=F(v)+v+[v\text{ odd}]$ for positive $v$ in the last rule. The usual reflected upper count is obtained by adding 2. The polynomial is the actual count given by that construction; below, no residue value is decoded as if it were an actual count.

Define $S=u+w=G(32x-16z+10)-1$, $d=u-w=cG+2r+1$, where $c=16z$ at the source and $16z+8$ at the destination. Then $N=T(S)-64d^2+5(F(u)+F(w))$ for $T(S)=64S^2+142S+81$.

## A cancellation in the two complementary words

Write $t=\operatorname{rev}_k(r)$, let $r_j$ be the $j$-th bit of $r$ from the least significant end, and set
$
J(r)=\sum_{j\ge0}(1-2r_j)\left(\left\lfloor\frac r{2^{j+1}}\right\rfloor+r_{j+1}\right).
$
This compensated coefficient depends on the actual suffix integer, not its padded length. The unsigned sum in this expression is exactly $r-r_0$: the floor-sum is $r-s_2(r)$, where $s_2$ is the binary digit sum, and the bit-sum is $s_2(r)-r_0$. In particular $|J(r)|\le r-1$ on our odd suffixes.

For any word $q$ of width $k$, the defining recurrence reads
$
F(q)-1=\sum_{j:q_j=0}\left(\left\lfloor\frac q{2^{j+1}}\right\rfloor+q_{j+1}\right).
$
The out-of-range bit is zero. In the complementary word $s$, a zero bit is a one bit of $r$, and for $j<k-1$
$
\left\lfloor\frac s{2^{j+1}}\right\rfloor+s_{j+1}
=2^{k-j-1}-\left\lfloor\frac r{2^{j+1}}\right\rfloor-r_{j+1}.
$
At $j=k-1$, when relevant, subtract one more because the bit beyond the word is zero on both words. Subtracting $t=\sum_j 2^{k-1-j}r_j$ therefore gives the identity
$
F(r)+F(s)-t+[r<G/2]=3+J(r)-2[r\ge G/2].
$
The usual concatenation rule from the same recurrence is
$
F(Ga+r)=F(a)+F(r)-1+(G-1-\operatorname{rev}_k(r))a+[a\text{ odd}][r<G/2].
$
The zero-bit sum proves it directly: the contribution from the prefix at suffix position $j$ is $2^{k-1-j}a$, with the only simultaneous parity correction at $j=k-1$. Apply it to each flank, using $\operatorname{rev}_k(s)=G-1-t$ and $a$ odd, $b$ even. Then
$
F(u)+F(w)=F(a)+F(b)+(G-1)a+1-c t+J(r)-2[r\ge G/2].
\tag{1}
$
In particular $c=0$ **cancels the entire high-scale reverse-digit contribution**, not merely its average. The remaining error is bounded by the unpadded suffix height $r$, rather than $G$.

## Absolute inverse for fixed signed imbalance

**Theorem.** For prescribed $z$, prescribed source/destination, and fixed $G$, at most one legal triple $(x,y,r)$ gives the supplied lower absolute count $N$. It can be found or excluded with $O(\log N)$ bounded-size integer comparisons and one rounded integer square-root computation, followed by a direct exact count check. Enumerating every $G$ is polynomial in the bit lengths of $N$ and $z$; there are at most $O(\log N)$ choices. Upper counts use $N-2$.

**Recover the height $S$.** The admissible values of $S$, as $x$ runs over even numbers with $x,y=x-z>0$, form a grid of step $64G$. For every candidate put $A=256(|c|+1)G^2$. By $|d-cG|<2G$ and $1\le F(v)\le v+1$, the exact count lies in its band
$
T(S)-64c^2G^2-A
\ \le N\le\
T(S)-64c^2G^2+A+5(S+2).
\tag{2}
$
Since $S=G(16(x+y)+10)-1$ and $x+y>|z|$, we have $|c|\le S/G$. The increase in $T(S)$ on a grid step is
$
T(S+64G)-T(S)=8192GS+262144G^2+9088G
>512(|c|+1)G^2+5(S+2).
$
So the bands are ordered and disjoint. A binary search on grid points $S\le N$ finds the sole possible $S$ or excludes it; this yields $x,y,a,b$.

**Recover the integer suffix.** Let
$
C=T(S)+5\bigl(F(a)+F(b)+(G-1)a+1\bigr).
$
By (1),
$
C-N=64d^2+5c t-5J(r)+10[r\ge G/2].
\tag{3}
$
For $c\ne0$ we have $|c|\ge8$; the sign of $d$ is the sign of $c$, and $|d|\ge(|c|-2)G$. The absolute correction to $64d^2$ is at most $5((|c|+1)G+1)$, so on a legal solution $(C-N)/64>0$ and
$
\left|\sqrt{(C-N)/64}-|d|\right|
\le \frac{5((|c|+1)G+1)}{64(|c|-2)G}<1/8.
$
For $c=0$, $d=2r+1\ge3$ and the correction has absolute value at most $5r+5=(5/2)(d+1)$. After division by $64d$ the error is at most $5/96<1/16$. Thus in every case exact nearest-integer square-root rounding determines $|d|$, and the sign determines $d$. Now $r=(d-cG-1)/2$ is the only possible suffix. Reject a failed integer, range, parity, count, or required displacement test. This final replay is crucial: rounding locates a possible witness, not a promise of existence.

The digit evaluation of $F$ and reverse bits uses time polynomial in the bit lengths. The physical addresses on both sides have $u,w\ge G$, so $N\ge256G^2$ bounds the number of scales. This covers the entire $z=0$ chamber without a shift bound; the other signed slices are global as well.

### Native bounded shifts without prescribing $z$

For the earlier native-range question, suppose $c_1G^6\le N\le c_2G^6$ and the physical displacement satisfies $|\Delta|\le c_3G^2$, for fixed positive $c_i$. The inherited displacement formula is
$
\Delta=15y-(16384G^2-5)z-B_r,\qquad
B_r=4096G^2+1004G+30+2048Gr+40\operatorname{rev}_k(r)
$
with $0<B_r<6144G^2$. If $z<0$, the formula bounds its magnitude by $(c_3+6144)/16383$. If $z\ge0$, use $N\ge256(16G)^2xy$ and $x\ge y$ to bound $y\le \sqrt{c_2}G^2/256$; then $z\le(c_3+15\sqrt{c_2}/256)/16383$. Thus testing a bounded list of $z$ completes the prescribed-single-endpoint inverse in this whole native branch, and recovers its shift by forward evaluation. Unlike the old nonzero-$z$ decoder this works even in the central chamber. Nothing here bounds the number of imbalances across completely unrestricted two-flank inputs.

## Exact height count for a whole slice

**Corollary.** For each fixed signed integer $z$ and fixed source/destination choice, with **no** displacement bound, the number of distinct lower counts $N\le X$ is $\Theta_z(\sqrt X\log X)$.

For the upper bound, $N\ge256(16G)^2 xy$ and $|x-y|=|z|$ imply $x=O_z(1+\sqrt X/G)$. The $O(G)$ suffixes give $O_z(\sqrt X)$ pairs at each feasible $G$, and $G\le\sqrt{X/256}$.

For the lower bound, use the subfamily
$
(|c|+2)^2G^2\le S/8.
$
Here $|N-T(S)|\le 64(|c|+2)^2G^2+5(S+2)\le23S$. Thus, even with $G$ unknown, the positive root of $T(q)=N$ rounds to $S$: movement by $1/2$ in $q$ changes $T$ by at least $64S+55$ on its lower side. Now
$
S+1=2G(16x-8z+5),\qquad v_2(S+1)=1+\log_2 G.
$
The count itself first recovers the height, then an independent exact valuation recovers the width. The fixed-$G,z$ theorem recovers $x,y,r$. Consequently all these outputs are distinct, across different scales as well as within one.

Indeed choose $S$ in $[\sqrt X/20,\sqrt X/16]$ on its grid and let every odd suffix occur. For all dyadic $4\le G\le \eta_z X^{1/4}$, with fixed sufficiently small $\eta_z>0$, these choices obey the subfamily condition, $x,y>0$ for large $X$, and $N<X$. There are $\asymp_z\sqrt X/G$ such grid points and $G/2$ odd suffixes, yielding $\asymp_z\sqrt X$ distinct values for each of $\asymp\log X$ scales.

This is an exact height-counting statement for one fixed-dimensional component, not positive density for the full spectrum. The structure is the arithmetic separation of an Archimedean height coordinate, a finite-place valuation, and an exact residual digit check. It is a literal instance of using both height and local information, inspired by the local-global and rank/height distinction in arithmetic geometry; no elliptic rank conjecture, Riemann Hypothesis claim, or equidistribution assumption enters. Here the complete height count can be proved because the local information is combined with an exact inverse and disjoint height bands, not merely because residue classes look compatible.

### A sub-square-root imbalance budget has density zero

One can allow the imbalance to vary with the target and still get a sharp obstruction to a near-balanced strategy. For $0\le Z\le \sqrt X$, the number of distinct absolute outputs up to $X$ with $|z|\le Z$, including **every width and suffix and both reflected tracks**, is
$
 O\!\left(\sqrt X\,(Z+1)
       \left(1+\log\!\left(2+\frac{\sqrt X}{Z+1}\right)\right)\right).
\tag{4}
$
In particular every $Z=o(\sqrt X)$ gives $o(X)$ absolute counts, even if we drop the restriction on the displacement completely.

To prove (4), at fixed $G,z$ write $v=\min(x,y)>0$, so $xy=v(v+|z|)$. The bound $256(16G)^2xy\le X$ puts the number of admissible $v$ at
$
O\!\left(\min\!\left(\frac{\sqrt X}{G},\frac{X}{G^2|z|}\right)\right)
$
when $z\ne0$, and $O(\sqrt X/G)$ when $z=0$. Multiplying by the number $G/2$ of suffixes and summing dyadic $G$, the contribution of a fixed $z$ is $O(\sqrt X(1+\log(2+\sqrt X/(1+|z|))))$, using the first term up to $G\asymp\sqrt X/|z|$ and a geometric tail afterward. Finally, with $Y=\sqrt X$ and $K=\lfloor Z\rfloor+1$, use
$
 \sum_{j=1}^{K}\log(2+Y/j)
 \le K\log(2+Y/K)+\sum_{j=1}^{K}\log(K/j)
 \le K\log(2+Y/K)+K.
$
This is a genuine height obstruction on the absolute image, not a modular deficit: enlarging the exact local read by any number of imbalances uniformly below the square-root scale cannot make this branch dense.

### Ordinal multiplication does not restore the missing height exponent

For any nonempty finite set of fixed imbalance values, let $\mathcal B$ be the union of its output banks (allowing both endpoints and their reflected partners). All counts here exceed 1. The exact count above gives abscissa of convergence $1/2$ for the nonnegative Dirichlet series $\sum_{b\in\mathcal B}b^{-s}$, with divergence at $s=1/2$ by partial summation. Ordinal sums multiply Hamilton counts; let $\mathcal M$ be the multiplicative monoid generated by the bank. For $s>1/2$,
$
 \sum_{n\in\mathcal M}n^{-s}
 \le \prod_{b\in\mathcal B}(1-b^{-s})^{-1}<\infty.
$
The inequality overcounts multiplicative factorizations and therefore makes no factorization uniqueness assumption. The abscissa for $\mathcal M$ is again $1/2$ since $\mathcal B\subseteq\mathcal M$; in particular the monoid has $O_\epsilon(X^{1/2+\epsilon})$ distinct counts through $X$. An exact inverse at fixed rank can support arbitrarily many CRT local tests and still fail to produce the absolute height growth necessary for a dense bank.
