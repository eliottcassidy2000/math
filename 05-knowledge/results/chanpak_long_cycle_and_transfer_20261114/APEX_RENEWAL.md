# One apex cannot supply all absolute targets

Let $W=w_1\cdots w_L$ encode one apex $v$ against a transitive spine in forward order, with $w_i=1$ for $v\to i$. Write $H(W)$ for the number of Hamiltonian paths. All one-apex/totally-ordered-spine tournaments are included; arbitrary edge reversals involving the apex and any contact-word restrictions are subfamilies.

**Theorem.** Let $A(X)$ be the number of distinct positive Hamiltonian path counts at most $X$ in this whole family, without a bound on tournament order. Then
$$
A(X) \le \frac{\log 2}{2}X+O(\sqrt{X}).
$$
In particular this entire family, not just a fixed-gap slice of it, misses at least $1-\log 2$ (about $30.7\%$) of positive odd targets asymptotically (in the liminf sense). This is an obstruction to absolute coverage for one-spine/one-apex controls, not for general tournaments or higher-interface controls.

## 1. Exact contact-word count and normalization

Deleting the apex from a path leaves at most two increasing runs in the spine. For an interior apex, let $i$ be the first vertex after the apex and $j$ the last vertex before it, so $w_i=1,w_j=0$. If $i<j$, every spine position smaller than $i$ is forced into the run before the apex, every one larger than $j$ is forced after the apex, and $j-i-1$ positions freely choose either run. These inversions contribute $2^{j-i-1}$. If $i>j$, a path is possible only for $i=j+1$, in which case it is the ordinary insertion of $v$ in the natural order. All ordinary insertions, including the two endpoints, count as $1+\#10$: the number of consecutive $01$ transitions plus the events “first bit is 1” and “last bit is 0” is $1+\#10$. Thus
$$
 H(W)=1+P_W(2)+[z^0]P_W(z),\qquad
 P_W(z)=\sum_{i<j}w_i(1-w_j)z^{j-i-1}.
$$
It follows that all counts are odd (the distance-one term is doubled). Leading zeros and trailing ones affect no inversion and can be stripped. Without an inversion the count is $1$; otherwise the normalized word is unique, starts with 1, ends with 0 and has length $L\ge 2$. This normalization can remove source and sink vertices without changing the count.

A useful exact implementation, for the intermediate prefix with integer value $p$, is
$$
 H(2p+1)=H(p),\qquad H(2p)=H(p)+p+(p\bmod 2),
 \qquad H(0)=1.
$$
The second increment is the sum of the newly created inversion weights into the terminal 0; the least significant bit's adjacency receives the extra copy. These are word-building recurrences, not a target inverse preserving existence.

## 2. Product limit with a quantitative error

For a normalized word set
$$
U_L=\sum_{i=1}^L w_i2^{1-i},\qquad
V_L=\sum_{j=1}^L(1-w_j)2^{j-L}.
$$
Each is in $[1,2]$. The unrestricted product $2^{L-2}U_LV_L$ has summand $w_i(1-w_j)2^{j-i-1}$. Retaining $i<j$ yields $P_W(2)$; discarded $i\ge j$ terms sum to $O(L)$ by the geometric series; the adjacency correction and constant are also $O(L)$. Hence uniformly in words,
$$
 \frac{H(W)}{2^{L-2}}=U_LV_L+O(L2^{-L}).
$$

Choose $W$ uniformly at width $L$, and keep only the first $m=\lfloor L/2\rfloor$ digits of $U_L$ and the first $m$ digits of the reversed expression $V_L$. Their two supports are disjoint, so the truncated variables are independent uniform grids on $[1,2]$, of mesh $2^{1-m}$. Coupling each grid to one independent uniform random variable by adding independent mesh-sized jitter shows that the normalized $H(W)$ differs from the product of independent $[1,2]$ uniforms by $O(2^{-L/2})$ almost surely under this coupling. The latter product has bounded density, so the distribution functions differ by $O(2^{-L/2})$, uniformly in the threshold.

The product distribution is
$$
F(z)=\begin{cases}
0&z\le 1,\
z\log z-z+1&1\le z\le 2,\
z-3+z\log(4/z)&2\le z\le 4,\
1&z\ge 4.
\end{cases}
$$
For example integrate $\max(0,\min(2,z/u)-1)$ over $1\le u\le 2$. Consequently the number of normalized words at width $L$ with $H\le X$ is $2^{L-2}F(X/2^{L-2})+O(2^{L/2})$.

## 3. Cancellation across repeated scales

Write $X=t2^k$, $1\le t<2$. The first-to-last inversion alone gives $H\ge 1+2^{L-2}$, so $L\ge k+3$ is excluded. For widths $2\le L\le k$, $F=1$, whereas the last two widths give $F(2t)$ and $F(t)$. The full word count $M(X)$ therefore satisfies
$$
 M(X)=2^k\bigl(\tfrac{1}{2}+\tfrac{1}{2} F(2t)+F(t)\bigr)+O(\sqrt{X})
     =X\log 2+O(\sqrt{X}).
$$
The second equality is exact for the main term: insert the two branches of $F$ and cancel $t\log t$, leaving $t\log 2$. Although the individual width contributions have a dyadic phase, the sum does not.

Converse corresponds to reverse-complement, an involution on normalized words fixing $H$. Odd widths have no fixed words and even widths have $2^{L/2-1}$; summing over relevant widths is $O(\sqrt X)$. Thus the number of word orbits with $H\le X$ is $M(X)/2+O(\sqrt X)$. Distinct counts can only be fewer; the count 1 adds one. This proves the theorem.

**Consequence for the proof strategy.** Both the nonlinear balanced swap and every other apex-incident switch inside this construction live in a universe with strictly insufficient count capacity. One must change the physical interface or prove that another independent family covers its complement. Surjectivity on every dyadic residue or tight local gap control inside this universe cannot imply absolute coverage. The coefficient $\log 2$ comes from merging widths, not a fitted exponent or finite spectrum extrapolation.
