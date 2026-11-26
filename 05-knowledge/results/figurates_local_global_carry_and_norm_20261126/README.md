# Polygonal versus polyhedral sums: carry, norm, and window

This continues the [joint Pollock moment test and exact selector](../pollock_joint_local_moments_and_selector_20261126/) and [June additive-basis regimes](../../../07-reflections/goldbach-polygonal-zeckendorf-additive-bases-s501.md). New results here: a minimal four-contact odd carry, a family of actual five-tetrahedral sums with size $\asymp X/\sqrt{\log X}$, and an infinite locally soluble but globally impossible **marked** five-term family. No unrestricted Pollock bound is claimed.

## 1. What the degree ladder actually transfers

Write $P_m(x)=x+(m-2)\binom{x}{2}$ and $F(x)=(x^3-x)/6$; nonnegative tetrahedral numbers are $F(x)$ for $x\ge1$. The [polygonal ladder](../../../07-reflections/bsd-hodge-the-polygonal-ladder-and-what-merges-kps7.md) and [Vandermonde truncation](../../../07-reflections/the-vandermonde-truncation-law-polygonal-vs-polyhedral-opus-S317.md) meet in the exact layer identity
\[
\sum_{i=1}^x P_m(i)= (m-2) F(x+1)+(3-m)P_3(x).
\]
Triangular layers give exactly one tetrahedron. Square layers give $2F(x+1)-P_3(x)=F(x)+F(x+1)$. This is a stacking relation, not an equivalence of sum problems. For a sum of polygonals, completing the square transforms each term by
\[
8(m-2)P_m(x)+(m-4)^2=(2(m-2)x-(m-4))^2.
\]
The quadratic lead is fixed, but the residues **and the nonnegative-index condition** on the recovered roots have to be retained.

The analogous cubic operation is to lock a pair's sum. Its exact correction is
\[
F(c+a)+F(c-a)=2F(c)+ca^2 .
\]
Thus the quadratic coefficient now **moves with the center**: choosing centers is choosing a quadratic form and its coverage window together. Three pairs plus one residue-adjusting term have room for a ternary form. This is the mechanism behind [Brady's large octahedral result](https://arxiv.org/abs/1509.04316), not a proof with only two pairs plus one term. The octahedral polynomial itself is $4F(x)+x$. The more recent [Basak–Dong–Saettone–Zaharescu results](https://arxiv.org/html/2608.08019v1) likewise separate a large-number reduction from the finite bottom range: their exact all-number bounds are 15 icosahedral and 22 dodecahedral, whereas the original 13 and 21 have genuine small counterexamples. Large-variable smoothing, fixed arity, and constructive carry are distinct resources, as emphasized in the earlier archive.

## 2. Four contacts are enough, and necessary, for an odd carry

Let $t,a,b,d$ be positive integers. Take the four even corners and the four odd corners of the additive cube with starting index $t$ and directions $a,b,d$:
\[
\begin{split}
 E&=(t,t+a+b,t+a+d,t+b+d),\\
 O&=(t+a,t+b,t+d,t+a+b+d).
\end{split}
\]
Their first and second index power sums agree, and
\[
 \sum_{x\in O}F(x)-\sum_{x\in E}F(x)=abd .
\]
Indeed this is the third mixed difference: degrees below three vanish and the $1/6$ cubic coefficient gives exactly the product. Hence setting $(a,b,d)=(1,1,g)$ makes **any prescribed positive gap** $g$, with positive indices and both lower moments preserved. One can append the same fifth index on both sides. For $g=1$ this is the four-term finite-difference carry $F(t+3)+3F(t+1)-3F(t+2)-F(t)=1$.

Four is minimal for an odd gap under the lower-moment guard. For two equal-length triples with the same sum, square sum, and hence $e_2$, their products agree modulo 4:
* The number of odd entries is determined by the square sum mod 4. With at most one odd entry, the product is divisible by 4.
* With two odds, their product mod 4 equals $e_2$ mod 4 (the even times the sum of the odds contributes zero). That product determines the sum of the odds mod 4, so the total sum determines the remaining even entry mod 4.
* With three odds the total sum mod 4 determines the parity of the number of residues equal to $-1$, hence the product.

Newton's identity gives a difference of tetrahedral totals equal to half the difference of products. It is even. Padding both shorter tuples by the same index reduces to the triple case. Thus three contacts cannot realize an odd gap with those guards, even before imposing positivity.

The cube carry has a second description in the earlier [tetrahedral matching-cut coordinates](../tetrahedral_pair_refinement_20261004.md). For four indices let $(S,U,V,W)$ be the Walsh transform, with $S$ the sum and the other three coordinates the pair cuts; the rows are $(1,1,-1,-1)$, $(1,-1,1,-1)$ and $(1,-1,-1,1)$. With $R=U^2+V^2+W^2$,
\[
 \sum x_i^2=(S^2+R)/4,\qquad
 \sum F(x_i)=(S^3+3SR+6UVW-16S)/96.
\]
Flipping one cut sign changes the total by $\pm UVW/8$ while keeping the lower moments. If $S$ is even, the cuts are all even and the sign flip preserves the integer lattice condition. Both tuples are positive if $S\ge4+|U|+|V|+|W|$; if both are positive that inequality is also necessary since across the two orientations all eight sign triples appear. For odd $S$ a single cut flip is not integral. With $(U,V,W)=(2a,2b,2d)$ and $S=4t+2(a+b+d)$ the two orientations are precisely the cube corners, up to permutation. This extends the existing affine cut-coordinate map by identifying its cubic channel.

The gap is not the starting value. For fixed nonzero $g$ an integer sign-flip carry has only finitely many cut choices by factoring $|g|$; both the center and any common fifth index are $O(X^{1/3})$ when totals are at most $X$. Thus a fixed gap from this sign-flip family can start at only $O_g(X^{2/3})$ targets up to $X$. This is a direct analogue of the [Chan–Pak fixed-address jump versus source distinction](../five_core_web_20261111/fivefold/DYADIC_SHELLS.md), now with a minimal contact count and an absolute sparsity bound. It is not a bound for all representations with the same gap.

## 3. A five-summand family that moves the mean and avoids repetitions

For each dyadic $c=2^k\ge2$, take two symmetric pairs about $c$ and a fifth term at $2u$:
\[
 N=4F(c)+F(2u)+cq,\qquad 1\le u\le c,\qquad
 q=a^2+b^2\le(c-1)^2 .
\]
Use indices $c\pm a,c\pm b,2u$, so all five are positive. The actual mean is not fixed when $u$ changes. The polynomial $F(2u)=(4u^3-u)/3$ permutes residues mod $c$: its difference at $u,v$ is $(u-v)[4(u^2+uv+v^2)-1]/3$, with odd bracket and odd denominator. Thus $N\bmod c$ chooses the fifth index uniquely. Recover $q$ by division and find $a,b$ by a square lookup or apply the two-square criterion. There are **no collisions within the dyadic block** across $(u,q)$, even if $q$ has several square decompositions.

If $R_2(Y)$ counts nonnegative sums of two squares at most $Y$, the block has exactly $cR_2((c-1)^2)$ distinct integer totals. They're all below $3c^3$. Taking the largest dyadic $c\le (X/3)^{1/3}$ gives $\gg X/\sqrt{\log X}$ actual targets up to $X$, by [Landau's two-square asymptotic](https://pure.mpg.de/pubman/item/item_3122284_3/component/file_3122285/brink_principal_oa._2011.pdf). Conversely a block occurring up to $X$ has $4F(c)\le X$, so $c=O(X^{1/3})$. Summing the block upper bounds over dyadic centers gives $O(X/\sqrt{\log X})$; the smaller centers at most $X^{1/6}$ contribute only $O(X^{1/2})$ by the trivial block bound. Thus this **restricted union**, not the unrestricted Pollock set, has order exactly $X/\sqrt{\log X}$. It gives a global construction far larger than the fixed-gap carries, but still density zero.

## 4. Unequal centers can pass all local tests and still fail globally

For two unequal paired centers $m,23m$, with fixed fifth index 1, the target
\[
 N_j=2F(m)+2F(23m)+m\cdot3\cdot25^j,\qquad m=2\cdot5^j+1
\]
asks exactly for
\[
 a^2+23 b^2=3\cdot25^j,\qquad |a|\le m-1,\quad |b|\le23m-1.
\]
It has a real solution in the positive window ($a=\sqrt3\,5^j,b=0$), and a solution over every $\mathbb Z_p$. For the latter it is enough to solve right-hand side 3 and multiply both coordinates by $5^j$. At 2 use $b=5$ and $a^2=-572=4(-143)$, whose unit part is $1\bmod8$. At 3 use $b=1$ and lift the square root of $-20$ from mod 3. At 23 lift $7^2\equiv3$ with $b=0$. At every other prime the norm form over $\mathbb F_p$ is split or an $\mathbb F_{p^2}$ norm; it represents every nonzero residue regularly, so Hensel applies.

There is no integer solution: for $j>0$, modulo 5 the nonsquare $-23$ forces $a,b$ both divisible by 5. Deflate repeatedly to $a^2+23b^2=3$, impossible by positivity. This is a failure of the **marked chart**, not a Pollock counterexample. For example $N_0=109497=F(1)+F(18)+F(43)+F(83)$, with room for a zero fifth tetrahedron. Removing the pair marking repairs this first target.

The practical next test is whether a *movable* quotient in §3, with several independent center choices, is guaranteed to hit a globally represented binary norm within its positive window. Section 4 says to test actual norm representability, not just congruences even to every modulus; §2 says a supply of small exact displacements does not by itself fill missing starting totals.

Replay: `python replay.py`. It checks 4,096 mixed carries, the triple obstruction mod 4, exact dyadic multiplicities through $c=256$, and the marked failure with an unmarked witness.
