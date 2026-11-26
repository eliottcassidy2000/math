# The only spine reversal invisible to the ordinary and cover counts

2026-11-26. A follow-on to [mixed spine reversals](TRANSFERRED-BOUNDARIES.md) and the [unique-terminal grammar](UNIQUE-TERMINAL-GRAMMAR.md).

**Theorem.** In a binary rooted left arm with its transitive spine, label the sites $0,1,\ldots,n-1$ and reverse the physical arc $i\to j$ with $i<j$. Let $F$ count spanning paths, $E$ spanning paths starting at the root, and $R$ ordered two-path covers whose first path ends at the root (the second may be empty). If $\Delta R=0$ and $\Delta E\ne0$, then for an integer $i\ge1$ the only possibilities are

    j = 2i+1,  address = 2^(2i+1)+1,
    (ΔF,ΔE,ΔR) = (0,2^i,0).

In particular the address-9 reversal $1\to3$ is the **only** change of this type that produces a unit jump with any right core. For the other members of the family the exact jump is $-2^i e$, where $e$ is the right terminal count.

**Proof.** For adjacent sites there is no change in the total number of paths on any spine subset. If $i>0$ the first site of the spanning spine path remains 0, so $\Delta E=0$. If $i=0$ the only loss or gain in $R$ is the change of last site for the subset $\{0,1\}$, giving the same difference $x_0-x_1$ as $\Delta E$.

For a nonadjacent pair write $m=j-i-1\ge1$, $A=2^m$, and $x_k$ for the root bit at site $k$. The mixed-subset expansion gives
$\Delta R=\sum_{S\supseteq\{i,j\}}[\Delta d(S)+\Delta k(S)d(T)]$,
where $T$ is the complementary subset, $d$ counts paths ending at a 1-bit and $k$ counts all paths; $d(\emptyset)=1$. The changed paths on $S$ arise by splitting its selected interior sites between the prefix and suffix at the new arc. Thus $\Delta k(S)$ is positive if an interior site is selected and zero otherwise. Write $\Delta R=A+\sum_k B_k x_k$. The expansion gives these exact properties:

* The constant term is $A$, from $T=\varnothing$.
* The only negative coefficient is $B_j=-2^i$. The lost increasing paths ending at $j$ have no selected interior or upper site, and have arbitrary choices below $i$. No added path ends at $j$. At a higher existing endpoint the loss, if any, is canceled by the replacement; all other contributions are nonnegative.
* For $k<i$, $B_k=A2^k$. An edited path on $S$ never ends that low, so only $\Delta k(S)d(T)$ contributes. Having $k=\max T$ forces the entire interval into $S$, with arbitrary $T$ below $k$.
* For $k\ge i$, $2^i$ divides $B_k$. With the top-bit endpoint condition fixed, membership of the $i$ sites strictly below $i$ is free; neither a charged path's last site nor the number of added paths depends on those choices. For $k=i$, at least $2^i$ contributions end at $i$ using the new arc. For every $k>i$ other than $j$, the subset $T=\{k\}$ gives a positive contribution if an interior site remains in $S$. In the exceptional case $m=1$ and $k$ its only interior site, the added path on $\{i,k,j\}$ with order $j,i,k$ contributes instead. Thus every coefficient other than $B_j$ is strictly positive.

If $i<m$, then $A>2^i$ and $\Delta R$ cannot vanish. If $i>m$, reducing the expansion modulo $2^{m+1}$ gives $\Delta R\equiv A(1+x_0)$, so vanishing requires $x_0=1$. With $i\ge1$, however, all new full-spine paths have first site 0 and $\Delta E=A(1-x_0)=0$. The remaining case is $i=m$. The constant and the single negative coefficient then cancel only if $x_j=1$ and all other bits are zero. For this word, $\Delta E=A$. Both endpoints of the edit belong to $S$; the complement cannot contain $j$, so its old $a(T)=1$ and its old $d(T)=0$ except at $T=\emptyset$. The mixed formula for $\Delta F$ therefore gives $-2^i+\Delta E=0$.

**Corollary (uniformity is rigid).** If the signed count shifts from a fixed left spine reversal are merely bounded as the odd binary right address $c$ varies, then $\Delta F=\Delta R=0$. The shift is $c\Delta R+F(c-1)\Delta F-\Delta E$. On $c=2^a+1$ the value of $F(c-1)$ is exactly $c$, forcing $\Delta F+\Delta R=0$. On $c=2^a-1$ it is $2^{a-1}+1$, so the unbounded difference $c-F(c-1)$ forces $\Delta R=0$ as well. In particular the only everywhere-unit reversal is the address-9 edit or its inverse, even without assuming the signs agree. This does not rule out other fixed left edits making a unit shift on a thinner infinite set of right addresses.

The address-$9$ lane has midpoint 7 modulo 9, so it cannot cross an odd-exponent even perfect number, as in the preceding checkpoint. The classification explains why changing the span of this universal mechanism alone cannot fix that residue mismatch. To hit an exact midpoint one must allow a right-dependent shift, a nonbinary left core, or a different physical edit.

The [replay](isolated_spine.js), run with `node isolated_spine.js 16`, expands the profile change as a polynomial in the root bits. A subset zeta pass evaluates each polynomial on all $2^{16}$ root words for each of the 120 physical spine arcs. The seven nonzero isolated-root hits are the exact theorem family at indices $i=1,\ldots,7$.
