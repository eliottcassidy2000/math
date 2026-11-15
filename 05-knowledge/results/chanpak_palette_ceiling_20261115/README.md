# The slope ceiling for matched tournament patches

**Proved, scoped to uniform marks on ordinal blocks.** This follows up
the [repeatable exact-digit construction](../chanpak_guard_depth_20261115/REPEATABLE_CLOCK.md)
and [physical mixed-radix transfer](../chanpak_outlet_blocks_20261115/MIXED_RADIX.md).
An independent background argument gives the ceiling for a fixed menu.
The weighted counting below also handles some branch-dependent menus
and arbitrary insertions of unmarked reservoirs.

The coordinates are an ordinary spine count $D$, a nonnegative marked
address $C$, the prescribed readout $N=H+2E$, and the last mark $z$.
For a patch (a nonempty sequence of the uniform marked blocks), on a
fixed incoming mark, its nonnegative integer transfer is
$$
D'=hD,\qquad C'=PC+\gamma D,\qquad
N'=hN+\beta C+\alpha D .
$$
Here $h=\prod_i h_i$, $P=\prod_i P_i$, with $h_i$ the ordinary count
and $P_i$ the count with two optional path slots. Always $P_i\ge2h_i$
because the paths putting everything in either single slot are included.
All admissible starting states below come from a fixed literal seed or
the empty spine.

## Every fixed matched palette misses density one

Let all choices in a fixed palette have the same $h,P,\gamma$ and
the same outgoing mark, which is the incoming mark for iteration.
The coefficient of the incoming $C$ is
$$
\beta=\sum_i(2-w_i)(P_i-h_i)
           \prod_{j<i}P_j\prod_{j>i}h_j .
$$
Indeed $C$ first passes through the preceding blocks, contributes to
$N$ in block $i$, and the contribution then passes through the later
blocks. Without $2-w_i$, the sum telescopes to $P-h$ (replace the
factors of $\prod h_j$ by $P_j$ in order). Consequently
$$
P-h\le\beta\le2(P-h).
$$
In particular the number $B$ of distinct integer slopes is at most
$P-h+1<P$ when $h>1$. This upper bound holds even if the palette
initially allows arbitrarily many different core profiles, locations
and marks. With $h,P$ fixed, the readouts give only finitely many
transfer descriptions anyway: there are at most $\log_2 P$ blocks,
and each pair $(h_i,P_i)$ is bounded.

If $h=1$, all the blocks are transitive. A block with one uniform
mark and order $\ell$ induces exactly the same transfer as $\ell$
singletons of that mark, from either incoming mark. This follows by
multiplying the singleton transfers, or from the mixed-radix formula.
After expansion the number of singleton marks is $\log_2 P$, and
$\gamma$ is the binary value of the expanded mark word. Hence a
matched palette with $h=1$ has just one transfer, not an exception
allowing $P$ slope choices.

For each fixed sequence of $k$ slopes, the $\alpha$ terms in $N_k$
contribute $D_0h^{k-1}\sum_{j=1}^k\alpha_j$. The bounded integer
intercepts give only $O(k+1)$ such sums. So there are $O((k+1)B^k)$
counts at depth $k$. If $C$ is positive initially or ever becomes
positive then $N_k$ grows on the scale $P^k$: use $\beta\ge P/2$,
$C'\ge PC$, and the nonnegativity of the remaining terms. Therefore
the number of distinct values at most $X$ is
$$
O(X^{\log_P B}\log X)\quad(B\ge2),\qquad
O((\log X)^2)\quad(B=1).
$$
If the address stays zero, all used marks are zero, and there is just
one baseline per depth. In that case changing patches in the same
palette has no effect on the scalar count.

The same proof applies to a predetermined schedule from a fixed
finite catalogue of matched palettes. For a type word with multipliers
$h_j$, the intercept contribution divided by $\prod_j h_j$ is
$D_0\sum_j\alpha_j/h_j$. A fixed common denominator again gives
$O(k+1)$ choices. After activation of the address the growth scale
is $\prod_j P_j$. The maximal $\log_{P_j} B_j$ over the catalogue
is strictly below 1. This case does not assume every type has the
same growth rate.

## A weighted bound when the branches choose different types

Let type $i$ denote a matched palette with $P_i$ its slot multiplier
and $B_i$ its number of slopes. A type change may change the whole
continuation map. Allow the branch to choose the type at each step,
subject to physical compatibility of boundary marks. If there is
an $s>0$ with
$$
\sum_i B_i P_i^{-s}<1,
$$
the type-dependent construction still has only $O(X^s)$ distinct
counts through $X$. More generally it suffices that the weighted
path sum for type words, with the extra factor $k+1$ at length $k$,
converges. For example, a finite automaton of admissible type words
can be handled by weighting each edge for type $i$ by $B_iP_i^{-s}$
and bounding the spectral radius on states that lie on admissible
paths by a number strictly below one.

To prove the sum criterion, fix an activated prefix and a type word
in the subsequent suffix. Its at most $O((k+1)\prod_j B_j)$ counts
each exceed a constant times $\prod_j P_j$. The Dirichlet sum of
these bounds over suffix words is finite under the displayed
condition, including the polynomial factor. If the initial address
is zero, any preactivation word consists only of unmarked patches
and scales both $D$ and $N$ by a product of $h_i$'s. The order of
these factors can be ignored. Its inverse-count mass is bounded
by a finite product of geometric sums for those $h_i>1$, with
$h_i=1$ inert. The first activation is a choice from a finite set.
Composing the prefix with any activated suffix preserves the bounds
by homogeneity of the linear transfer. This treatment counts
distinct scalars, not all redundant representations with zero padding.

## Concrete adaptive obstruction: arbitrary unmarked SCCs through six

Allow the previously published patches A and B to be chosen in any
order, and let the branch also insert arbitrarily many unmarked
ordinal blocks anywhere, subject only to every strong component
of an inserted block having order at most six. Each insertion
comes before, after, or between intact copies of A and B, not inside
their specified physical templates. The upper bound is
$$
O(X^{21/40}) .
$$
For this example, do not charge each unmarked insertion as a
separate arbitrary type. Unmarked transfers are exactly
$$
(D,C,N)\mapsto (hD,PC,hN+2(P-h)C).
$$
They commute, and their product is determined by the two products
$h,P$. Hence between consecutive A/B decisions an unmarked run
can be counted by multiplicities of the strong profiles
from the [literal component replay](../chanpak_guard_depth_20261115/replay_menu.txt),
including the singleton $(1,2)$. Use the profile bank by order
(with harmless repetitions of the same $P$). For each such gap
the inverse-count weight is bounded by
$$
U(s)=\prod_i\frac1{1-P_i^{-s}} .
$$
Every noninitial A/B decision contributes at most $2\cdot1024^{-s}$,
and an arbitrary terminal unmarked gap contributes one factor
$U(s)$. The first decision activates a zero address and can be
paid as a fixed constant; the unmarked initial prefix is handled
by its products of $h_i$ as above. The geometric sum converges
if $2\cdot1024^{-s}U(s)<1$. The [rational certificate](verify.py)
checks this inequality at $s=21/40$, without treating approximate
floating-point values as proof.

This does not exclude arbitrary adaptive palettes. It distinguishes
two real design changes: varying the coarse carrier enough to
exceed the weighted slope budget, or breaking the uniform ordinal
seam which imposes the slope ceiling in the first place. It is not
an obstruction on the full Chan–Pak spectrum.
