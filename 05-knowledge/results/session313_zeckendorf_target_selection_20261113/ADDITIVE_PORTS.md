# Independent physical switches without a transitive spine

This is an actual family of tournaments, not just a representation of the target integer. It starts from the [one-apex cut formula](../session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md) and its failure as a separable digit bank.

For a single apex over a transitive spine with bit word $b$, the inherited formula is $H-1=\sum_{i<j} b_i(1-b_j)c_{j-i}$, with $c_1=2$ and $c_h=2^{h-1}$ for $h\ge2$. If two disjoint equal-width edits are separated by at least one unchanged vertex, their mixed difference factors as
$$
-\left(\sum_{i\in A}\delta_i2^{-i}\right)\left(\sum_{j\in B}\epsilon_j2^{j-1}\right).
$$
Each factor is nonzero for two different binary blocks, by uniqueness of binary positional value. Thus even two such fixed-position genuine edits cannot be independent. The trick below replaces the transitive sites by modules with identical responses in just the directions where they interact.

## The one-sided response law

Take an ordinal chain of modules $T_1,\ldots,T_n$ and an apex $x$. Within $T_i$, let $U_i$ be the out-neighbors of $x$ (the permissible starts on the right of $x$). Write $h_i$ for the module's Hamiltonian count and $q_i$ for its ordered two-slot path-cover count allowing empty slots. For a subset $A$ in a module, let $S(A)$ count nonempty Hamiltonian paths of $A$ starting in $U_i$ and $E(A)$ those ending in its complement; both are zero for the empty set, while $H(\varnothing)=1$. Let
$$
s_i=S(T_i),\quad e_i=E(T_i),\quad
s_i^*=\sum_A S(A)H(T_i\setminus A),\quad
e_i^*=\sum_A E(A)H(T_i\setminus A),\quad
b_i=\sum_A S(A)E(T_i\setminus A).
$$

Split a full path at $x$. Each side visits modules in increasing order. If both sides are nonempty, write $i$ for the last module on the left and $j$ for the first on the right.
- If $i<j$, necessarily $j=i+1$; the contribution is $e_i s_j$ times the product of $h$ at all other sites.
- If $i=j$, the contribution is $b_i$ times the product of $h$ at all other sites.
- If $i>j$, the two end modules supply $e_i^*s_j^*$, the strict interval $(j,i)$ supplies a product of $q$, and modules outside $[j,i]$ supply a product of $h$.
- If the left or right side is empty, use $s_1$ or $e_n$ and the product of $h$ at all other sites.

These cases are disjoint and exhaust the paths. Therefore any switches that fix $(h,q,s,s^*)$ at every site can vary $(e,e^*,b)$ **independently and additively** across arbitrary numbers of sites. It is enough to fix the incoming factors of the products, not the entire boundary.

## Two literal switch types

Using zero-based outgoing masks, $(20,1,2,7,14)$ has $(h,q)=(9,84)$. The two apex cuts $U=\{0\}$ and $\{2\}$ both have $(s,s^*)=(1,11)$, and $(e,e^*,b)$ is respectively $(8,64,10)$ and $(6,52,10)$. In a chain of $n$ identical modules, switching from $\{2\}$ to $\{0\}$ at site $i$ (now one-indexed) increases $H$ by
$$
2\,9^{n-1}+132\sum_{j=1}^{i-1}84^{i-j-1}9^{n-i+j-1}.
$$
All subsets of switches are legal. The sum is absent at the first site. The displayed rule is affine independently in every bit, not just pairwise along a particular path in the space of cuts.

The other module $(24,17,3,6,12)$ has $(h,q)=(15,120)$. Cuts with masks 5 and 3 have the same four external responses $(s,s^*,e,e^*)=(6,42,9,63)$, but $b=28$ and 26. They yield a purely internal $+2$ times the product of all other $h$, in **any** ordinal chain of modules sharing an apex. In particular the change is exactly $+2$ whenever the other modules are transitive, regardless of their sizes or apex incidences. This is an exact offset with the four named responses unchanged. It does not assert that the second output is ready for every possible nonuniform continuation.

## Independence is not enough for coverage

For the repeated order-five binary bank, every count is at least $9^n$: every path of the ordinal chain admits an insertion of the apex, and deleting it recovers the original path. Thus through target $X$, all bank lengths together supply at most $\sum_{n\le\log_9 X}2^n=O(X^{\log_9 2})$ counts. This is a zero-density family. The switches are truly independent, but they do not generate offsets at sufficient density or anchor a prescribed target. The pure $+2$ switch does not by itself solve that problem either: it needs a family whose unshifted values hit the adjacent targets.

The next quantitative design test is to find interacting port responses that give enough choices per multiplicative cost, and verify overlap on the actual output ranges. The symbolic Fibonacci/Jacobsthal decoder tells us what such a range proof could look like; it is not a substitute for the physical anchors.
