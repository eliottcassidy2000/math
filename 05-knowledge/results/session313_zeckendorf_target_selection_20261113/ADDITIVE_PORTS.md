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
2\,9^{n-2}+132\sum_{j=1}^{i-1}84^{i-j-1}9^{n-i+j-1}.
$$
For the terminal site use $2\,9^{n-1}$ instead of $2\,9^{n-2}$. All subsets of switches are legal. The sum is absent at the first site. The displayed rule is affine independently in every bit, not just pairwise along a particular path in the space of cuts.

The other module $(24,17,3,6,12)$ has $(h,q)=(15,120)$. Cuts with masks 5 and 3 have the same four external responses $(s,s^*,e,e^*)=(6,42,9,63)$, but $b=28$ and 26. They yield a purely internal $+2$ times the product of all other $h$, in **any** ordinal chain of modules sharing an apex. In particular the change is exactly $+2$ whenever the other modules are transitive, regardless of their sizes or apex incidences. This is an exact offset with the four named responses unchanged. It does not assert that the second output is ready for every possible nonuniform continuation.

## Independence is not enough for coverage

For the repeated order-five binary bank, every count is at least $9^n$: every path of the ordinal chain admits an insertion of the apex, and deleting it recovers the original path. Thus through target $X$, all bank lengths together supply at most $\sum_{n\le\log_9 X}2^n=O(X^{\log_9 2})$ counts. This is a zero-density family. The switches are truly independent, but they do not generate offsets at sufficient density or anchor a prescribed target. The pure $+2$ switch does not by itself solve that problem either: it needs a family whose unshifted values hit the adjacent targets.

The next quantitative design test is to find interacting port responses that give enough choices per multiplicative cost, and verify overlap on the actual output ranges. The symbolic Fibonacci/Jacobsthal decoder tells us what such a range proof could look like; it is not a substitute for the physical anchors.

## A physical target decoder in a sparse 19-adic family

The order-six core with outgoing masks $(40,1,3,6,15,30)$ (also used in Session 309) has $(h,q)=(19,256)$. The out-cut masks 40 and 17 have the same $(s,s^*)=(8,94)$, with $(e,e^*,b)$ equal to $(15,175,72)$ and $(17,187,78)$. Switch from 40 to 17. At site $i$ the resulting increment in the repeated core bank is
$$
w_i=130\,19^{n-2}+1128\sum_{j=1}^{i-1}256^{i-j-1}19^{n-i+j-1}
\quad(i<n),
$$
and at $i=n$ use $8\,19^{n-1}$ instead of $130\,19^{n-2}$. The nonterminal contribution is $6\,19^{n-1}+2\cdot8\,19^{n-2}$ (the same-module split uses one block, an adjacent split uses two); the terminal contribution is $(6+2)19^{n-1}$. For $i\ge2$ the $j=1$ summand has valuation $n-i$ and coefficient $1128\cdot256^{i-2}$ a unit mod 19. For $i=2<n$, the constant term has the same valuation, but the sum of the two coefficients is $1258=4\bmod19$, still a unit. Thus for $n\ge2$, the valuations of sites $2,\ldots,n$ are $n-2,\ldots,0$ respectively (at $n=2$ the terminal coefficient is $1280=7\bmod19$). Site 1 also has valuation $n-2$.

Fix site 1. Subtract the all-40 baseline from a target. Reduce mod 19 to identify the final switch bit: the residue must be 0 or the unit of $w_n$. Subtract its contribution, divide by 19, and continue backward through site 2. Require the final integer remainder to be zero. This is a genuine exact selector for this physical bank. The latest differing site determines the valuation of two output counts' difference. Thus the bank realizes a binary prefix tree of target addresses in 19-adic distance, but only $2^{n-1}$ of $19^{n-1}$ residue classes. The terminal remainder test distinguishes exact targeting from a residue match. Site 1 can also be tried both ways; no claim of all odd counts follows.

## Exact independence test, including the adjacent seam

The one-sided guard is sufficient but not the only way to obtain independence. In fact the same formula gives a complete test when each site's two choices keep $h_i,q_i$ fixed. Let $\Delta$ denote the change from option zero to option one. For $j<i-1$, the mixed coefficient for the two sites $j,i$ is
$$
\Delta s_j^*\,\Delta e_i^*
\prod_{j<k<i}q_k\prod_{k<j\ \mathrm{or}\ k>i}h_k .
$$
For $i=j+1$ it is instead
$$
(\Delta s_j^*\,\Delta e_i^*+\Delta e_j\,\Delta s_i)
\prod_{k<j\ \mathrm{or}\ k>i}h_k .
$$
There are no higher-order terms, so the switches are independent exactly when all these coefficients vanish. Every product outside the parentheses is positive. In particular, if there is a first site changing $s^*$ and a last site changing $e^*$, the latter cannot lie more than one step to the right of the former. Before the first site all two-slot starts are fixed; after the last site all two-slot ends are fixed. Only adjacent sites can cancel a nonzero two-slot interaction against a one-slot interaction. This characterizes the obstruction without assuming the one-way guard up front.
