# Exact open-to-closed line laws

This note uses the **literal merged complement-line multigraph**, not the single-arc-flip graph or its fundamental paired-black quotient. One line is a complementary pair of labelled tilings; count it once even if its class endpoints coincide. The open/closed convention is inherited from [THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md); the score flux and Hodge viewpoint are inherited from [THM-4739](../../../01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md). Here the flux becomes an exact all-order generating function, including its colors.

Take the path $1\to\cdots\to n$ and orient each free edge $i<j$ by sign $x_{ij}=+1$ for forward. In the open representative $x_{1n}=+1$. Put $h=e_1-e_n$, and let
$$
z=\sum_{e=(i,j)\in D} x_{ij}(e_i-e_j),
$$
where $D$ is the edge set of $K_n$ with the whole undirected Hamilton cycle removed. The signed score vectors of the open and closed ends are $2h+z$ and $-z$. This isolates the endpoint correction from the divergence of the remaining free graph. For any signed score vector $u=2s-(n-1)\mathbf1$, the number of cyclic triangles is
$$
f(u)=\frac{n(n^2-1)}{24}-\frac{\|u\|^2}{8}.
$$
Thus the directed flux $\Delta=c_3(\mathrm{closed})-c_3(\mathrm{open})$ is $1+(z_1-z_n)/2$.

Let $m=\binom{n-1}{2}$ and $b=(m+\lfloor(n-1)/2\rfloor)/2$ be the total and mirror-fixed sign-cube dimensions. For every $n\ge3$ the Laurent polynomials counting all and blue lines by flux are respectively
$$
L_n(t)=2^{m-1-2(n-3)}\,t^{4-n}(1+t)^{2(n-3)},\qquad
B_n(t)=2^{b-1-(n-3)}\,t^{4-n}(1+t^2)^{n-3}.
$$
The black polynomial is $L_n-B_n$. For the first formula, exactly $2(n-3)$ freely varying signs enter $z_1-z_n$, each with coefficient one. The other signs supply the power of two. For blue lines those endpoint signs pair under mirror reflection, giving $n-3$ free signs of coefficient two. In particular every blue line has $\Delta\equiv n\bmod2$; zero-flux blue lines are impossible at odd orders. Both colors have mean directed flux 1 when nonempty (the black case follows by subtraction). Normalizing all lines gives variance $(n-3)/2$; normalizing blue gives variance $n-3$. This counts line occurrences, not the simple support, and it does not infer that any given score-compatible class pair is adjacent.

## Joint score-to-triangle generating function

Use Laurent monomials $X^v=\prod_i X_i^{v_i}$. For all lines take the product over singleton edge sets in $D$; for blue take the product over mirror orbits in $D$. In either case write $v_O=\sum_{(i,j)\in O}(e_i-e_j)$ and form
$$
P(X)=\prod_O(X^{v_O}+X^{-v_O}).
$$
The coefficient $[X^z]P$ counts the exact signed divergence $z$. Map each monomial $X^z$ to $a^{f(2h+z)}c^{f(z)}$ to obtain the full two-variable distribution by open and closed triangle levels. Multiplying one factor at a time is an exact vector recurrence; retaining the order of the coordinate labels is essential until both projections are applied. Difference of the two products after projection gives black. This is stronger than the separate vertex degree/level inventories, but weaker than class-pair incidence.

For all lines, $c_3(\mathrm{closed})$ has mean $n(n^2-1)/24-(m-1)/4$ and variance $n\binom{n-3}{2}/16$. In the squared norm of $z$, the centered terms are pairs of distinct edge signs with coefficient twice the inner product of the two incidence vectors. The inner product has absolute value one exactly for edges sharing a vertex. The graph $D$ has degree $n-3$ at every vertex, and distinct pair monomials are orthogonal. The linear flux is orthogonal to all these quadratic terms, so the covariance of the open and closed triangle levels is exactly the closed variance. For blue lines the same argument gives closed variance $\sum_{O<O'}(v_O\cdot v_{O'})^2/16$; the flux remains uncorrelated with the closed level.

Direct sign-cube controls through $n=7$ verified both binomial laws and the closed variance. The proof uses all-order sign factors, not extrapolation.
