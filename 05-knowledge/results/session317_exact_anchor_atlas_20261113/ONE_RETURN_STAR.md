# One-return-star surgery: exact separated path-count coordinates

**Proved for the cut-one class.** This is an inverse-design reduction, not a proof of Chan–Pak or of an eventual threshold.

## Historical transfer

A reproducible shuffle of top-level results Markdown using repeating bytes `317\n` sampled [golden difference diagonals](../difference_diagonals_golden_phase_recursion_20261004.md), the adjacent [denominator realization gate](../denominator_arithmetic_filter_20261004.md), [composite Pell and prime-power conductor layers](../tree_resonance_composite_pell_20261013b.md), [Gold prefix attachments](../gold_prefix_side_poset_20261019.md), [the next day's debt cones and idoneal dictionary](../debtcones_20261020_endpoint_layers_zero_heavy_cones_annulus_idoneal.md), [the four-node trijection](../trijection_four_nodes_three_edges_20261020.md), and [the mixed two-jet audit](../overnight_20260906_moments_smith_audit.md). Follow-up reads included [Collatz radix transport](../collatz_radix_families_20261018p.md) and [fair-attachment path-cover transforms](../session26_random_bridges_20261022/session26-attachments/proof.md).

The Gold construction separates free internal shuffles from the guard at the first touched boundary. The golden phase atlas still needs an integer-realization predicate after finding a formal itinerary. The radix note retains both ends of a word because the front controls size and the back controls the next branch. These transfer here as a signed exposed-boundary correction, physical admissibility of affine charts, and an exact digit-reversal coordinate.

Sessions 303 and 315 give the maximum visit degree of an inflated root as one plus the matching number of the bipartite graph of arcs from its out-neighbors to its in-neighbors. Every interior excursion consumes a crossing arc. Conversely matched arcs furnish direct excursions, with unused in-neighbors and out-neighbors traversed as prefix and suffix. This matching observation is inherited.

## Exact transfer

A nonempty bipartite graph of matching number one is a star. Taking converse if necessary, every such root chart has the following form: start with the ordinal join $I\to U\to W$, where $U$ is an arbitrary nonempty tournament replacing the root; then add apex $v$, dominated by all of $U$, beating a nonempty set $S\subseteq I$, with unrestricted orientations to $W$. Conversely these charts have matching number one. If the star is a single edge, either endpoint can serve as its centre after choosing orientation.

Write $h_X=H(X)$, using $h_\varnothing=1$ for empty side blocks. Let $K(X)$ count ordered two-slot directed path covers allowing empty slots. Define:

- $\alpha$: ordered covers of $I$ with nonempty suffix starting in $S$, prefix optional.
- $\lambda$: those counted by $\alpha$ whose prefix is nonempty and ends in $S$.
- $\beta$: ordered covers of $W$ with prefix either empty or ending at a vertex beating $v$; suffix unrestricted.
- $d=H(W+v)$, on the actual induced tournament.

Then
$$
H=\alpha\beta K(U)+(h_I d-h_W\lambda)H(U).
$$

**Proof.** Delete $v$ from a Hamiltonian path and keep the ordered pieces before and after it. Each piece meets the blocks in ordinal order. There are three cases.

1. Both pieces meet $U$: the post-$v$ piece must begin in $I$, starting in $S$; the pre-$v$ piece may end in $U$, or in $W$ at a vertex beating $v$. Its contribution is $\alpha\beta(K(U)-2H(U))$.
2. $U$ entirely before $v$: a nonempty $I$-suffix gives $\alpha\beta H(U)$ as before. If the $I$-suffix is empty, the $W$-pieces together with $v$ form any Hamiltonian path of $W+v$, giving $h_IdH(U)$.
3. $U$ entirely after $v$: a nonempty $I$-suffix is again necessary. If the $W$-prefix is nonempty, the two independent guards suffice. If it is empty, the $I$-prefix must be empty or end outside $S$; subtract the illegal $h_W\lambda$ splits. The contribution is $(\alpha\beta-h_W\lambda)H(U)$.

Adding proves the formula. The affine chart for $U=TT_m$ has slope $B=\alpha\beta>0$ and intercept $C=h_Id-h_W\lambda$. Here $\beta\ge h_W$; to see $\alpha>0$, pick a singleton suffix in $S$ and a Hamiltonian prefix on its complement.

## Every slope, but a coupled offset

Take $I=TT_n$, $W=\varnothing$. Let the ordinary binary membership mask of $S$ be $B$, with the first vertex supplying its high bit. A suffix starting at position $j$ gives $2^{n-j}$ choices, so $\alpha=B$. Write $L(B)=\lambda$, $L(0)=0$. Prepending zeros does not change $L$: each such new vertex must be in the prefix and cannot be its terminal vertex in a counted split. Appending one bit yields
$$
L(2B)=L(B),\qquad L(2B+1)=L(B)+B+(B\bmod2).
$$
A terminal zero must be in the suffix. A terminal one can be in the prefix and paired with an older selected suffix start, accounting for $B$, or be the sole suffix, contributing one iff the old terminal bit is selected. Hence
$$
H(B,m)=2^mB+1-L(B),\qquad m\ge1.
$$
In particular $L$ is even. The absolute spectrum with arbitrary $m\ge1$ is already the $m=1$ spectrum by replacing $B$ with $2^{m-1}B$. Indeed every nontransitive one-apex tournament reduces to this form by choosing the last vertex beating the apex as $U$ and dropping the sink tail. The one-apex spectrum's absence of 13 and its positive-density omissions despite modular completeness are inherited from [Session 26](../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md). The present claim concerns the rooted slope–offset pair.

For fixed intercept $1-D$, the admissible slopes form a finite union of dyadic rays: write $B=2^ju$ with $u$ odd, then $D=L(B)=L(u)\ge(u-1)/2$. Equality holds only for $u=1$ or $u=2^a+1$, $a\ge2$. Indeed for $u=2q+1$ equality requires $q$ even and $L(q)=0$; the recurrence says $L(q)=0$ iff $q=0$ or $q$ is a power of two. Thus the rays are precisely $2^ju$ with odd $u\le2D+1$ and $L(u)=D$.

Let $b_j$ be the bits from the right and $V(B)=\sum_j b_j2^{-j-1}$ the binary radical inverse. Unwinding the recurrence gives
$$
L(B)=\sum_j b_j(\lfloor B/2^{j+1}\rfloor+b_{j+1}),\qquad
|L(B)-BV(B)|\le1+\lfloor\log_2B\rfloor.
$$
For uniform $2^r\le B<2^{r+1}$, $V(B)$ is a shifted mesh of step $2^{-r}$. Thus $L(B)/B$ converges in distribution to uniform on $[0,1]$, with CDF error at most $(r+2)/2^r$. The continuous limiting range does not certify exact targets.

**Exact radix seam.** Write $B=2^rq+\ell$, $r\ge1$, $0\le\ell<2^r$, and let $R_r(\ell)$ reverse exactly those $r$ bits. Then
$$
L(B)=L(q)+qR_r(\ell)+L(\ell)+(q\bmod2)\lfloor\ell/2^{r-1}\rfloor.
$$
The low-bit floor terms give $qR_r(\ell)$; the sole adjacent-bit term across the seam is the final product. Consequently
$$
H(B,m)=(2^{m+r}-R_r(\ell))q+1+2^m\ell-L(\ell)-L(q)
 -(q\bmod2)\lfloor\ell/2^{r-1}\rfloor .
$$
The inverse needs both the residual $L(q)$ and its parity guard.

## Cancellation and a two-sided exact construction

If $\lambda=0$ and the intercept is one, $h_Id=1$ forces $I,W+v$ transitive. A cyclic triple provides three distinct paths that can be extended by insertion preserving their orders, ruling out count one in the nontransitive case. Equation $L(B)=0$ then forces $\alpha$ to be a power of two. In $W+v$, the vertices beating $v$ form a prefix of length $j$, and its subsets give exactly $\beta=2^j$. Thus the cancellation-free intercept-one subclass yields only $2^k+1$, $k\ge1$, and all those counts occur. More generally a prime positive intercept with $\lambda=0$ forces one of $I,W+v$ transitive.

With cancellation, for $k\ge1$ take $I=TT_{k+2}$ and select its first and last vertices. Then $\alpha=2^{k+1}+1$, $\lambda=2^k$. Take $W=TT_{k+1}$ with $v$ beating the first $k$ vertices and losing to its last vertex. The last is the only allowed nonempty prefix endpoint, giving $\beta=2^k+1$. Also $W+v$ is the cyclic join of the first $k$-chain and the two singleton blocks, with $d=2^k+1$ (three one-run itineraries, one two-run itinerary). Hence the exact intercept is $d-\lambda=1$, and
$$
H=(2^{k+1}+1)(2^k+1)2^m+1,\qquad k,m\ge1.
$$
The slopes are odd and grow unboundedly, escaping one-sided fixed-offset rigidity. A smaller witness uses $I=TT_2$, $S=I$, and $W+v=C_3$, giving $H=9\cdot2^m+1$.

## Remaining inverse problem

For transitive $I$ one can select any $\alpha>0$ by its binary digits. Prescribing intercept $c$ then asks for $H(W+v)-H(W)L(\alpha)=c$, but the prescribed slope also needs $\beta$. It is not enough to solve the intercept equation abstractly or modulo arbitrarily many primes. A changing-core construction would need a return that refreshes these guarded two-slot counts. Nothing here proves that return or an absolute coverage threshold.

Replay: `python3 replay_one_return.py`. It independently evaluates literal expanded tournaments by subset DP, evaluates slot counts by permutations, and tests the seam in a separate finite arithmetic universe.
