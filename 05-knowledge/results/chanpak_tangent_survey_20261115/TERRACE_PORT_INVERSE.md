# Exact inverse for two exposed terrace gaps

The older [square-sum port compiler](../square_sum_spectrum_overnight_20261031g/PORT-ORDERS.md) retains common labelled endpoints and shared bit constraints. The [prescribed-edge absorbers](../square_sum_prescribed_edge_absorbers_20261019.md) similarly separate a local replacement from its eligibility in a prescribed host. Here that distinction guides actual tournament count selection: expose two terminal widths but retain the same entire interior word. This produces a complete inverse inside a physical family, not a conclusion from two separately feasible boundaries.

## Physical model and transfer

Use the backward-chain tournaments of [the first triangle terrace](../session337_rooted_line_surfaces_20261114/TERRACE.md). Start with a forward tournament, choose endpoint-inclusive anchors, and reverse just the edges from each later anchor to the previous anchor. A gap is the number of nonanchor sites between successive anchors. The cited theorem identifies the Hamiltonian count with a selection sum on gaps: a maximal selected run containing a positive gap contributes the product of its positive weights $2^a$; a pure-zero selected run has weight 2 at even length and 0 at odd length. Our words have positive first and final gaps, and these are strong tournaments.

Let the positive gaps be $a,c_1,\dots,c_j,b$ with intervening zero corridors $q_0,\dots,q_j$ (each may be zero). Fix all $c_i,q_i$; only $a,b$ vary. Put
$$
K_q=\begin{pmatrix}A_q&B_q\\ B_q&D_q\end{pmatrix},
$$
where $K_q(s,t)$ is the zero-corridor selection sum given the neighboring positive gap selection bits $s,t$. The [inherited corridor formulas](../chanpak_target_first_terrace_density_20261114/TERRACE_DENSITY.md) give the coefficient generating functions $(1+z^2)/\Delta,(1+z)/\Delta,(1+z+z^2)/\Delta$, respectively, for $\Delta=1-z-z^2-z^3$. The complete interior profile is
$$
M=K_{q_0}\operatorname{diag}(1,2^{c_1})K_{q_1}\cdots
  \operatorname{diag}(1,2^{c_j})K_{q_j}
 =\begin{pmatrix}A&C\\B&D\end{pmatrix}.
$$
Condition on the two terminal selection bits to get the literal path count
$$
F(x,y)=(1,x)M(1,y)^{\mathsf T}=A+Bx+Cy+Dxy,\qquad x=2^a,\quad y=2^b.
$$
This needs no assumptions of low rank, symmetric context, or factoring of the whole count. It applies when no corridor has length 0 or 13.

## Complete prescribed-count decoder

Each $K_q$ is entrywise positive and increasing in each of its two bits. Indeed, the differences have generating functions $z(1-z)/\Delta$ and $z^2/\Delta$, with nonnegative coefficients. Positive linear combinations preserve the two end monotonicities; thus $1\le A,B,C\le D$.

Unless $a=b=1$,
$$
2^{a+b}\le F(2^a,2^b)/D
\le 2^{a+b}(1+2^{-a})(1+2^{-b})
\le (15/8)2^{a+b}<2^{a+b+1}.
$$
The count at $(a,b)=(1,1)$ is a strict minimum by positivity. If a prescribed $N$ equals this minimum it is the only answer; if $N$ is smaller there is none. Otherwise any answer must have
$$
k=a+b=\lfloor\log_2\lfloor N/D\rfloor\rfloor,\qquad
E=N-A-D2^k,\qquad
Bx^2-Ex+C2^k=0.
$$
Compute the discriminant $E^2-4BC2^k$, require a nonnegative integer square, and test the at most two integral roots $x$ for being powers of two with $2\le x\le 2^{k-1}$. The corresponding $y=2^k/x$ is forced. A failed test proves absence from the chosen interior family, not absence from the tournament spectrum. A success directly constructs the ordered tournament from the two exponents. This is exact inversion with no factorization and no enumeration of terminal gap lengths.

**Only possible alias: swapping the terminals.** Each corridor has adjacent row and column ratios in $[1,2]$, and its ratio $D_q/B_q$ is strictly below 2, because
$$
2A(z)-B(z)=(1-z+2z^2)/\Delta,\qquad
2B(z)-D(z)=(1+z-z^2)/\Delta
$$
have respectively nonnegative and strictly positive coefficients. Use the nondecreasing coefficients of $1/\Delta$ in each difference. Taking positive weighted averages in the first and last kernels shows $1\le B/A<2$ and $1\le C/A<2$, so $1/2<B/C<2$. If two distinct solutions exist at the same $k$, the partner root must be $x'=Cy/B$; it is a power of two only if $B/C$ is a dyadic power, hence only if $B=C$. Then $x'=y$ and $y'=x$, exactly the swap. For $B\ne C$ an ordered answer is unique for every target.

With one fixed interior positive gap of weight $w=2^c$ and corridors 1 and 2, the matrix is
$$
\begin{pmatrix}3+6w&3+8w\\6+6w&6+8w\end{pmatrix}.
$$
Here $3/4<B/C<1$ for every $c\ge1$, giving ordered uniqueness at every fixed middle width. When there is no interior positive, the result extends the previously worked single-zero shifted-product inverse to every zero corridor. One can scan all candidate corridors without factoring: $D_q\ge A_q\ge3^{\lfloor q/3\rfloor}$ by allowing disjoint selected even pairs separated by a skipped site. Since $N\ge4D_q$, only $O(\log N)$ lengths can work. For the explicit asymmetric middle with corridors fixed at $(1,2)$, there are likewise only $O(\log N)$ middle widths to test by the $D$ coefficient.

This removes terminal address search but does not supply an interior for each target. A fixed finite palette of interiors gives only $O((\log X)^2)$ distinct counts through $X$ since $a+b\le\log_2 X$. The inherited theorem that the entire first terrace is sparse still applies. The point is an exact part of an inverse obtainable without a heuristic congruence lift; target-driven supply of the interior must come from a richer construction.

Run the attached Python replay for 334 direct subset-DP controls, as well as exhaustive small-window comparisons of the complete address fibers.

