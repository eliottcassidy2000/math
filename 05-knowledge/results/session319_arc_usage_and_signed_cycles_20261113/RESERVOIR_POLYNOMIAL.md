# One neutral arc creates irreversible cover readiness

This extends the equal-count dyadic pair from [Session 318](../session318_reservoir_descent_20261113/DYADIC_SWITCH.md) to arbitrary internal tournaments, not just a transitive run. It also gives a sharp obstruction to using the switch for a homogeneous exact return.

Label the seven quotient vertices `a,b,x,t,c,d,e`. The two outgoing-mask lists are
`A=(54,64,122,115,98,66,1)` and `B=(54,32,122,115,98,64,3)`.
Replace $t$ by any common nonempty tournament $U$, with uniform inherited orientations. Write $A[U],B[U]$ for the two results, $n=|U|$, $H(U)$ for Hamilton paths, $F_j(U)$ for ordered covers by $j$ nonempty paths, and $Q_U(z)$ for ordered covers with empty slots allowed, evaluated at integer $z$ and extended polynomially.

## Literal neutral flip and preserved arc marginals

First relabel $A$ by swapping $b,d$, fixing $U$ pointwise. This reverses only $d\to b$ in the quotient: both vertices have identical other incidences. Now reverse $b\to e$ to $e\to b$; the result is $B$. This is a single arc flip up to the indicated relabelling.

The only vertex outside $U$ that beats any of $U$ is $x$. In a Hamilton path there can be at most two runs of $U$, and a two-run path must start in $U$. The quotient-word counts for both $A,B$ are 27 with one run and 9 with two. In the one-run class there are 9 words starting at $t$, 1 ending at $t$, and 18 using $xt$. In the two-run class all 9 start at $t$, exactly 1 ends there, and all 9 enter their second $t$-run from $x$. The replay enumerates these small quotient-word sets. Therefore
$$
H(A[U])=H(B[U])=27H(U)+9F_2(U)=9(H(U)+Q_U(2)).
$$
Both outputs are strong. Indeed there is exactly one cyclic quotient word rooted at $t$ through the other six vertices, `tcdbeax` or `tcbdeax`, respectively. Every Hamilton cycle through $U$ has one run, since each run needs its own incoming contact $x$. Thus each graph has exactly $H(U)$ directed Hamilton cycles (modulo cyclic shift).

Moreover every named vertex of $U$ has the same Hamilton start and end counts in both outputs; every internal arc of $U$ and every entry arc $x\to u$ has the same number of paths using it. This follows by filling the 27 one-run and 9 two-run words with the same internal paths or covers and using the incidence table just stated. When $U=TT_L$ with increasing order $t_1,\dots,t_L$, these common marginals are explicit:
$$
\begin{aligned}
\mathrm{start}(t_j)&=9\cdot2^{L-j},&
\mathrm{end}(t_j)&=2^{j-1},\\
w(x,t_j)&=9\cdot2^{L-j}+9[j=1],&
w(t_i,t_j)&=9\cdot2^{L-(j-i)}+9[j=i+1]\quad(i<j).
\end{aligned}
$$
For a fixed first vertex in a given run, earlier labels are forced into the other run and later labels are free, with a correction if that other run can be empty; ending is dual. For an internal edge $t_i\to t_j$, assign both endpoints to a chosen run, all strictly intervening labels to the other, and all remaining labels freely. The two-run count is $2^{L-(j-i)}-2[j=i+1]$; adding the 27 one-run copies of consecutive edges gives the formula.

## Exact polynomial dominance in every ordered-cover coordinate

Put $R_U(z)=Q_U(z)/z$. For all $U\ne\varnothing$,
$$
\boxed{Q_{B[U]}(z)-Q_{A[U]}(z)
=z^2(z^2-1)[R_U(z+1)-R_U(z-1)].}
$$
To prove it, use the odd-cycle packing expansion (each selected odd directed cycle has weight $2z$; an uncovered vertex contributes $z$). A mixed directed cycle has exactly one run in $U$ and uses $x$, so at most one such cycle occurs in a packing. A quotient cycle with $d$ outside vertices, rooted at $t$, permits any nonempty directed path on a subset $S\subseteq U$ as its run. Its actual order is $d+|S|$. After deleting it, the remaining block $U-S$ cannot participate in a mixed cycle, and independently contributes $Q_{U-S}$.

Let $Q_0$ be the packing polynomial on the six outside vertices, and $G_0,G_1$ the sums over quotient cycles rooted at $t$ of the outside-remnant packing polynomials, grouped by even/odd $d$. The quotient table is

| | $Q_0$ | $G_0$ | $G_1$ |
|---|---|---|---|
| $A$ | $z^6+8z^4+8z^2$ | $z^4+3z^2+1$ | $z^3+3z$ |
| $B$ | $z^6+8z^4+6z^2$ | $z^4+4z^2+1$ | $z^3+2z$ |

Slot convolution $Q_U(z+y)=\sum_{S\subseteq U}Q_{U[S]}(y)Q_{U-S}(z)$, and $Q_{U[S]}(-1)=(-1)^{|S|}H(U[S])$, give
$$
Q_{A[U]}=Q_0Q_U(z)+zG_0[Q_U(z+1)-Q_U(z-1)]
 +zG_1[Q_U(z+1)+Q_U(z-1)-2Q_U(z)]
$$
and the analogous expression for $B$. Subtracting the table cancels $Q_U(z)$ and gives the claimed factorization. The replay verifies the quotient table by explicit directed cycles and packings, independently of literal subset DP on expanded instances.

For $n\ge2$, $R_U$ has nonnegative coefficients and leading term $z^{n-1}$, by the same packing expansion. Thus the difference polynomial can be written $z g(z-1)$, with $g(u)=(u+1)u(u+2)[R_U(u+2)-R_U(u)]$, a polynomial with nonnegative coefficients, zero constant, and positive leading coefficient of degree $n+1$. Inclusion-exclusion for all slots nonempty gives
$$
F_j(B[U])-F_j(A[U])
=\sum_{i=0}^j(-1)^{j-i}\binom ji i\,g(i-1)
=j\,\Delta^{j-1}g(0).
$$
For a monomial $u^d$, $\Delta^k u^d(0)$ counts onto maps from $d$ objects to $k$ labelled boxes, hence is positive for $1\le k\le d$ and zero for $k>d$. The leading term of $g$ therefore proves strict dominance at every $2\le j\le n+2$, with equality at $j=1$ and $j>n+2$. For $n=1$ the entire polynomial difference vanishes.

**No common homogeneous return.** For any common substitution context, each ordered output cover projects to an ordered list of quotient walks, allowing repetitions across and within slots. Filling the marked block’s ordered runs contributes its corresponding ordered nonempty cover count. Thus every output cover coordinate has a nonnegative itinerary expression in the marked factor’s cover counts, with the other factors fixed. Hence the $B[U]$ side never decreases relative to the $A[U]$ side. Any strong quotient with at least three vertices admits a Hamilton itinerary using the marked block twice and every other block once: start with a directed triangle through the marked vertex, repeat that vertex at the end, and insert all unused quotient labels. This strictly separates the output Hamilton counts for $n\ge2$. Every subsequent common homogeneous substitution retains that strict inequality because its one-run Hamilton term has positive coefficient. The result does not constrain target-dependent heterogeneous or different contexts.

## More individual cycles, fewer covers

For $U=TT_L$, let $C(y)$ count directed cycles by order, up to cyclic shift. The same finite quotient inventory gives
$$
 C_{A[TT_L]}(y)-C_{B[TT_L]}(y)=y^4[1+(1+y)^{L+1}].
$$
The difference on outside cycles is $2y^4+y^5$, and on mixed cycles it is $(y^4+y^5)[(1+y)^L-1]$, since each nonempty subset of the transitive run has one directed path. In fact replacing $(1+y)^L$ by $\sum_{S\subseteq U}H(U[S])y^{|S|}$, with $H(\varnothing)=1$, gives the general cycle difference as well (cycles wholly in $U$ cancel).

In the transitive-run case $A$ has no disjoint pairs of odd cycles; $B$ has the fixed triangle `bde` disjoint from the cycle `tax` with any odd-size nonempty subset of $TT_L$ filling $t$. This is the only pair type: two mixed cycles would both use $x$, the outside core has no disjoint pair, and the outside vertices left by a mixed cycle can contain a directed triangle only for a two- or three-outside-vertex quotient cycle. Thus adding individual cycles does not imply adding covers: it is the distinct packing compatibility that reverses the order.

No prime or $7p$ coverage conclusion follows: the equal counts in this particular switch are always multiples of 9.
