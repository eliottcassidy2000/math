# A graph-polynomial certificate for arbitrary three-lists

**PROVED**, by finite interpolation below. The certificate is sufficient;
failure to find it is not a proof of a list obstruction. The finite census
and its independent coefficient checks have separate **FINITE-EXACT** status.

## L1. The coefficient identity

Let $t_1,\ldots,t_n$ be nonnegative integers. For each $i$, let $S_i$
be a set of $t_i+1$ distinct rational numbers and put

\[
D_i(a)=\prod_{b\in S_i\setminus\{a\}}(a-b).
\]

For every $P\in\mathbb Q[x_1,\ldots,x_n]$ of total degree at most

\[
T=\sum_i t_i,
\qquad
[x_1^{t_1}\cdots x_n^{t_n}]P
=\sum_{a\in\prod_i S_i}\frac{P(a)}{\prod_i D_i(a_i)}.                 \tag{1}
\]

**Proof.** Interpolating the univariate polynomial $x^k$, for
$0\le k\le t_i$, on $S_i$ and comparing the coefficient of $x^{t_i}$
gives

\[
\sum_{a\in S_i}\frac{a^k}{D_i(a)}=
\begin{cases}0&k<t_i,\\1&k=t_i.\end{cases}                         \tag{2}
\]

Both sides of (1) are linear in $P$. Check a monomial $x^u$ of total
degree at most $T$. If $u=t$, its weighted sum is 1 by (2). Otherwise
some $u_i<t_i$: if every $u_i\ge t_i$, at least one strict inequality
would make the total degree exceed $T$. The weighted sum factors over
the variables, and this coordinate makes it zero by (2). This proves
(1), even when another individual exponent is greater than its $t_i$.
If $t_i=0$, $S_i$ is a singleton and its denominator is the empty
product 1. The zero-variable case is the same empty-product identity. QED.

The total-degree bound matters. For one variable, $S=\{1,2\}$, $t=1$
and $P=x^2$, the right side is 3 whereas the coefficient of $x$ is 0.

## L2. The list-colouring consequence

Let $G$ be any finite simple graph on labelled vertices $1,\ldots,n$,
with $m$ edges, and use the fixed polynomial

\[
P_G(x)=\prod_{\{i,j\}\in E(G),\ i<j}(x_i-x_j).
\]

If $[x^t]P_G\ne0$ and $0\le t_i\le2$ for every $i$, then $G$ is
three-choosable: every assignment of lists of at least three colours has
a proper colouring. More precisely, it is colourable from any lists of
sizes at least $t_i+1$.

**Proof.** The polynomial is homogeneous of degree $m$, so a nonzero
coefficient has $\sum_i t_i=m$. From each list choose $t_i+1$ colours.
The union of these finite selected lists is finite; embed the actual colour
names injectively into $\mathbb Q$, preserving equality of names across
different lists. Call the resulting sets $S_i$. Identity (1) says a
nonzero number is a weighted sum of values $P_G(a)$, so at least one
value is nonzero. All its edge factors are nonzero, exactly saying that
the selected colours give a proper colouring. Return to their original
names. No bound or common palette is assumed. QED.

Changing edge-factor orientations multiplies the entire polynomial by one
global sign. Relabelling vertices permutes exponent vectors and may change
this sign, so certificate existence is an isomorphism invariant. A certificate
does not assert anything about the geometry of a graph realization. It says
that **if this graph is the private graph and each private vertex has at least
three allowed colours, all the lists extend**, regardless of their overlap.

The edgeless graph has the certificate $t=0$, coefficient 1. The triangle
has certificates, although its coefficient of $x_1x_2x_3$ is zero. The
complete graph $K_4$ has no certificate with all exponents at most 2:
its Vandermonde polynomial has exponent permutations of $0,1,2,3$.
The common list $\{0,1,2\}$ is its direct negative control.

## L3. What the census computation certifies

The producer in [../cores/census.cpp](../cores/census.cpp) multiplies the
edge factors in the quotient by the monomial ideal
$(x_1^3,\ldots,x_n^3)$. Discarding a term whose exponent exceeds 2 is
safe: further multiplication can only increase exponents, so that term can
never contribute to a retained coefficient. Every retained nonzero term
is therefore a coefficient of the original homogeneous polynomial.

The independently implemented [check.py](check.py) instead chooses one
endpoint from each edge, exhaustively enumerating all $2^m$ choices.
A choice contributes its endpoint multiplicities and the sign
$(-1)^{\#\text{larger endpoints chosen}}$. Only at a completed choice does
the independent code discard exponent vectors outside $\{0,1,2\}^n$.
It compares the complete remaining support size, the first coefficient in
the producer's base-three order, and the selected certificate, for every
one of the 38 order-eight representatives. Thus no intermediate capped
polynomial recurrence is shared with the producer.

The independent receipt confirms 30 certificates and zero capped support
for the eight residual masks

```
5675891 5683827 5691763 5692274
6215909 6216436 14587370 48675698
```

These eight require additional arguments. Neither zero support nor passing
the necessary metric filters implies geometric realizability, failure of
three-choosability, or a boundary obstruction. The exhaustive graph-universe
claim belongs to the census and its own orbit/universe audit; this checker
independently audits all coefficients of its supplied representative bank.

The listed graphs have 12 to 16 edges. Literal expansion contributes at
most $2^{16}$ signed terms to a coefficient. Even for an arbitrary
eight-vertex simple graph, $2^{28}$ bounds every intermediate coefficient
in absolute value, so the producer's signed 64-bit arithmetic is ample.

Reproduce from the worktree root:

```powershell
python 05-knowledge/results/two_port_eight_core_20261103/recovery/check.py
python -O 05-knowledge/results/two_port_eight_core_20261103/recovery/check.py
```

The checker also verifies (1) on all monomials under its stated degree bound
for small exponent vectors, checks the necessary degree-bound hostile, and
checks the graph-polynomial interpolation directly on all labelled graphs
through four vertices. Its tests use explicit exceptions, not Python
assertions, so optimized execution retains the checks.
