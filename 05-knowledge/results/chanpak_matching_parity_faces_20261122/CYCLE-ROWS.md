# A marked cycle reads the relative permutations of Latin rows

**PROVED.** Give each labeled tournament a weight equal to its number of directed Hamiltonian cycles. The orientation parities of disjoint perfect matchings then record exactly which pairs of the matchings form Hamiltonian cycles. When the rows come from a Latin square, the correlation graph is its row-Hamiltonian relation. The [path-marked law](README.md) on the same space has independent row parities.

## The cycle-marked law

Fix an even $n\ge4$, a reference orientation of $K_n$, and edge-disjoint perfect matchings $M_1,\ldots,M_k$. For a tournament $T$ write $p_i(T)$ for the product of its orientation signs on $M_i$. Draw a uniform directed cyclic order of the $n$ vertices up to rotation, and orient every other pair independently and fairly. This yields a uniform pair $(T,C)$ with $C$ a marked directed Hamiltonian cycle of $T$.

Define a graph $G$ on the $k$ matching rows. Join $i$ and $j$ if $M_i\cup M_j$ is an undirected Hamiltonian cycle. On each such edge let $\epsilon_{ij}$ be the product of the orientation signs along a traversal of that cycle. Reversing the traversal changes $n$ signs and leaves the product unchanged.

**Theorem.** For every $u\in\{\pm1\}^k$,
$$
\Pr(p=u)=2^{-k}\left(1+\frac{2}{(n-1)!}
            \sum_{ij\in E(G)}\epsilon_{ij} u_i u_j\right).
$$
In particular $G$ is recovered exactly as the support of pairwise parity correlations. The row law is uniform precisely when every pair of rows has disconnected union (the anti-perfect case for a factorization). A perfect factorization gives complete correlation support. The only moments here are squarefree parity characters, with each row used zero or one time.

**Proof.** Conditional on a planted cyclic order, the mean of a parity character is zero whenever a coordinate in the character is still free. Three or more rows cannot fit in a single cycle because their union has degree at least three. A fixed single row has mean zero by exchanging the vertex labels of one of its matched pairs. This bijection of marked pairs reverses the sign of that row and permutes its other orientations identically. Two rows form a degree-two spanning graph, so their union can fit only if it is the planted cycle itself. In the connected case two of the $(n-1)!$ directed cyclic orders fit, with the same sign. These are all nonconstant Fourier moments. Invert the Fourier transform on the $k$ row bits. ∎

**Magnitude and sampling.** With $U$ uniform on the parity cube the exact chi-square divergence is
$$
\mathbb E_U\left[(2^k\Pr(p=U)-1)^2\right]=4|E(G)|/((n-1)!)^2,
$$
by character orthogonality. Moreover $k\le n-1$, so every sector mass relative to uniform is within $1/(n-3)!$ of one. At $n\ge6$ every target sector can therefore be sampled with at most $4/(1-1/(n-3)!)$ independent trials in expectation. Draw a cyclic order and let $s$ be its number of saturated rows. Reject if a saturated parity is wrong. Otherwise retain the order with probability $2^{s-2}$ and fill all unsaturated row parities by free pivots as in the path sampler. Since $s\le2$, this is a valid gate. For a compatible order the proportion of completions with the target is $2^{s-k}$, so the gate weights orders proportionally to the number of available completions. Its mean acceptance is $2^{k-2}\Pr(p=u)$, giving the bound. At $n=4$ the lower bound is zero, and the full matching factorization has sectors of zero cycle-marked mass.

## Universal mixed flips for the cycle count

**Theorem.** Let $C(T)$ count directed Hamiltonian cycles up to rotation on $n\ge3$ vertices. The mixed flip difference $\Delta_D C$ is nonzero for some outside orientation precisely in the following cases. At odd order $D$ is a linear forest. At even order $D$ is a linear forest or a spanning cycle.

**Proof.** The elementary Fourier criterion in [Theorem 1 of the path note](README.md) applies to every potential. A nonzero coefficient of $C$ can use only edges of one undirected Hamiltonian cycle, by expanding each directed cyclic-order indicator. A proper subset of such a cycle is a linear forest. The coefficient at a full cycle is nonzero exactly at even order, since the two cyclic traversals contribute opposite signs at odd order and equal signs at even order. Every linear forest can be extended to a spanning path by joining endpoints of its components, and this path closes in $K_n$. At even order the coefficient of the resulting full cycle witnesses visibility. At odd order take its spanning path alone. There are just two directed cyclic orders containing that many specified undirected edges, and their contributions agree because the path has even length. Thus that coefficient witnesses visibility at odd order. ∎

At odd order the two count potentials have identical families of detectable flip sets. At even order the mixed cycle difference can see every spanning perfect matching and every full Hamiltonian cycle, both of which the path count annihilates identically.

## Which path information was discarded?

In the uniform path-marked experiment let $B$ be the event that the final vertex points back to the first. Since $n\ge4$, this arc was free when the path was planted, so $\Pr(B)=1/2$. Each marked cycle yields exactly $n$ such path marks by choosing the first vertex. The row law conditional on $B$ therefore equals the cycle-marked law. The [path law](README.md) is uniform on the parity sectors, so Bayes' rule gives
$$
\Pr(B\mid p=u)=\tfrac12+\frac{1}{(n-1)!}\sum_{ij\in E(G)}\epsilon_{ij}u_i u_j .
$$
Consequently closing and nonclosing paths contribute opposite row correlations, which cancel when the endpoint-closure bit is forgotten. Keeping just this one bit recovers the Hamiltonian unions of the rows. On every parity sector at $n\ge6$ the fraction of marked paths that close differs from one half by at most $1/(2(n-3)!)$.

## A literal Latin-square specialization

Use column and symbol labels as separate sets of $m\ge2$ vertices. Each row of a Latin square is a perfect matching, represented by a permutation $\sigma_i$ from columns to symbols. All rows are edge-disjoint. Embed these matching rows in $K_{2m}$, referencing their arcs from columns to symbols. Traversing $M_i$ forward and $M_j$ backward induces $\sigma_j^{-1}\sigma_i$ on the columns. It follows that the two-row cycle-marked correlation is
$$
\mathbb E[p_i p_j]=
\begin{cases}
2(-1)^m/(2m-1)!,&\sigma_j^{-1}\sigma_i\text{ is an }m\text{-cycle},\\
0,&\text{otherwise}.
\end{cases}
$$
There are exactly $m$ backward steps on the undirected alternating Hamiltonian cycle. Thus this tournament measure recovers the row pairs forming a single cycle and distinguishes perfect from anti-perfect factorizations. Row labels and the two shores are retained. The law loses the lengths of cycles in a disconnected row union, as well as all information about the Latin square beyond its row-Hamiltonian relation. The relative-permutation characterization is classical. [Allsop and Wanless, *Row-Hamiltonian Latin squares and Falconer varieties*](https://arxiv.org/abs/2211.13826) · [Mira Cult, pinned matching-row sampling writeup](https://github.com/Mira-Cult/math/blob/30ddaaad27b48df2002225c3e3f40cffd02c5253/notes/combinatorics/exact-latin-square-sampling/README.md)

Two predictions follow directly. When $m$ is even a full $m$-cycle is an odd permutation, so the correlation graph is bipartite by row permutation sign. For the addition table modulo $m$, the correlation graph joins $i$ and $j$ exactly when $\gcd(i-j,m)=1$, since the relative permutation is translation by their difference. In prime order every pair is detected. These statements concern a fixed Latin square embedded in the tournament measure. They require no assumption on the complexity of sampling Latin squares.

## Finite checks

Run `python3 05-knowledge/results/chanpak_matching_parity_faces_20261122/cycle_probe.py` from the repository root. The [script](cycle_probe.py) directly enumerates all cycle-marked pairs for the complete one-factorizations at orders four and six. For cyclic Latin squares of orders two through four, it enumerates all directed cyclic orders of the larger tournament and integrates their remaining fair orientations. [Exact output](cycle_probe.out)
