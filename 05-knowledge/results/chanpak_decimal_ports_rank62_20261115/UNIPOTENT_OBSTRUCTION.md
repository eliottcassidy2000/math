# The sharp no-go bound for an absolute append-only counter

Write $H(T)$ for the directed Hamiltonian-path count. Let
$$T_0\subseteq_{\mathrm{ind}}T_1\subseteq_{\mathrm{ind}}\cdots\subseteq_{\mathrm{ind}}T_r$$
be any chain of induced tournament extensions. A stage may append arbitrarily many vertices and choose contacts to any number of ports. Let $k$ be the number of stages at which $H$ strictly increases. Then, for $k\ge1$,
$$\boxed{H(T_r)\ge f(k+2),}$$
where $f(n)$ is the Moon–Busch minimum for a strong tournament of order $n$. The bound is sharp for every $k$. In particular a nonconstant polynomial absolute path-count readout is impossible under a repeatable append-only operation, whether its internal multi-port matrix is nonnegative, signed, or even finite-dimensional. A signed difference between two physical counts is not subject to this obstruction.

## Detecting strict stages

In the inherited odd-cycle packing expansion $H(T)=\sum_{\mathcal C}2^{|\mathcal C|}$, the sum is over collections of vertex-disjoint directed odd cycles, including the empty collection. An old packing persists under induced extension. Every vertex of a nontrivial strong component is in a directed triangle: a shortest directed cycle through it must have length three, since in a longer such cycle the chord from the chosen vertex to the second successor shortens one of the two possible resulting cycles through that vertex. A stage is therefore strict exactly when at least one new vertex lies on a triangle in the enlarged tournament. A new triangle gives a new positive packing contribution; if instead all new vertices are in no triangle, then they are singleton SCCs and contribute unit factors in the ordered strong-component product, without changing the old components.

Choose one fresh triangle vertex from each strict stage. In any final component $C$ hit by these representatives, look at the earliest chosen representative in $C$. Its witness triangle is already present and lies inside $C$. The two other vertices of that triangle were not chosen in an earlier stage (by minimality), cannot be chosen at the same stage (only one choice is made), and cannot be chosen later (they are already present). Thus $|C|\ge m_C+2$ when $m_C$ representatives lie in $C$.

## The floor arithmetic and attainment

The [Moon–Busch floor and interval chains](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md) give
$$f(3q)=3\cdot5^{q-1},\quad f(3q+1)=5^q,\quad f(3q+2)=9\cdot5^{q-1}\qquad(n\ge3).$$
Hence $f$ is strictly increasing. Let $g(m)=f(m+2)$ for $m\ge1$. The ratios $g(m)/5^{m/3}$ all lie in $[5^{1/3},9/5]$ by the three residue formulas. Since $(5^{1/3})^2>9/5$, we get $g(a+b)<g(a)g(b)$ for $a,b\ge1$. Therefore
$$H(T_r)\ge\prod_{C\text{ hit}} f(|C|)\ge\prod_{C\text{ hit}}g(m_C)\ge g(k).$$
Equality requires exactly one hit SCC, of order $k+2$ and at the strong floor, with all other SCCs singleton.

To attain the bound, choose a strong Moon–Busch minimizer on $k+2$ vertices. In a strong tournament, a directed cycle of length one smaller than the order induces a strong subgraph; use tournament pancyclicity repeatedly to delete down to a directed triangle, then delete any one vertex to obtain a two-vertex tournament. Reverse this deletion chain. Every newly added vertex lies on a triangle at its stage, so all $k$ stages are strict and the final count is $f(k+2)$.

This complements the narrower [no physical Heisenberg orbit in the fixed-determinant root-change grammar](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5172-hosts-and-cores-no-heisenberg-orbit-profile-congruences-and-the-cyclic-join-as-a-z-sqrt-minus-two-inner-product.md). An absolute polynomial counter must rewire or replace old arcs (or retreat to a difference or normalized port measurement). This result does not rule out an absolute exponential radix counter, e.g. a potential operation sending $H-1$ to $1000(H-1)$.

## Stationary repeated ribbons, including arbitrary span-three words

Fix a backward span $b$, and fixed words $P,U,S$ with $U\ne\varnothing$. Let $T_m$ be the ordered tournament encoded by $PU^mS$ with the bounded-backward-arc convention. Either $H(T_m)$ is eventually constant, or $H(T_m)\ge\exp(\Omega(m))$. This works at any bounded span and without a bottleneck-free hypothesis.

Indeed, a triangle on $i<j<\ell$ is either the forward-forward-backward cycle, requiring $\ell-i\le b$, or the backward-backward-forward cycle, requiring both gaps at most $b$. Every triangle thus fits in a window of width $2b$. If the bi-infinite word $U^{\mathbb Z}$ contains any triangle, it repeats in linearly many translates, so a long $U^m$ contains linearly many triangle-active vertices. Strong components of size at least three have $H(C)\ge3^{|C|/3}$ by $f$, so the SCC product is exponential in their total size. Otherwise every triangle lies in a bounded fringe at the two ends. Their induced tournament stabilizes for sufficiently large $m$: the internal words near both ends stabilize, and long cross arcs are forward. Every vertex outside both fringes is a singleton SCC and can be deleted without affecting $H$. This rules out any nonconstant polynomial contribution as the largest visible growth of the full physical count under stationary ribbon iteration. Signed differences can cancel a faster visible term and are not ruled out.


## The decimal ladder is impossible for a stationary single unmarked module power

Fix a pointed host $(S,z)$ and one nonempty tournament $U$. Replace $z$ with the ordinal sum of $k$ copies of $U$, for every $k\ge1$. An exact all-height identity
$$H(S[z\gets U^{\to k}])=1+B^k,\quad B\ge2,$$
is possible within this architecture if and only if $B$ is a positive power of two. In particular the architecture cannot give the base-1000 ladder.

Indeed, let $a_r$ count quotient words that visit the hole $r$ times. The [unmarked substitution interface](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_symmetry_interfaces_20261115/OBSERVABLES.md) gives $(r-1)!\mid a_r$, as the $r-1$ intermediate excursions have distinct actual exterior vertices and may be permuted freely. For the ordered nonempty slot counts $c_r$, the exact output is $\sum_r a_r c_r(U^{\to k})$. The optional slot polynomial satisfies $P_{U^{\to k}}(j)=P_U(j)^k$ by ordinal concatenation in each color. Using binomial inversion,
$$c_r(U^{\to k})=\sum_{j=1}^r(-1)^{r-j}\binom{r}{j}P_U(j)^k.$$
The bases $P_U(j)$ are positive and strictly increasing for $j\ge1$: all coefficients of the odd-cycle-packing formula for $P_U(x)$ are nonnegative and the exponent of the empty-packing term is at least one. If $R$ is the largest index with $a_R\ne0$, the largest base has coefficient $a_R$ with no other contributions. Distinct positive exponentials are linearly independent (or successively divide by the largest base and take a limit), so $R$ must have $a_R=1$, and the factorial implies $R\le2$. An $R=1$ output has only one exponential and cannot work. At $R=2$ the output is
$$P_U(2)^k+(a_1-2)H(U)^k.$$
Hence $P_U(2)=B$, $a_1=3$, and $H(U)=1$. A tournament with one Hamiltonian path is transitive by the odd-cycle packing formula, so $P_U(2)=2^{|U|}$. Conversely take the pointed triangle host and $U=TT_d$; it produces $1+2^{dk}$.

There is nevertheless an exact normalized observable for $1000^k$: $P_{TT_{3k}}(10)=1000^k$ is just ten-way coloring of a transitive spine. In any fixed signed **integer** linear combination of unmarked host outputs on this module, the coefficient of the $10^{3k}=1000^k$ mode is divisible by $9!$. Only $r\ge10$ can contribute to this mode, and every $a_r$ is divisible by $(r-1)!$. Thus the slot mode cannot be simply relabeled as an unweighted Hamiltonian readout, even when negative combinations of finitely many unmarked hosts are allowed. A marked construction would have to defeat the excursion permutation tax, or a different architecture would have to realize the ladder.

