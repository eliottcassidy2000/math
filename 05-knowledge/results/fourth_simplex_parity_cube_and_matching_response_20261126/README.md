# Fourth-simplex parity cubes and the marked-cycle response law

Continues the [seven-vertex Pollock–Chan–Pak parity cube](../pollock_chanpak_literal_parity_cube_20261126/). This packet proves sharp matching bounds in all response orders and an explicit formula for the highest surviving path-cover coefficient. It gives a minimal literal strong four-switch cube and three exact fourth-simplex endpoint matches. It does not assert a pure parity construction in every order or solve either conjecture.

## 1. Nine vertices and four disjoint switches

Use $A_i,B_i$ for $i\bmod4$, and $s$. Each block has the directed square $i\to i+1$, with opposite chords $0\to2$ and $1\to3$ at the base. Those four chords are the marked switches. Use $A_i\to B_j$ precisely when $i,j$ have the same parity, and set $s\to A_i$, $B_i\to s$. All 16 switches give strong tournaments: each block contains its directed square and the three parts form a cyclic chain of reachability.

Rotation by two in either block toggles both switches in that block. Simultaneously rotating both blocks by one preserves every unmarked direction and permutes the two switches in each block while changing the parity within each block. These operations generate eight permutations. For each fixed total parity they can reach either pair of block parities, and rotations by two act transitively on the two choices within each block of that parity. Thus the two parity sides are two literal isomorphism orbits of the switch family.

Subset dynamic programming gives **2917 for even parity and 2919 for odd parity**. All ordered path cover counts other than the one-path count agree. In order from two through nine runs they are
\[
38622,\ 232422,\ 821040,\ 1857000,\ 2751840,\ 2610720,\ 1451520,\ 362880.
\]
The agreement is forced by the general response bound below, not accidental in this example.

## 2. The sharp matching bound

Let $A_j(T)$ count ordered covers by $j$ nonempty directed paths and put $P_T(z)=\sum_j A_j(T)\binom{z}{j}$. Choosing the occupied colors in increasing order identifies this with the principal specialization of [Grinberg–Stanley's Rédei–Berge symmetric function, Theorem 1.38](https://arxiv.org/html/2307.05569). Under the convention where the weak color runs follow complement arcs, reverse each path in its own run. Thus the odd directed cycle expansion specializes to
\[
 P_T(z)=\sum_{\mathcal C}2^{|\mathcal C|}z^{n-\sum_{C\in\mathcal C}(|C|-1)},
\]
over packs of pairwise vertex-disjoint directed odd cycles of length at least three; vertices outside the pack are singletons.

Mark a matching $E$ of $r\ge1$ arcs and write $\Delta_E f=\sum_{S\subseteq E}(-1)^{|S|}f(T^S)$, where $T^S$ reverses exactly $S$. A pack not using a particular marked arc cancels in the full contrast, pairing the two states of that switch. If a cycle uses $m$ marks it needs at least $2m$ vertices because the marks are disjoint, hence at least $2m+1$ because the cycle is odd. Therefore $\Delta_E P_T$ vanishes identically on $n\le2r$, and otherwise has degree at most $n-2r$. By the binomial basis, $\Delta_E A_j=0$ for $j>n-2r$. At $n=2r+1$ only Hamiltonian odd cycles containing all marks can contribute, and $\Delta_E P_T=z\,\Delta_E H$. In a pure parity cube the alternating contrast is $2^{r-1}$ times the difference of the two sides, so all the higher path-cover layers of the nine-vertex example really are preserved.

The threshold $2r+1$ is **sharp for a nonzero contrast** for every $r$. At that order, thread the marked pairs as blocks along an undirected Hamiltonian path, with the one extra vertex at an end. Expand $H$ in orientation signs $x_e\in\{1,-1\}$ on the complete graph. The degree-$2r$ monomial on this path occurs in its unique undirected Hamiltonian cycle completion. Both orientations of that odd cycle contribute the same nonzero sign to the monomial, since its degree is even. Every non-Hamiltonian pack contributes only degrees at most $2r-2$: for each odd cycle sum its two orientation indicators, which cancels odd-degree terms and leaves degree at most its length minus one. No second Hamiltonian cycle contains all edges of the path. Hence the chosen monomial survives, contains every marked arc, and guarantees a nonzero full contrast in some orientation (by uniqueness of the multilinear sign expansion). Ordered join with a transitive block multiplies $P_T$ by a power of $z$ and keeps the degree bound sharp. This proof does not guarantee only two parity values or strong connectivity for every $r$; the explicit construction does both for $r=4$.

### Exact highest layer by distributing unmarked vertices

Let $U$ be the vertices not incident to the matching. For a nonempty subset $J$ of marks and $u\in U$, let $h(u,J)$ be the full matching contrast of $H$ in the subtournament induced by the endpoints of $J$ and $u$, with the reference directions inherited from $T$. It has $2|J|+1$ vertices. Then
\[
 [z^{n-2r}]\Delta_E P_T(z)
 =\sum_{\pi\text{ a partition of }E}
   \ \sum_{\iota:\pi\hookrightarrow U}
   \prod_{J\in\pi}h(\iota(J),J).
\]
Equality in the degree estimate forces every contributing cycle to use exactly one unmarked vertex beyond its marked endpoints, and forbids any additional nontrivial cycle without marks. The cycles therefore partition the marks and use distinct extra vertices. On each support the minimal-order result says that $h$ counts precisely twice its signed Hamiltonian-cycle response; switch states on disjoint blocks factor. This proves the formula including the required no-collision condition on unmarked vertices.

For two marks $e,f$, the answer is the sum of two different geometries: $\sum_u h(u,\{e,f\})$, the one-extra-vertex joint term, plus $\sum_{u\ne v}h(u,\{e\})h(v,\{f\})$, the disjoint triangle term. Dropping $u\ne v$ would create fictitious cycle packs. The top binomial coefficient change $\Delta_E A_{n-2r}$ is $(n-2r)!$ times the displayed power coefficient. The replay checks this independently by ordered subset partitions on 108 matching backgrounds. The formula is an explicit way for lower, local marked objects to control a surviving macro coefficient, not an assertion about absolute Hamiltonian totals.

## 3. Same module gap, different starting scale

Replace $s$ by any nonempty tournament $M$ with common external directions. Marked switches and strong parity orbits persist. The count on a parity side is a sum over placements of runs of $M$ in an ordering of the eight outside vertices, weighted by the ordered path covers $A_j(M)$. A backward adjacency outside requires a module slot; regard it as compulsory if it has the proper entry and exit arcs, and reject that order otherwise. If there are $b$ compulsory and $f$ admissible discretionary slots, add $\binom{f}{j-b}$ to run coefficient $j$. Enumerating the $8!$ orders gives the even coefficient row
\[
 (2917,4638,3274,984,96)
\]
for $A_1$ through $A_5$; the odd row adds $(2,1,0,0,0)$. At most five runs can occur because each nonfinal run must exit to a distinct $A$ vertex. The exact low count in the $P_M$ basis is
\[
7P_M(1)-240P_M(2)+298P_M(3)+504P_M(4)+96P_M(5).
\]
The absolute gap is $2A_1(M)+A_2(M)=P_M(2)$, exactly the same as in the three-switch family. For a transitive $k$-vertex module $P_M(j)=j^k$ (assign each vertex a color), so both gaps are $2^k$ while the low counts have respective leading scales $5^k$ and $3^k$. Thus identical local response functions on **every** module need not agree on where the response begins.

## 4. Fourth-simplex endpoints and a sharp odd-carry arity law

Write $G(x)=\binom{x}{4}$. The even and odd index lists of $t+\sum_{i\in S}d_i$, for subsets of four directions, each have eight entries and agree in the first three moments. The total gap in $G$ is $\prod_i d_i$ by finite differences. Adding common indices leaves all this data intact. The following are exact matches under $H=2N+1$:

| Module | Tournament endpoints | Fourth-simplex endpoints | Base, directions; common indices |
|---|---|---|---|
| Singleton | 2917, 2919 | 1458, 1459 | $7;(1,1,1,1);(4,4,9,10)$ |
| Cyclic triangle | 56223, 56235 | 28111, 28117 | $2;(1,1,1,6);(4,30)$ |
| Transitive four | 209329, 209345 | 104664, 104672 | $7;(1,2,2,2);(24,40)$ |

With four directions the even subset list gives the high simplex value, so the parity labels interchange between the two problems. At the first match, the changing simplexes alone give 1120 and 1121. Direct subset DP checks all 48 literal tournament orientations and exact comparison checks all three simplex pairs. These are fourth-simplex sums, not five-tetrahedral Pollock representations.

The difference in resources is sharp. For degree $d\ge1$, suppose two equal-length lists have equal power sums in degrees 1 through $d-1$, but different sums of $\binom{x}{d}$. Cancel any common entries; put $W(z)=\sum z^{x_i}-\sum z^{y_i}$ for the resulting lists of length $m$. Lower binomial moments imply $W(z)=(z-1)^d R(z)$ in $\mathbb Z[z]$, and the gap is $R(1)$. If it is odd then $R(-1)$ is odd. Thus $2^d\le |W(-1)|\le2m$. This proves necessity of **$2^{d-1}$ changing indices per side**, and at equality one whole list has even indices and the other odd. A $d$-direction parity carry with $d-1$ unit directions and one direction equal to the desired positive odd gap attains the bound. This generalizes the earlier cubic four-contact obstruction and gives eight for degree four.

More generally for an integer-valued degree-$d$ polynomial $f$, with leading Newton coefficient $a_d=d![x^d]f$, a lower-moment-preserving gap is $a_d R(1)$. Odd gaps are impossible if $a_d$ is even. If $a_d$ is odd, every odd multiple of $|a_d|$ has exactly the same sharp arity (exchange the sides for the sign). For ordinary $m$-gonal numbers $a_2=m-2$, for tetrahedra $a_3=1$, and for octahedra $a_3=4$. In particular the same guard permits odd carries for odd polygons and tetrahedra but forbids them for octahedra; it is the Newton coefficient, not geometric dimension alone, that explains the difference. Nonzero tournament responses require a linear vertex threshold $2r+1$, whereas moment-matched odd figurate carries cost exponentially many active indices. Matching gap units and matching absolute targets are therefore useful comparisons, not interchangeable constructions.

Replay: `node replay.js`.
