# Moore count coincidences, binary projections and a degree-57 target

**Status.** The tournament encodings and binary equivalences are exact. A degree-57 Moore graph is not constructed here. The known cases are calibrated by finite exact replays. Counts from different hosts are compared as integers, not as evidence of isomorphism.

## Two exact odd-count representations and which counts were already available

Use any fixed labels for a graph. For a subset S of labels, give each label a three-vertex cell: transitive if outside S, cyclic if inside S. Both cells share a prescribed two-edge spine and differ by one free arc. An ordinal join of the cells has Hamilton-path count $3^{|S|}$. The ordered, marked word (or a fixed-spine edit deck) retains the positions, while the scalar retains only the weight. For a graph representation use S as an open neighborhood, or as the XOR of two closed neighborhoods. In a Moore graph that XOR has size $2d-2$ on edges and $2d$ on nonedges. This uses literal Hamilton counts where the [previous prime-marker deck](../moore_path_cover_prime_markers_20261124/) used exponents of 47 in ordered two-cover counts at constant Hamilton count.

A second, strong readout of a positive weight $a$ uses a transitive $a$-chain as one block of a directed triangle, with two singleton blocks. Equivalently reverse only the extreme arc of a transitive tournament of order $a+2$. It has $2^a+1$ paths: each directed odd cycle contains both singletons and an odd subset of the chain, and any two such cycles intersect, so the odd-cycle packing identity gives $1+2\cdot2^{a-1}$. Positions can still be carried by the marked original deck; this strong scalar tournament by itself knows only $a$.

| Moore degree | Open row as product | Closed XOR as product: edge / nonedge | Open row as one strong tournament |
|---:|---:|---:|---:|
| 2 | 9 | 9 / 81 | 5 |
| 3 | 27 | 81 / 729 | 9 |
| 7 | 2,187 | 531,441 / 4,782,969 | 129 |
| 57 (conditional read) | $3^{57}=1,570,042,899,082,081,611,640,534,563$ | $3^{112}$ / $3^{114}$ | $2^{57}+1=144,115,188,075,855,873$ |

All product reads already have Chan–Pak witnesses using triangle joins, without any Moore input. For strong readout collisions, height 9 also comes from two cyclic triangles in ordinal join; height 129 comes from a cyclic triangle joined to a height-43 tournament (six-site mask 92, checked by [counts.py](counts.py)). These are strong versus nonstrong constructions of the **same integer**. At $a=57$,
$$
 2^{57}+1=3^2\cdot571\cdot174763\cdot160465489
$$
so even this large strong height is already constructible by taking the multiplicative closure of the [old finite Chan–Pak bank through $2^{35}$](../../../01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md). All the prime factors lie in range and are not 7 or 21. The larger strong pair reads $2^{112}+1$ and $2^{114}+1$ contain primes $358429848460993$ and $275415303169$, respectively, and thus cannot be assembled from factors drawn solely from that finite interval. Both are nevertheless built by the same strong extreme-arc family with no Moore graph. This distinguishes the limitations of **multiplication of that finite bank** from the reach of the all-width chain family.

For comparison, the earlier [compact and direct chart session](../moore_chanpak_counts_20261124/) found C5 heights 9, 11 and 15 and Petersen heights 5,213, 6,961 and 8,265. The new code-based values are additional chart contracts, not revised counts for those charts. No full compact Hoffman–Singleton height is inferred here.

### A collision that only paired readings can see through

For three distinct centers let $e$ be their number of induced edges and XOR their three open neighborhoods and, separately, their three closed neighborhoods. Write $a,b$ for the resulting weights. The two ordinal decks can be combined as one height $3^a5^b$: use cyclic triangles for the first read and any fixed height-5 four-site tournament (e.g. the extreme-arc construction with $a=2$) for the second. Each product is already in the Chan–Pak multiplicative closure.

| Triple | Open exponent $a$ | Closed exponent $b$ | Number of triples |
|---|---:|---:|---:|
| Independent, no common neighbor | $3d-6$ | $3d-3$ | remainder |
| Independent, with common neighbor | $3d-2$ | $3d+1$ | $n\binom d3$ |
| Exactly one edge | $3d-4$ | $3d-5$ | $(nd/2)(n-2d)$ |
| Two-edge path | $3d-2$ | $3d-3$ | $n\binom d2$ |

The equality behind the table is that XOR weight is the sum of the three set sizes minus twice the three pair intersections plus four times the triple intersection. Open pair intersections have size 0 or 1 and closed pair intersections size 2 or 1. A triple of independent vertices has at most one common neighbor; for a two-edge path the closed triple intersection is the middle vertex. The one-edge triples count each edge and the $n-2d$ vertices adjacent to neither endpoint. In particular **the closed readout alone confuses a two-edge path with an independent triple without a witness**; the open readout separates them, but has a different collision. At degree 57 the counts in the table order are 5,325,320,000; 95,095,000; 290,472,000; and 5,187,000. This is an exact histogram target on jointly indexed reads, not a new list of scalar heights conjecturally missing from Chan–Pak.

## A construction equivalence: lock the parity, then saturate one degree moment

Work over $\mathbb F_2$ at an **even** number $n$ of vertices. Let $U$ be an even binary code with $U\cap U^\perp=0$ (an even LCD code). Its orthogonal projector $P_U$ exists and is symmetric and idempotent. Its diagonal is zero: the squared norm of row $i$ equals its parity, zero, and also equals the $i$-th diagonal of $P_U^2=P_U$. It kills the all-one vector $\mathbf1$. Consequently
$$
 B=I+J+P_U
$$
is the adjacency matrix of a simple graph with odd degrees, and $B^2+B=J$ in characteristic two. Conversely every symmetric zero-diagonal solution of this last equation arises this way. Its diagonal equation forces $B\mathbf1=\mathbf1$; then $P=I+J+B$ is a symmetric idempotent with zero diagonal, its image is even LCD, and its image equals $\ker B$ (on $\ker B$, $J=0$ by the same equation, and $BP=0$). This is the complement form of the [known even-LCD projector correspondence](https://arxiv.org/abs/2407.07689), specialized to the Moore parity equations. With a generator matrix $G$ whose rows form a basis of $U$, one can build $P_U=G^{\mathsf T}(GG^{\mathsf T})^{-1}G$ directly.

Let $s_i$ be the ordinary integer degrees of $B$, and let $c_{ij}$ be the common-neighbor counts. The binary constraint says $c_{ij}\equiv1-B_{ij}\pmod2$ for distinct labels. Thus
$$
 \sum_i s_i^2-n(n-1)=4\sum_{i<j}\left\lfloor\frac{c_{ij}}2\right\rfloor\ge0.
$$
This follows by summing both the common-neighbor counts and their parity lower bounds over ordered distinct pairs, and using $\sum_i s_i(s_i-1)$ for total two-step counts. In particular:
- **Maximum-degree route.** At $n=d^2+1$, if all $s_i\le d$, the inequality forces all $s_i=d$ and every common-neighbor count to be its parity lower bound. The Moore property follows. In a regular root-branch search, proving parity is enough; there is no separate exact-one witness check once every degree is already correct.
- **Single moment route.** Equality forces triangle-freeness and exactly one common neighbor for every nonedge. The only nonregular equality graph is a star. Indeed writing the integer identity $B^2+B=D-I+J$, with $D=\operatorname{diag}(s_i)$, and commuting both sides with $B$, yields $(s_i-s_j)(1-B_{ij})=0$. Distinct degree classes are completely joined; triangle-freeness then forces exactly two independent degree classes, and the common-neighbor condition forces one of the two classes to have size one.

The star has binary adjacency rank two. A degree-57 graph has binary nullity 1520: the integral eigenvalues are $57,7,-8$ with multiplicities $1,1729,1520$; modulo 2 the zero-primary component is semisimple by $B(B+I)^2=0$. Thus at 3,250 sites the search is exactly: **find an even LCD code of dimension 1,520 whose complemented projector graph has $\sum s_i^2=10,559,250$**. This is equivalent to the Moore construction, not an existence claim. The moment can be read as the exponent of 3 in the product-cell tournament encoding all ordered pairs of neighbors at each center. Its raw height is not special: even the impostor star gives it. Rank plus parity plus the exponent is the certificate. The pentagon, of odd order and even degree, is not in this binary projector class.

### Constraint on a natural parity-preserving search move

For an even support vector $v$, $T=I+vv^{\mathsf T}$ is orthogonal and fixes $\mathbf1$. Hence $B'=TBT^{\mathsf T}$ preserves the binary equation and the LCD dimension. A weight-two $v$ just permutes two labels. For support $S$ of size four, a vertex outside $S$ with $t$ neighbors in $S$ has all four contacts toggled exactly when $t$ is odd. Its degree changes by $+2$ for $t=1$ and $-2$ for $t=3$.

At a Moore solution write $N_j=\#\{z\notin S:|N(z)\cap S|=j\}$ and $C=\sum_{z\notin S}\binom{|N(z)\cap S|}{2}$. The total contacts leaving $S$ are $4d-2e(S)$, and $C\le6-e(S)$ by counting pairs of $S$. Subtracting $2C$ gives $N_1-3N_3-8N_4\ge4d-12$; therefore $N_1+N_3\ge4d-12$ and $N_1-N_3\ge4d-12$. The move changes at least **216 outside degrees** at a putative degree-57 solution, and adds at least $16d-54=858$ to total degree (each boundary change counts twice and the four-site induced forest can lose at most three internal edges). It cannot be a solution-to-solution move there. The exact small replays find 20 neutral four-site supports for Petersen but none for Hoffman–Singleton, where the least observed squared-degree excess is 1,216. This warns against using elementary LCD transvections as fine endgame moves: coordinated larger changes would have to cancel degree changes, much as in the two-port Chan–Pak work.

See [the minimum-weight constraints](MINIMUM.md) for another filter on any code-driven search. Replay with **python3 counts.py** and **c++ -O3 -mpopcnt -std=c++17 census.cpp -o census && ./census**. The latter rebuilds both known odd-degree solutions from their codes, completely enumerates their weight distributions and checks every four-site transvection. Stored exact outputs are included.
