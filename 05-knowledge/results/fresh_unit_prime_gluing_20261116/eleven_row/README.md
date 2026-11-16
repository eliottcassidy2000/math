# Induced prime gluing: the entire 11 row fails

This extends the [prime-gluing and fresh-deposit packet](https://github.com/eliottcassidy2000/math/blob/9792e615c40d8dc3316f06de950fcad546a024e5/05-knowledge/results/fresh_unit_prime_gluing_20261116/README.md). The new result concerns gluing **with both summands still present as induced subtournaments**. It is not an obstruction to count realization or to a construction that rewires a summand.

**Theorem.** For every prime $q\ge11$, no tournament of count $q+10$ contains induced tournaments of counts $q$ and $11$.

Write $H$ for Hamiltonian-path count and $\tau$ for directed-triangle count. Use the [Grinberg–Stanley positive packing expansion](https://arxiv.org/abs/2307.05569) of $H$ over disjoint directed odd cycles and the earlier [triangle-budget theorem](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md): if $H(T+x)-H(T)=2r>0$, then $\tau(T)\le r$. Any new directed triangle gives a strictly positive increase, and all counts are induced-monotone. Tournament strong components form an ordinal sum, with counts multiplying. A tournament with all directed triangles contained in a specified set has only singleton strong components outside that set, because every vertex of a nontrivial strong tournament lies in a directed triangle.

Start with the induced $q$-witness and add the remaining vertices in any order. Since $q$ is prime there is one nontrivial strong component $C$ of count $q$. Up to the first strict addition $x$, the only extra components are singleton pads: a neutral addition adds no directed triangles. At the strict addition $\tau(C)\le5$. As a strong tournament on $s\ge3$ vertices has at least $s-2$ directed triangles, $C$ has order at most seven.

*Uniform contact.* If $x$ sees $C$ uniformly, $C$ is still a module. Contract it to obtain a quotient with a directed triangle, since every new triangle uses $x$ and one cannot use two vertices in $C$. Every Hamilton path of the quotient lifts to $q$ paths by putting a path through $C$ in its quotient position. The quotient count is at least three, so the new count is at least $3q>q+10$.

*Mixed contact.* A mixed contact creates a directed triangle through $x$ inside $C+x$: a path in $C$ from an outneighbor of $x$ to an inneighbor of $x$ crosses the cut. The positive packing expansion shows that the local increase $H(C+x)-q$ is at most the full increase at this step. Count all strong cores on 3–7 vertices with at most five directed triangles, and enumerate every mixed orientation of the new vertex. The following are **all** prime core counts at least 11, with their minimum local half-increases:

| Core count $q$ | 11 | 13 | 17 | 19 | 23 | 29 | 31 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Minimum local half-increase | 4 | 5 | 5 | 6 | 9 | 8 | 10 |

The case $q=11$ would give the already forbidden count 21. For $q=13,17$, the local increase spends the entire budget of ten, so there is no new triangle outside $C+x$, and every later insertion is neutral. Thus the host is $C+x$ with singleton pads: here $C+x$ is strong by mixed contact. A standalone scan of every induced subset of each of the qualifying $C+x$ gives no count 11 (12 degree-ordered core/contact pairs for each $q$). Padding cannot change an induced count, finishing the proof.

## Reproduction

```sh
c++ -std=c++17 -O3 check_eleven.cpp -o /tmp/check_eleven
/tmp/check_eleven
```

Only cores labeled in nondecreasing outdegree order are kept, which retains a labeling of every isomorphism class. The checker uses direct Hamilton-path subset DP for every contact, not the rooted-pack update used to discover the result. It checks the table, the scope of the table, and both final subset exclusions. No stored atlas is needed.

## A stronger triangle record bound

In the prior [triangle decoder](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md), recording the first vertex of every directed triangle injects the Hamilton paths into the record space. The record coordinates are not independent. Group any $k$ distinct triangles sharing an arc $a\to b$. One cannot select $a$ in one record and $b$ in another. The two allowed families (no $a$ chosen and no $b$ chosen) each have $2^k$ records, with precisely one in common. For any **partition** of the triangles into such books of sizes $k_i$, this proves
$$
                          H\le\prod_i (2^{k_i+1}-1).
$$
One can always choose books of size one or two; larger books give an additional tool to use where cycles are concentrated on a few arcs. A one-backward-arc tournament with $m$ intervening vertices has all $m$ triangles in a single book and is bounded by $2^{m+1}-1$ instead of $3^m$.

For a uniform bound on strong tournaments, match triangles that share an arc in the graph of triangles. If a maximum matching has size $M$, then $H\le3^\tau(7/9)^M$. Each matched edge sees at most two *distinct unmatched* triangles. Indeed, if each endpoint has an unmatched neighbor, distinct choices at opposite ends give a length-three augmentation. If only one endpoint has unmatched neighbors, its two arcs not shared with its partner can support at most one each, because unmatched vertices are independent.

An isolated triangle in the overlap graph is a module: mixed contact to a directed triangle creates a new triangle on one of its arcs. These isolated triangles are disjoint. If two intersected in exactly one vertex, substituting another vertex of the first module for that common vertex would give a neighbor of the second sharing its other two vertices; intersection in two already shares an arc. In a strong tournament of order at least four their simultaneous contraction leaves a nontrivial strong quotient. For every contracted vertex choose a quotient triangle through it. A distinct chosen quotient triangle with $t$ contracted vertices lifts to $3^t\ge3t$ nonisolated triangles. Distinct quotient triangles have disjoint lift sets, and the choices cover the contracted vertices, so at most $\tau/4$ triangles were isolated. All nonisolated unmatched vertices see some matching edge, implying $4M\ge3\tau/4$. Consequently
$$
 H(T)\le \bigl(3(7/9)^{3/16}\bigr)^{\tau(T)} \qquad(\text{strong }T,\ |T|\ge4).
$$

For the general fixed-summand problem, induced $p+q-1$ with prime $q\ge p\ge11$ forces $q\le[3(7/9)^{3/16}]^{(p-1)/2}$. This follows by applying the triangle budget at the first strict insertion and using the strong core of the old $q$-witness. More fundamentally the same budget bounds the possible host cores: for $r=(p-1)/2$ there are at most $2r$ triangles in the host, because new triangles cost at least one per unit of half-increase. Keep their at most $6r$ vertices and remove the singleton pads. The induced witnesses retain their counts on the retained vertex set. So every fixed induced row is decidable by a finite search; more targets on that row cannot fix an induced obstruction. The positive direction still needs lifts that change the source while controlling the continuation data. None of these bounds obstructs such a lift.
