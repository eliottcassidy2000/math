# When two independently chosen module widths must follow one perfect-number ray

Write $E_p=2^{p-1}(2^p-1)$ for any integer $p\ge2$. A strong quotient is a fixed strong tournament in which two distinct vertices are replaced by transitive chains of independently chosen positive widths $m,k$. A visit means a maximal contiguous run in one chain along a directed Hamilton path. The **two-visit condition** requires at most two visits to each chain for all choices of widths.

**PROVED.** For a fixed strong quotient satisfying the two-visit condition and any fixed integer $e$, every sufficiently large exact solution of $H(m,k)=E_p+e$ uses a unique pair of constant shifts whenever the quotient admits a ray. The shifts are computed from at most three candidates, with $m=p-1+i, k=p-1+j$. The successful pair realizes the entire ray once both widths are positive. The threshold is effective.

**FINITE-EXACT with PROVED continuation.** For strong quotients on at most ten vertices under the two-visit condition, the lower neighbor occurs exactly at $2\le p\le12$ and the upper exactly at $3\le p\le13$. The upper endpoint at exponent 13 uses ports 1 and 9 with widths 6 and 9 in the ten-vertex quotient with upper triangle 101011011111111111100101111001111001111111111, giving 33,550,337 on 23 vertices. In particular the upper and lower neighbors at $p=31,61$ both fall outside this grammar. For every $p\ge14$, an exact straddling flip would have to leave this grammar at both endpoints. The earlier complete search through eight vertices allows arbitrary visit multiplicities at $p=61$. [Earlier unrestricted-width census](../chanpak_decimal_kernel_20261122/TWO-PORT-CENSUS.md) · [Earlier fixed-shift coefficient criterion](../chanpak_absolute_liouville_leech_20261122/MODULE-OFFSETS.md)

## A finite test for independent widths

The four counts $c_{rs}$ enumerate spanning words through the quotient with $r,s\in\{1,2\}$ visits to the two ports, all ordinary vertices used once, and no consecutive copies of either mark. A transitive chain has one way to fill a single run and $2^m-2$ ways to fill two runs. Thus
$$
 H=A2^{m+k}+B2^m+C2^k+D,\qquad
 (A,B,C,D)=(c_{22},c_{21}-2A,c_{12}-2A,c_{11}-2c_{21}-2c_{12}+4A).
$$
Both slopes with a fixed positive other width are positive. Cutting a Hamilton cycle of the strong quotient at either port gives a word visiting it twice and the other port once.

Put $K=\max(1,A,|B|,|C|,|D-e|)$, take the least positive $t$ with $2^t>12K$, and put $L=K(2^{t-1}+1)$. It suffices to take any $p$ with
$$
 p\ge t+1,\qquad 2^p>4L+2,\qquad p\ge2\lfloor\log_2K\rfloor+4.
$$

First suppose a width, say $m$, is less than $t$. Then $H-e=U2^k+V$, with $1\le U\le L, |V|\le L$. Equality gives
$$
 2^{2p-1-k}=U+\frac{2^{p-1}+V}{2^k}.
$$
The fraction is strictly positive because $2^{p-1}>L$. Also $2^k>2^{p-1}(2^p-2)/U$, so the fraction is less than $2U/(2^p-2)<1/2$. No integer-exponent power of two lies strictly between $U$ and $U+1/2$, even when negative exponents are allowed.

Now $m,k\ge t$. Reduction modulo $2^t$ forces $D=e$. If $A=0$, then $B,C\ge0$ and the binary digit sum of $B2^m+C2^k$ is at most the sum of the binary digit sums of $B,C$, at most $2(1+\lfloor\log_2 K\rfloor)$. But $E_p$ has $p$ ones. If $A>0$, put $\delta=B2^{-k}+C2^{-m}$, so $|\delta|<1/4$. The equality now gives
$$
 2^{2p-1-m-k}=A+\delta+\frac{A+\delta}{2^p-1}.
$$
The last term is positive and less than $1/4$. Thus the power of two lies strictly within $1/2$ of the positive integer $A$, forcing equality with $A$ itself. This is the step that turns separately selected widths into an exact fixed binary scale.

With $u=2^{p-1}$, write $2^m=2^i u, 2^k=2^j u$. We have
$$
 A=2^a,\qquad i+j=s=1-a,\qquad B2^i+C2^j=-1,\qquad D=e.
$$
When the valuations of the two nonzero summands differ, the smaller valuation is zero. When they agree, equality of valuations already determines $i$. Hence the candidate first shifts are
$$
 i=-v_2(B),\qquad i=s+v_2(C),\qquad i=\tfrac12(s+v_2(C)-v_2(B)).
$$
Omit expressions using a zero coefficient and discard nonintegral candidates. Substitution tests the remaining pairs exactly. If both coefficients vanish, there is no candidate. At most one pair survives. With $z=2^i$ the linear-scale equation is
$$
 Bz^2+z+C2^s=0.
$$
If $B=0$ this equation has at most one root. For $B\ne0$, two distinct positive dyadic-power roots $2^i,2^{i+h}$, $h>0$, would have sum $-1/B$. Thus $-1=B2^i(1+2^h)$, impossible modulo the odd divisor $1+2^h$ in $\mathbb Z[1/2]$. This proves uniqueness for the ordered pair of ports. Conversely a successful pair has the requested polynomial identity for every positive pair of widths on its ray. This also proves the infinite-exponent assertion for an arbitrary sign or integer offset.

## Detecting two visits without generating wider tournaments

For a port $v$, each directed triangle through it consists of $v\to o\to i\to v$, with $o,i$ in its out- and in-neighborhoods. Intersect the pairs $\{o,i\}$ over these triangles. For the two ports $v,w$ in a strong quotient, **every path has at most two visits to $v$ at all widths if and only if this intersection contains an ordinary vertex** meaning a vertex different from $v,w$. Apply the test in both directions. This specializes the single-port matching theorem to a second, independently repeatable vertex. [Matching and retention theorem](../chanpak_reader_retention_20261115/MATCHING_AND_RETENTION.md)

Indeed each return from $v$ to $v$ contains a transition from an out-neighbor to an in-neighbor, so an ordinary vertex in the intersection is consumed on any return. For the converse duplicate $w$ as a transitive pair. The bipartite graph of arcs $o\to i$ now has no vertex meeting all its edges. An ordinary common vertex would have been common before duplication. A common copy of $w$ is impossible since the two copies have identical neighbors across this bipartition. A bipartite graph with edges and no common endpoint has two disjoint edges because an intersecting bipartite edge family is a star. The resulting two triangles give a path with three copies of $v$ separated by two excursions. Insert all unused vertices by the usual tournament path insertion, keeping the original subsequence. The three separated visits persist.

## Exhaustive quotient search

The counts below use all strong isomorphism classes generated by nauty. Each pair denotes an unordered pair of vertex labels in one emitted representative, before identifying automorphic pairs. The dynamic program appends a legal next symbol to an ordinary subset, two visit counters, and a last label. The candidate $K$ for either sign is bounded by $\max(A,|B|,|C|,|D|+1)$. At order nine its maximum is 825, giving threshold 25. At order ten it is 3,735, giving threshold 29. For smaller exponents solve at each fixed smaller width $m$ by testing whether $(E_p+e-D-B2^m)/(A2^m+C)$ is a positive power of two. The cycle cuts give $H\ge2^m-2,2^k-2$, so each positive width is at most $2p$. Also $\min(m,k)\le v_2(E_p+e-D)$ whenever the difference is nonzero. Both orderings of the ports are tested. [Nauty and Traces](https://pallini.di.uniroma1.it/) · [Enumeration and exact integer replay](port_rays.cpp)

| Quotient vertices | Strong classes | Classes with an eligible pair | Eligible pairs | Lower exponents | Upper exponents |
|---:|---:|---:|---:|:---|:---|
|3|1|1|3|2|none|
|4|1|1|5|2|none|
|5|6|4|22|3|none|
|6|35|23|105|3, 4|3|
|7|353|154|573|3–5|3–5|
|8|6,008|1,511|4,026|4–6, 8|4–7|
|9|178,133|22,184|41,994|4–9|4–10, 12|
|10|9,355,949|511,366|743,104|5–12|5–13|

The last upper hit reveals how a width can correct a constant at one exponent. Its polynomial is
$$
 H(m,k)=1024\cdot2^{m+k}-207\cdot2^m+18\cdot2^k-63.
$$
At widths $m=p-7,k=p-4$ its excess over $E_p+1$ equals $2^{p-7}-64$. At $p=13$ the entire secondary scale is supplied by $63(2^6+1)=2^{12}-1$. A transfer to $p=61$ with the same difference-of-squares identity would require a construction producing $2^{30}-1$ in the role of 63. The theorem proves this cannot be achieved simply by choosing two new widths on any of these quotients.

For uniform offsets with $|e|\le100$, the number of marked-pair and shift solutions at order nine is two apiece for $-5,-9,-15$. The order-ten totals are eight apiece for $-5,-9,-15$, and two for $-41$. All smaller orders give none. The $-5,-9,-15$ templates each save one vertex over the previous ten-vertex constructions. The first source strings below give the upper triangle in row-major order from vertex zero, with 1 for an arc pointing from the lower label to the higher one.

| $e$ | Upper triangle | Ports | $(i,j)$ in widths $p-1+i,p-1+j$ | $(A,B,C)$ |
|---:|:---|:---|:---|:---|
|−5|111110111111111111011101011111111101|1, 8|−2, −3|64, −23, 38|
|−9|110111111101111110111111111110101111|0, 8|−4, −2|128, −72, 14|
|−15|110111111101101111010111111111111111|0, 8|−3, −3|128, −47, 39|
|−41|101011001111111111111111110111111111111111110|2, 8|−3, −3|128, −85, 77|

From the repository root the following reproduces the entire table and prints the first source, ports, and widths for each attained exponent and sign. The nine- and ten-vertex results concern the specified visit condition. Larger quotients, a third moving module, and a third visit fall outside this search. The exact Chan–Pak spectrum and an exponent-61 one-arc straddle remain open.

~~~sh
c++ -O3 -std=c++17 05-knowledge/results/chanpak_square25_dilation_20261123/port_rays.cpp -o /tmp/port_rays
for n in 3 4 5 6 7 8 9 10; do gentourng -cq "$n" | /tmp/port_rays "$n"; done
~~~
