# Two visits per module at the exponent-61 perfect number

Put $E=2^{60}(2^{61}-1)$. Mark a set of distinct sites in a strong tournament and replace each marked site by a transitive chain of independently chosen positive width. Between two blocks keep the source orientation. A *visit* is a maximal consecutive run of one marked chain in a Hamilton path. Call the marked source *two-visit* if every choice of widths allows at most two visits to each chain.

**PROVED.** A single strong ten-vertex source realizes a product of three factors by independently replacing four sites with chains. Let its ordinary order be $0,\ldots,9$, with all arcs directed forward except $3\to0$, $5\to3$, $8\to5$ and $9\to7$. Set $A(X,Y)=4XY+3X-3Y-1$ and use $X,Y,Z,W$ for the four powers of two. The two following profiles are exact for all positive widths.

| Marked sites, in order | Height |
|---|---|
| $0,4,6,9$ | $A(X,Y)(1+Z)(1+W)$ |
| $0,4,6,7$ | $A(X,Y)(1-Z+W+2ZW)$ |

The ten-site source has a directed Hamilton cycle $0,1,2,4,6,9,7,8,5,3,0$ and every inflation is strong. The ordinary triangle guards give at most two visits to each mark. Applying the recurrence below gives coefficient vector $(-1,3,-3,4)$ in $1,X,Y,XY$ for $A$. The four consecutive blocks for $1,Z,W,ZW$ are respectively $A,A,A,A$ and $A,-A,A,2A$. This proves both identities with an integer recurrence on a ten-site source.

For every $p\equiv1\pmod4$ with $p\ge9$, set $a=(p-5)/4$ and $q=(p-1)/2$. With the first marking put widths $(a,a,q,p)$ or $(a,a,p,q)$. The diagonal cancels the signed linear terms, so the common height becomes
\[
 (2^q-1)(2^q+1)(2^p+1)
 =(2^{p-1}-1)(2^p+1)
 =2^{p-1}(2^p-1)-1 .
\]
Both orders are $2p+3$. The construction uses the actual powers in the integer factorization and works whether or not $2^p-1$ is prime. The second marking, with widths $((p-3)/2,(p-3)/2,(p-1)/2,(p-1)/2)$, gives the same height on $2p+2$ vertices for every odd $p\ge5$. At $p=61$ these half-width hosts are the earlier reflected lower and its converse. The quarter-width hosts and their converses yield four additional isomorphism classes, all of order 125.

**FINITE-EXACT.** Allow any number of independent two-visit marks and any strong source of order at most ten. Exactly eight isomorphism classes in this family realize $E-1$, four on 124 vertices and four on 125. None realizes $E+1$. Starting at any of the eight lower classes, and making arbitrarily many height-preserving flips, reaches ten classes in total. Every single physical arc reversal from every one of the ten misses $E+1$. Larger sources and placements admitting three visits remain open.

**FINITE-EXACT.** Among all two-visit sources with exactly three marked sites and at most ten vertices, the only tournaments of height $E-1$, up to isomorphism, are the inherited 124-vertex strong lower and its converse. No such source realizes $E+1$. Furthermore, reversing **any** single physical arc of **any** lower in this family misses $E+1$, including reversals that destroy all of its marked modules. The earlier all-arc obstruction gives this last conclusion once the entire lower fiber has been identified. [Strong lower construction](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md) · [All-arc census for the inherited lower](../chanpak_e61_fibonacci_target_20261122/PLATEAU.md) · [Earlier two-port census through ten](../chanpak_square25_dilation_20261123/TWO-PORT-RIGIDITY.md)

| Source order | Strong source classes | Admissible three-site placements | Placements and widths giving the lower | Upper |
|---:|---:|---:|---:|---:|
| 4 | 1 | 2 | 0 | 0 |
| 5 | 6 | 15 | 0 | 0 |
| 6 | 35 | 84 | 0 | 0 |
| 7 | 353 | 458 | 0 | 0 |
| 8 | 6,008 | 2,674 | 2 | 0 |
| 9 | 178,133 | 18,649 | 12 | 0 |
| 10 | 9,355,949 | 180,852 | 42 | 0 |

A placement is an unordered subset of sites on one labelled representative per source isomorphism class. Widths belong to the individual sites. At each of orders 8, 9 and 10, exactly half of the lower placements yield each of the two converse classes. The directed triangle with all three sites marked fails the two-visit condition.

## Which labelled placements suffice?

**PROVED.** For any number of marked sites, the two-visit condition holds exactly when, for each mark $v$, the directed triangles through $v$ have a common *ordinary* vertex apart from $v$. In a strong source every site is on a directed triangle.

Indeed, a return to $v$ traverses an edge from an outneighbor of $v$ to an inneighbor of $v$; its endpoints complete a directed triangle. Two successive returns have disjoint intermediate vertices in the inflation. A common ordinary vertex therefore allows at most one return. Conversely, double each other mark. If the bipartite graph of return edges in this expansion had no two independent edges, its edges would all share one endpoint. Such an endpoint must be an ordinary vertex, since every marked site has an interchangeable second copy. With no ordinary guard there are therefore two disjoint return edges. Keep their four endpoints aside, contract the selected chain to the vertex $v$, and take a Hamilton path of the remaining tournament. At its occurrence of $v$ insert the two disjoint triangular detours, using three copies of $v$. This gives a spanning path with three visits. The same test handles one, two, or any number of sites, and requires only the source triangle incidences. [Earlier shared-endpoint test](../chanpak_passive_pubsub_20261117/reader/PROOF.md)

For an admissible set $P$, let $c_S$ count spanning directed words on the quotient in which exactly the marks in $S\subseteq P$ occur twice. Every ordinary site and every other mark occurs once. A two-run chain of width $w$ contributes $2^w-2$ ordered nonempty splits of its ranks, and a one-run chain contributes one. Thus the full height profile is the multilinear polynomial
\[
 F(X_1,\ldots,X_k)=\sum_{S\subseteq P} c_S\prod_{j\in S}(X_j-2),\qquad X_j=2^{w_j}.
\]
The recurrence appends an unused ordinary site or an available mark after an outgoing arc, with a binary mask for ordinary sites and ternary visit counters for marks. Each coefficient is computed by integer addition. Cutting a Hamilton cycle of the strong quotient at a mark gives a spanning word visiting that mark twice and each other site once. Consequently the profile strictly increases with any chain width and satisfies $F\ge 2^{w_i}-2$. The possible widths form a finite, explicitly bounded set even though the problem imposes no size bound on the inflated tournament. [Word recurrence](seam_profiles.cpp)

## Exact inversion and class identification

Sort the unknown widths temporarily. If the least width is $m$, then $2^m$ divides $N-F(0,\ldots,0)$ for the desired exact count $N$. Substituting that width leaves the identical test for the next width with the new constant term. A common divisor of the remaining nonconstant monomial coefficients also divides the residual. Whenever the diagonal value for all remaining widths exceeds $N$, positivity of the visit expansion rules out that width and all greater choices. In the last variable, exact division leaves a power-of-two test. The search includes ties. Variables that occur together in precisely the same monomials enter solely as one product, so their widths can be merged and then split over all positive compositions. None of these divisibility tests converts a congruence into a solution without evaluating the final equality. [Width solver](seam_widths.py)

The source polynomial determines the attainable heights but loses the physical location of the flipped arc. Retain **every** source placement carrying a winning coefficient vector. In each winning weighted tournament, contract a pair of sites that have identical relations with every third site, adding the two chain widths. The pair forms a larger transitive module. Maximal transitive modules partition the vertices: overlapping transitive modules have transitive union since the three parts (left difference, intersection, right difference) are linearly ordered, and any nontrivial transitive module contains an adjacent homogeneous pair. Hence the weighted irreducible quotient is an isomorphism invariant. Once reduced, sort by chain width and minimize the triangle word over permutations within each equal-width block. The two resulting canonical encodings agree with the source having outgoing masks $(2,8,3,37,15,151,63,95)$, widths $(1,29,1,1,29,1,61,1)$, and with its converse, respectively. Taking converses reverses each directed path and commutes with arc reversal. [Class identification](seam_hosts.py)

## Exhausting any number of marked sites

Two source twins agree on all third vertices and hence form a transitive pair. If at least one twin is marked, contract both to one marked site and add their positive widths (an ordinary vertex contributes width one). Every resulting tournament remains exactly the same. The smaller marked source still satisfies the ordinary-guard criterion. No triangle uses both twins; a triangle using one can substitute the other. Thus an ordinary twin removed by contraction cannot have served as a guard for a third site, and the merged site inherits the guard of its marked summand. Induction permits us to enumerate only sources in which every marked site is free of twins. This reduction leaves two ordinary twins available as separate singletons.

The empty marked set has height at most $10!$. For a single mark the profile is $a2^w+d$ with $1\le a\le11!$ and $|d|\le 10!+2\cdot11!$. If $w\le60$ even the largest such value is below $E-1$. If $w\ge61$, either target requires $d\equiv2^{60}\pm1\pmod{2^{61}}$, outside that bound. Two marks were already excluded at this exponent by the [two-port rigidity theorem](../chanpak_square25_dilation_20261123/TWO-PORT-RIGIDITY.md). All sites cannot be marked because a strong source requires an ordinary guard. This leaves three through nine marks. For exactly three we used every placement in the table above; for four or more, induction uses the twin reduction. At order ten the numbers of eligible placements after reduction for four through nine marks are respectively $8023,1602,312,44,3,0$. The first two arities produced lower hits, and all six missed the upper. Replaying every winning coefficient at all orders yields exactly the eight classes described above.

## Which physical arcs survive the width calculation?

The all-arc recurrence splits a Hamilton path at the arc into two paths with pinned boundary ranks. For each ordered join of the two paths it assigns the transitive ranks to the prescribed visits by inclusion-exclusion. Subtracting the two orders gives the exact physical change. Breadth-first search retains each distinct canonical tournament reached by a neutral flip. For the three classes displayed, their converses give identical statistics. The half-width reference class and its neutral partner exhaust the four classes of order 124 with their converses. At order 125, quarter A has one partner outside the initial eight, and quarter B is isolated up to isomorphism. The partners of quarter A and its converse are distinct: contraction gives their neutral degree one in the graph of distinct classes, and their sole neighbors are distinct. [All-arc recurrence](../chanpak_e61_fibonacci_target_20261122/rail_arcs.py) · [Neutral-component replay](seam_quarter.py)

| Widths and state | Physical arcs | Neutral flips | Distinct other neutral isotypes | Least positive change | Flips of $+2$ |
|---|---:|---:|---:|---:|---:|
| $(29,29,30,30)$ | 7,626 | 116 | 1 | 1,237,940,043,897,066,293,863,383,042 | 0 |
| Its neutral partner, the prior lower up to converse | 7,626 | 118 | 1 | 47,223,664,846,288,638,115,840 | 0 |
| Quarter A, $(14,14,30,61)$ | 7,750 | 117 | 1 | 35,184,372,056,064 | 0 |
| Quarter A's neutral partner | 7,750 | 115 | 1 | 1,729,382,255,299,657,728 | 0 |
| Quarter B, $(14,14,61,30)$ | 7,750 | 116 | 0 | 35,184,372,056,064 | 0 |

On the two quarter-width sources the smallest positive change is exactly $2^{15}(2^{30}-1)$. In A it is obtained by flipping the contact from site 8 to rank 15 of site 9. In B it reverses the contact between rank 45 of site 6 and site 8. Ranks start at zero. The neutral escape from A reverses site 7's contact with the center of the 61-chain at site 9.

## An exact rank derivative at every width

**PROVED.** In the first four-site marking, write $(x,y,z,w)$ for the widths and $(X,Y,Z,W)$ for their powers of two. Let $r$ be an interior rank in port 9 and $t=w-1-r$. Reversing its contact with ordinary site 8 changes the height by
\[
 \Delta_9=A(X,Y)\big((1+Z)2^r-2^t\big).
\]
At an interior rank $r$ of port 6 take $t=z-1-r$ instead. The change on that contact with site 8 is
\[
 \Delta_6=A(X,Y)\big((1+W)2^t-2^r\big).
\]
In particular these two types of contact cannot give a unit gap of two at any positive widths, since $A(X,Y)\ge15$. On the quarter-width lower the choice $r=(p-1)/4$ in port 9, or its reflected rank in port 6 in the swapped realization, cancels two equal powers exactly. The resulting positive change is
\[
 2^{(p-1)/4}\big(2^{(p-1)/2}-1\big)
\]
for every $p\equiv1\bmod4$, $p\ge9$.

To derive the two formulas, pin the distinguished rank at the joining end of one of two path fragments and enumerate the two fragment orders. A chain of width $u$ assigned to $k$ nonempty visits has $\sum_{s=1}^k(-1)^{k-s}\binom{k}{s}s^u$ assignments. If one visit must end at rank $r$, the other $k-1$ visits are free; after removing empty free visits by inclusion-exclusion, each lower rank has $s+1$ choices and each upper rank has $s$. This replaces the filling by $\sum_{s=1}^{k-1}(-1)^{k-1-s}\binom{k-1}{s}(s+1)^r s^t$. When the run must begin at the rank, exchange $r,t$. Interiority ensures both exponents are positive, allowing the zero base to be dropped. Collect the two signed orders in variables $2^r,2^t$. For port 9 the three surviving terms are $A2^r+ZA2^r-A2^t$. For port 6 the corresponding terms are $A2^t+WA2^t-A2^r$. This coefficient calculation uses the ten-site integer recurrence and holds for all the stated widths and ranks. [Symbolic coefficient replay](seam_derivative.py)

## Physical flips on the nearby offset families

**FINITE-EXACT.** The complete earlier catalogue of two-visit, two-port uniform families through ten sites with $|e|\le100$ has offsets $e=-5,-9,-15,-41$. At $p=61$ its 32 presentations produce eight actual tournaments, two mutually converse at each offset. Across their 58,564 physical arc reversals, none with nonzero change produces another count between $E-200$ and $E+200$. Every count-preserving reversal just exchanges consecutive ranks in one of the two chains. Thus every one of the eight is isolated under neutral moves up to isomorphism. [Uniform scalar families and census](../chanpak_square25_dilation_20261123/TWO-PORT-RIGIDITY.md) · [Labelled sources and all-arc replay](seam_offsets.py)

| Offset | Order | Arcs in either converse host | Neutral | Closest negative change | Closest positive change |
|---:|---:|---:|---:|---:|---:|
| $-5$ | 122 | 7,381 | 113 | $-1,873,497,444,986,126,344$ | $1,008,806,316,530,991,106$ |
| $-9$ | 121 | 7,260 | 112 | $-1,083,197,534,374,707,738,657,685,504$ | $792,633,534,417,207,294$ |
| $-15$ | 121 | 7,260 | 112 | $-1,702,167,554,017,397,880,536,432,640$ | $309,485,009,821,345,071,476,244,480$ |
| $-41$ | 122 | 7,381 | 112 | $-1,873,497,444,986,126,318$ | $1,008,806,316,530,991,098$ |

To exhaust the presentations, the published nine-site source for each of the first three offsets and its converse give the two nine-site placements. Split either marked site into a transitive twin pair, retaining either vertex as the mark and subtracting one from its inflated width. These eight ten-site placements give the same two physical hosts. They are distinct even when the two marks have the same color, by minimizing the labelled word over all permutations preserving the mark set. This accounts for all eight ten-site placements at each offset in the published census. For $-41$ the two marked ten-site sources are already distinct at the source level and exhaust that offset. Canonical comparison of the actual weighted twins also distinguishes each physical host from its converse. The recurrence that counts two pinned paths then gives the table.

**PROVED application to ordinal products.** All seven inherited triangular-offset constructions at exponent 61, with $e=-3,-15,-21,-45,-55,-91,-105$, also miss $[E-200,E+200]$ after any height-changing single flip. An internal change has absolute value at least twice the other component's Hamilton count, since heights are odd. Each of the two factors has height at least $2^{60}-15$. Across the cut of strong components of orders $a,b$, the gain is at least $2\lfloor ab/2\rfloor$. Indeed, a full path in the first component starting at the old tail and a full path in the second ending at the old head give one new odd cycle for every prefix/suffix pair of opposite parities. Every old odd-cycle packing survives and each new odd cycle contributes at least two. For the explicit signed factors each order exceeds 45, since their chain widths are at least $60-11$ and the pure $+1$ factor is the transitive chain with its extreme arc reversed. Thus even after the largest lower offset of 105, a cross reversal exceeds $E+200$. [Triangular factor representatives](../chanpak_absolute_liouville_leech_20261122/TRIANGULAR-OFFSETS.md) · [Cross-component gain bound](../chanpak_arc_fiber_grammar_20261122/README.md#the-only-unit-flips-crossing-a-component-boundary)

## Reproduction

From this directory, with nauty's `gentourng` and a C++17 compiler on the path, the following exhausts the three-site problem. The first pass deduplicates polynomials; the second deliberately keeps every winning placement.

~~~sh
work=$(mktemp -d)
c++ -O3 -std=c++17 seam_profiles.cpp -o "$work/profiles"
for n in 4 5 6 7 8 9 10
do
    gentourng -cq "$n" | "$work/profiles" "$n" 3 > "$work/p$n"
    python3 -B seam_widths.py 61 --winning-coefficients < "$work/p$n" > "$work/w$n"
    gentourng -cq "$n" | "$work/profiles" "$n" 3 "$work/w$n" > "$work/h$n"
done
cat "$work"/h* | python3 -B seam_hosts.py
~~~

The same binary with fourth argument <code>--reduced</code> enables twin reduction. For every $5\le n\le10$ and $4\le k<n$, replace the inner loop commands above with the following using distinct names for the chosen $n,k$. Concatenate the generated <code>hits-*</code> files for class identification. The selected-polynomial pass retains all placements.

~~~sh
gentourng -cq "$n" | "$work/profiles" "$n" "$k" --reduced > "$work/p-$n-$k"
python3 -B seam_widths.py 61 --winning-coefficients < "$work/p-$n-$k" > "$work/w-$n-$k"
gentourng -cq "$n" | "$work/profiles" "$n" "$k" --reduced "$work/w-$n-$k" > "$work/hits-$n-$k"
~~~

Running <code>python3 -B seam_quarter.py</code> reproduces the physical-arc table and the neutral components.
