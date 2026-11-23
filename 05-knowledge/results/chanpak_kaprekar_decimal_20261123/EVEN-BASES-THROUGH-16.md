# All even radices through 16, including the denominator-nine gap

The smallest even base admitting exact uniform dilation of the gap cycles from **every even base 2 through 16** is **765764**, and every one of its nine distinct imported cycles has a literal tournament walk. The base has a 36-vertex shared layout; no optimality of this layout or the lifted edit lengths is claimed. The two small results below explain why this is not just an arithmetic least common multiple.

## Base 6: labels cost one edit, but save one vertex at the optimum

The sole nonzero width-two orbit is \(5\to25\to15\to5\).

- **Labelled: four edits, minimum order six.** Start with outgoing masks \((54,60,56,49,32,0)\), and reverse \((0,4),(3,5),(0,4),(3,5)\). The exact trace is \(5\to11\to25\to15\to5\). A return on identical labels uses an even number of reversals, with at least three distinct-height transitions. Height 25 exceeds the order-five maximum 15.
- **Up to relabelling: three edits, minimum order seven at this length.** Start with \((118,124,120,113,96,64,0)\) and reverse \((3,6),(5,6),(3,6)\). The trace is \(5\to25\to15\to5\), with endpoint isomorphism \((5\ 6)\). Three is the number of required transitions. A height-5 host has the unique strong four-vertex core plus ordinal pads. At order six, all three padded starts and their 45 arcs give no height-25 neighbor. Smaller hosts pad to six.

These statements minimize vertex order **at the stated edit optimum**. The labelled six-vertex witness is also of course valid when relabelling is allowed.

## Base 8: dilate the gaps, not the forbidden Hamiltonian counts

The two base-8 gap orbits are \((1,7,5)\) and \((3)\), whose ordinary output heights would be \(7,49,35\) and 21. There is no height 7: a split would give at least \(3\cdot3=9\), a strong order of at most four gives at most 5, and [Busch's exact lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) gives at least 9 from order five. A split of height 21 would require height 7. A strong 21 could only have order six (order five max 15, strong order seven min 25), and the 1024 Hamilton-framed order-six orientations give no height 21.

A different source packet uses **five times the gaps**, not the impossible seven times the gaps. On the seven-vertex start of the preceding section, reverse \((0,4),(4,6),(3,5),(4,6),(0,4),(3,5)\). It returns on identical labels with trace
\[
5\to11\to35\to55\to25\to15\to5.
\]
At target base 26 the gap dilation is \(d=(26+1)/9=3\). Ordinally multiply this packet by \(15=(26-1)d/5\), using a triangle and the height-5 core. The actual trace is
\[
75\to165\to525\to825\to375\to225\to75.
\]
The detour at 225 is the lifted fixed phase, so both arithmetic orbit types are realized. This is the earliest larger even base for uniform dilation of the denominator-nine gaps, since \(B+1\) must be an odd multiple of nine. The six-edit walk is an **upper bound for the unrestricted base-26 problem**.

The normalized packet, with required heights \(5,35,25,5\), does have exact labelled minimum six. The gain for a proposed \(5\to35\) edge has valuation one and nonexceptional ratio, so the [singleton bound](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget) reduces it to order six; no edge exists in the 45 tests there. If the whole walk uses at most four edits, the first segment must take two and the other two must take one each. Retain the union of the active supports at heights 5 and 35, at most \(4+7=11\) vertices; 35 has one strong core since 7 is missing. Any removed vertex is singleton on both sides, and neither of the two first flips can touch it, since that would leave a forbidden one-edge \(5\to35\). The strict descending suffix uses only internal arcs and cannot bring back any discarded singleton. At order eleven, the eight padded sources give 50 two-edit peaks of height 35 in 11880 pair tests, and zero 25 neighbors among their 2750 subsequent arcs. Thus at least five edits are needed even without a return isomorphism. Labelled parity forces six, attaining the given normalized upper bound. This labelled lower bound also applies to an ordinally scaled copy with edits confined to the normalized block, not to every way of realizing the lifted heights.


Up to relabelling the normalized optimum is **five edits**. On that same seven-vertex start, reverse \((0,4),(4,6),(0,2),(4,5),(0,4)\). The trace is \(5\to11\to35\to29\to25\to5\), with endpoint map \((1,2,3,0,5,6,4)\) from initial to final labels. It scales at base 26 to \(75\to165\to525\to435\to375\to75\). This five-edit packet does not visit the fixed phase; the six-edit labelled packet above does. Seven vertices are minimal for both normalized optima since no order-six Hamilton-framed orientation has height 35.

## The common radix and the missing multiplier

Put \(M=B+1\), \(h=M-2\), and \(L=109395\) as in [the earlier seven-cycle construction](COMMON-BASE.md). Every even source \(b\in\{2,4,\ldots,16\}\) has a gap-one phase. Preserving its whole orbit by uniform gap dilation forces \(b+1\mid M\). Thus even just the arithmetic gap requirement gives
\[
\operatorname{lcm}(3,5,7,9,11,13,15,17)=765765=7L\mid M.
\]
At \(M=765765\) there are nine distinct imported cycles, consisting of all sign-quotient doubling cycles in the eight reduced denominators \(3,5,7,9,11,13,15,17\). Denominators 3 and 5 already occur in the earlier base-14 packet bank. Denominator 17 contributes two disjoint cycles; each of the other seven contributes one.

The earlier sufficient factorization would ask for a tournament of count 7. That is impossible but unnecessary: attain the whole fixed factor
\[
F=7h=7\cdot765763=5360341.
\]
The following small certificates are all **transitive spines with two apices**. The table gives the spine length, the outgoing masks of the two apices into spine coordinates \(0,\ldots,n-1\), and their mutual arc. This specifies every orientation.

| Exact height | Spine length \(n\) | \(a\) outgoing mask | \(b\) outgoing mask | Mutual arc |
|---:|---:|---:|---:|---|
| 5360341 | 16 | 9064 | 24227 | \(b\to a\) |
| 765763 | 13 | 829 | 533 | \(b\to a\) |
| 21879 | 11 | 1096 | 275 | \(a\to b\) |

There is also an independent signed-shell witness of height \(F\), with core rows \((214,76,72,49,102,199,136,286,127)\), tag mask 460, and a nine-vertex transitive filler. Its exact coefficients are \(D=881,G=10478,L_{\rm tag}=5276\), so the [shell formula](BASE16-CYCLE.md#transporting-a-literal-orbit-to-larger-bases) gives \(10478\cdot512+881-5276=5360341\). The replay checks these certificates with full subset DP.

**A small automaton for exact readings.** For any fixed order \(x,y\) of the apices, a path breaks the transitive spine into three slots before, between and after them. Within a slot the order is forced. Scan the spine in its increasing orientation and assign each vertex to a slot. In the first slot track whether it is empty, or its last vertex can or cannot lead into \(x\). In the middle slot require \(x\) to beat its first vertex, and track the analogous last-vertex guard for \(y\). An empty middle slot is allowed at the end exactly when \(x\to y\) exists. In the last slot track just occupancy, requiring \(y\) to beat its first vertex. This gives 18 states per apex order and evaluates the **literal integer count**, not its residues. The two apex orders are disjoint and exhaust all paths. More generally, with \(k\) apices an order of the apices needs at most \(2\cdot 3^k\) states for the \(k+1\) spine slots; sum over all apex orders.

## A shared program at this radix

The prior seven packets use their original selector factors \(1105,765,561,429\) times \(F\). The normalized base-8 packet takes multiplier
\[
\frac{h(M/9)}5=2431F.
\]
Its factor 2431 is reached on the **same nine selector vertices**. Start from \((444,401,138,322,76,414,423,24,148)\), and on this branch reverse \((1,2),(1,4),(0,2),(5,7),(5,8)\). Its selector readings are \(1105,765,955,1517,2265,2431\). This branch applies also to every earlier common base \(B=Ly-1\) with an attained fixed factor \(y(B-1)\), so all six previous certified bases pick up the denominator-nine cycle with no added vertices.

For base 6 the required multiplier is \(K_6=Mh/35=21879\cdot765763\), not a multiple of \(F\). The 11-vertex certificate for 21879 in the table below is smaller than its two-apex certificate. The labelled base-6 packet scales to
\[
83770643385\to184295415447\to418853216925\to251311930155\to83770643385.
\]

**One fixed multiplier core, 36 vertices total.** Keep the same 15-vertex \(h=765763\) core for all nine packets. The remaining factor is \(7C_b\) for the earlier and denominator-nine packets, or 21879 for the base-6 packet. A programmable 11-vertex slot carries these six states (pad the ten-vertex witnesses by an ordinal singleton).

| Exact slot height | Outgoing masks |
|---:|---|
| 7735 | \((670,200,322,468,166,847,17,356,595,222)\) |
| 5355 | \((584,509,913,580,649,733,276,841,569,66)\) |
| 3927 | \((160,245,625,199,681,328,17,100,735,235)\) |
| 3003 | \((602,564,977,246,128,661,562,67,763,152)\) |
| 17017 | \((1494,2028,48,997,682,1921,564,580,724,5,988)\) |
| 21879 | \((1024,1017,539,1377,265,789,821,381,1541,1177,246)\) |

Together with the ten-vertex moving packet, the layout has order \(15+11+10=36\). Change only differing internal arcs when switching programs; do not change ordinal cross arcs. No short program is claimed for this 11-vertex slot. If instead the original nine-vertex selector is wanted, a 37-vertex alternative has a 27-vertex multiplier block holding either \(F\) and that nine-vertex selector, or the 15+11 base-6 factor and one pad.

The earlier 45-vertex layout also accommodates this new radix along with the six previously documented bases: keep ten moving vertices and allow the 35-vertex multiplier portion to be fully reprogrammed. No programming edits are counted as part of a replayed cycle.

[Runnable exact replay](kaprekar-even-radices.js) verifies the two closure conventions, the small order bounds, the normalized four-edit obstruction, both factor certificates, the two-apex automaton, and the exact scaled traces.
