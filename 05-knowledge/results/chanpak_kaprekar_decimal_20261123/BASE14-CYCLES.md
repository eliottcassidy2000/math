# Base 14: all three width-two arc-walk optima

The base-14 gap map is \(g\mapsto|2g-15|\). Its positive gap orbits are \((1,13,11,7)\), \((3,9)\) and \((5)\), and the literal Hamiltonian-path heights are the gaps times 13.

| Prescribed orbit | Minimum reversals | Minimum order at that length |
|---|---:|---:|
| \(13\to169\to143\to91\to13\) | 8 | 8 |
| \(39\to117\to39\) | 2 | 9 |
| \(65\to65\) | 0 | 7 |

Both columns hold whether the final tournament must be identical or only isomorphic. For the long orbit the edit bound holds even if the last tournament merely has height 13.

## Literal witnesses

Rows below are outgoing bitmasks, with vertices numbered from zero.

- **Long orbit:** \((234,244,249,242,225,192,128,0)\), and arcs \((3,7),(4,6),(1,3),(1,4),(1,3),(4,6),(1,4),(3,7)\) in that order. The exact trace is \(13\to85\to169\to99\to143\to203\to91\to85\to13\). Every marked arc is reversed twice. It detours through the base-12 phase 99; the [base-12 optimum](BASE12-CYCLE.md) in turn detoured through the base-14 phase 169.
- **Two-phase orbit:** \((430,508,504,496,481,448,385,256,0)\), flip \((6,8)\) to reach 117 and undo it. The cycle must visit two distinct heights, so two is optimal. The padded order-eight height-39 bank contains six shapes and has no hit among 168 marked arcs. At order nine the twelve shapes give two hits. The bank is obtained from the two strong height-39 classes of order seven or from a triangle and the unique order-five height-13 core, in both orders.
- **Fixed orbit:** \((18,20,9,19,36,79,31)\) has height 65. No Hamilton-framed orientation among the \(2^{10}\) order-six frames has height 65. Padding rules out all lower orders.

## Three forced detours for the long orbit

The unique nontrivial SCC of a height-13 source is the strong order-five core, e.g. \((10,20,25,18,1)\). The gains to 169 and 91 have valuations two and one; neither ratio is exceptional in the [singleton budget](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget). The four order-eight padded shapes give zero \(13\to169\) hits in 112 arcs; the three order-seven shapes give zero \(13\to91\) hits in 63 arcs. The first and last segments each require two edits.

The chain \(169\to143\to91\) cannot consist of two single-arc steps. Read it backwards. In a [strictly increasing walk](BASE16-CYCLE.md#optimal-even-up-to-relabelling) terminal singleton SCCs are untouched and removable. Height 169 has at most ten nontrivial vertices (one strong core or two order-five height-13 cores), by [Busch's exact minima](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). Height 91 has just one strong core of order at most nine. A split would require height 7, which is excluded by the order-four maximum of 5 and the strong order-five minimum of 9. The anchored strong-91 bank has 15 order-seven classes, 12 order-eight classes and two order-nine classes (182 order-nine frames before deduplication). Pad to ten. From its 100 shapes there are ten arcs to 143, zero then to 169. This forces a third, disjoint extra edit.

Four transitions plus three extras leave only two seven-edit allocations:
\[
A=(2,1,2,2),\qquad B=(2,2,1,2).
\]

## Schedule A forces a split 143, then fails

Retain the union of the initial 13-support and the nontrivial support at 169 in the first two steps. It has at most \(5+10=15\) vertices. Deleting a vertex singleton at both points preserves both heights. Neither first-segment arc can disappear in the restriction, or else it would imply a direct \(13\to169\) arc. The descent \(169\to143\) must be internal, so the third arc survives as well. Pad to order 15.

From the eleven starting shapes, the two-edit peaks of height 169 are 16 strong order-eight, 28 strong order-nine, and 42 split into two order-five height-13 cores, with multiplicities. Only the split peaks admit the next 143, giving 168 arrivals, all split into heights 11 and 13.

Now \((11,13)\) cannot reach 91 in two arcs. Both arcs inside the original cores cannot produce the lone strong 91-core, because the cores stay ordinal and each has order five. If neither arc changes an induced core, the height remains at least \(11\cdot13\): [Busch's path-merging lemma](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) preserves the path orders of two disjoint induced blocks, yielding the product as a lower bound. Thus exactly one arc is internal to an original core. A single edit of the unique height-11 or height-13 core always leaves at least four nonsingleton vertices (eight edits leave five, two leave four in each). The intact other core supplies five, so the old and new supports overlap in at least nine vertices. Their union fits on \(10+9-9=10\) vertices. Deleting other singleton pads cannot erase an arc, because a single arc between the split 143 and 91 is impossible: its reverse, from the sole 91-core, must either leave a factor 91 or put that core into one larger SCC, not split it into two factors. The exact order-ten scan of the two split source orders gives zero hits in 1980 unordered arc pairs. The more conservative order-eleven scan has six padded shapes and zero in 8910 pairs.

## Schedule B has only four possible 143 cores

Here \(143\to91\) must be a single reversal. The standalone increasing \(91\to143\) also contracts to order ten: either 143 is strong of order at most ten, or it splits as \(11\cdot13\), each factor on an order-five core. The ten hits above are all internal flips preserving the core order (six at order eight, four at nine). They give exactly four strong 143 classes:
\[
\begin{split}
&(2,4,201,147,39,207,155,19),\\
&(2,196,9,19,231,207,141,13),\\
&(322,4,393,499,487,455,390,259,2),\\
&(450,4,393,499,487,391,422,258,2).
\end{split}
\]
A witness would therefore require a reverse prefix \(143\to169\to13\), with two arcs per segment, starting from one of these cores and pads.

The union of the 143 and 169 supports has at most \(9+10-6=13\) vertices. In fact an internal flip of any of the four 143 cores leaves at least six nonsingleton vertices. If neither of the two arcs is internal to the original 143 core, the whole core survives; if exactly one is internal, its induced support persists through the other arc. If both are internal, all the 169 nontrivial support lies within the original 143 vertices, and it must have at least six vertices (order at most five has maximum 15). When reading \(13\to169\) backwards in this argument, two arcs inside the order-five 13-core cannot produce 169; zero or one internal arcs retain at least four of its original five vertices in the 169 support. Therefore the final 13 adds at most one vertex to the union, giving an **order-14 cap**. At these checkpoints any other singleton vertex can be deleted without losing either arc of either segment: \(169\leftrightarrow13\) is not a single edge, and a shorter \(143\leftrightarrow169\) would give the forbidden middle chain using the original internal 91 edge inside the retained 143 core.

Pad to 14. From 26 padded 143 sources, the 106470 unordered first pairs give 80 distinct labelled height-169 peaks. Testing all 180540 second pairs that internally edit at least one of the original peak SCCs gives zero count 13. Pairs with no internal arc preserve the product of induced SCC counts as a lower bound. A separate check at the unrefined order-18 cap gives 42 starts, 488376 first pairs, 128 distinct peaks, and zero hits among 531880 permitted next pairs. Thus neither seven-edit schedule works even with no closure constraint.


The same finite test can be *contracted to order eleven* by keeping the rooted information that support size alone discards:

- For a strong 169 peak of order at least eight, the final 13 cannot acquire a new vertex. If it did, only one of the two arcs could edit the original 169 core. That induced reading would have support at most four (it must fit inside the new five-vertex core with a vertex outside the old core), and therefore height 1, 3, or 5. The reverse one-arc move to 169 is impossible from 1, and the singleton budget bounds its strong order by five from 3 and seven from 5. For a split 169 peak with two 13 cores, both arcs instead have to act inside one core (use the same induced-product and exact ordinal-primality arguments as above). Thus the final 13 can add a vertex only when peak support is at most seven.
- One internal edit leaving exactly six nonsingleton vertices in any listed 143 core always leaves a strong height-23 block and singleton pads. Commute that internal edit first. The remaining \(23\to169\) gain 146 has valuation one and nonexceptional ratio. Since \(23\nmid169\), the 23 core must merge, and at most two singleton vertices can enter the merged support. Thus the peak has at most eight nonsingleton vertices.
- The remaining possibly large case is an order-nine 143 core, one internal edit leaving seven vertices, and a peak on ten nontrivial vertices. In both ninth-order cores the edit leaves a strong order-seven height-29 block and two singleton pads, both on the same side. The cross move to a strong order-ten 169 has gain 140. A pair of singleton endpoints would give the intermediate height-29 core and at least one more intermediate singleton, of total valuation at least three, impossible. The rooted coordinate at the core endpoint must therefore be 35, since the two remaining intermediate singleton pads supply a factor four. In the first order-nine core the rooted vectors on the 29 block contain 35 only as a starting coordinate, giving a merge to the right toward its two old singleton sinks. In the other they contain 35 only as an ending coordinate, giving a merge to the left toward its two old singleton sources. In either case the merge absorbs both old pads and hence retains all nine peak-source vertices.

With these refinements the union is at most eleven: overlap six allows peak mass at most eight; overlap seven in an order-nine source allows peak mass at most nine, unless the actual overlap is nine; and overlap at least eight already gives the bound. A possible newly added final vertex at peak mass at most seven also satisfies it. The order-eleven replay has 14 sources, 20790 first pairs, 44 distinct peaks and zero hits among 51120 permitted next pairs.

At eight edits the three segments after the first cost at least five in total. At order seven, none of the three padded 13 sources attains 169 in two or three arcs (630 pairs and 3990 triples). Padding excludes all smaller hosts. This proves the long orbit's exact minimum order.

[Core and witness replay](kaprekar-base14-cores.js) and [schedule replay](kaprekar-base14-schedules.js) use exact subset path counts and ordinal factorization; the largest padded shapes need subset tables only within induced SCCs.

## A common infinite lift of all three orbits

The [transport rule](BASE16-SECOND-CYCLE.md#both-base-16-orbits-lift-together) takes attained factors \(d=13x+1,\ e=15x+1\), and lifts base 14 to \(B=195x+14\) by the ordinal multiplier \(de=d(B-1)/13\). For \(k\ge6\), take \(x=2^k-2\). The signed shell \(G2^w+D-L\) gives both factors:

| Factor | Core rows | Tag | \((D,G,L)\) | Filler width |
|---|---|---:|---|---:|
| \(13\cdot2^k-25\) | \((6,4,8,19,39,15)\) | 37 | \((15,104,40)\) | \(k-3\) |
| \(15\cdot2^k-29\) | \((34,28,9,17,37,78,31)\) | 56 | \((89,480,118)\) | \(k-5\) |

Ordinally join the two fixed shells to each of the three witnesses. This transports all three cycles to each base \(195\cdot2^k-376\), starting at 12104, at the same edit lengths. The long, two-phase and fixed constructions use \(2k+15,\ 2k+16,\ 2k+14\) vertices respectively. These are lifted constructions, not minimum claims in the larger bases.
