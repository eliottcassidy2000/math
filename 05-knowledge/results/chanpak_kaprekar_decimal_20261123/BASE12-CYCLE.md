# Base 12: ten arc reversals, with a one-vertex gap between closure conventions

The width-two base-12 Kaprekar gap map is \(g\mapsto|2g-13|\). The unique positive gap orbit is \(1,11,9,5,3,7,1\), with required ordinary Hamiltonian-path counts
\[
11,\ 121,\ 99,\ 55,\ 33,\ 77,\ 11.
\]
The optimum is **ten literal reversals**. At this length the minimum host order is **nine up to relabelling** and **ten with identical labelled return**. The lower bound on length actually holds without any endpoint requirement beyond returning to height 11.

## Exact witnesses

Rows are outgoing masks and arcs use zero-based labels.

**Identical return (order 10).** Start at \((1022,972,984,976,994,974,896,768,512,0)\). Reverse
\[
(3,7),(3,9),(3,6),(1,3),(3,9),(1,3),(0,3),(3,7),(0,3),(3,6).
\]
The trace is \(11\to35\to121\to99\to169\to55\to33\to155\to77\to19\to11\). Each arc occurs twice and all labels return. The detour 169 is itself the square phase of the first base-14 orbit.

**Return up to relabelling (order 9).** Start at \((486,492,488,497,487,448,384,256,0)\). Reverse
\[
(2,6),(2,8),(2,5),(0,2),(2,8),(0,2),(1,4),(2,7),(1,5),(2,3).
\]
The trace is \(11\to35\to121\to99\to169\to55\to33\to41\to77\to119\to11\). The map from initial to final labels is \((4,0,5,1,3,6,7,2,8)\). It is not the identity. In fact it is unique: automorphisms of the initial tournament act freely on its eleven Hamiltonian paths, so their order divides both 11 and \(9!\).

## Three disjoint detours

1. The [first jump \(11\to121\)](FIRST-JUMP-SQUARES.md) requires two arcs at every order.
2. No arc joins 11 and 77. On the 11 side the positive gain 66 has valuation one and ratio seven. The [singleton budget](FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget) therefore keeps at most two pads over the unique order-five height-11 core. The three order-seven shapes give 63 arc tests and zero hits. The first \(11\to77\) at distance two has order seven: \((102,108,104,113,103,64,0)\) with arcs \((0,5),(1,6)\) gives \(11,29,77\); neither of the two order-six shapes reaches 77 in two arcs.
3. No single height-33 host has both a 55-neighbor and a 77-neighbor. Its possible nontrivial SCC heights are \((33)\) or \((3,11)\) in either order, with support at most eight. Gains 22 and 44 have valuations one and two, respectively, with neither exceptional ratio. Preserve the union of the two required pad sets, at most \(2+3\); the order-13 rooted scan of 103 shapes gives 194 arcs to 55 and 12 to 77, with no common source. The order-six and order-seven strong 33 classes number five and three respectively.

The three extras lie in disjoint parts of the six-transition walk. At nine edits there are only two possible segment schedules
\[
A=(2,1,1,1,2,2),\qquad B=(2,1,1,2,1,2).
\]

## Why nine edits fail at every order

**The common peak.** Both schedules contain consecutive one-arc descents \(121\to99\to55\). Read backwards. In a [strictly increasing walk](BASE16-CYCLE.md#optimal-even-up-to-relabelling), all terminal singleton SCCs are untouched and may be deleted. At 121 the terminal core is either strong of order at most nine, or two order-five height-11 cores. Thus pad to ten. The complete 22-shape height-55 bank consists of four strong order-seven cores or ordinal joins of heights 5 and 11 in either order, with singleton placements. Eight first arcs reach 99, and twelve following arcs reach 121; every peak is two order-five height-11 cores. There is also a divisibility explanation for the only possible split: the quotient \(99/55=9/5\) adds a factor nine to the edited terminal block, so the next arc must edit that same block because \(9\nmid121\). The factors outside it divide 11.

Both arcs of the first segment \(11\to121\) must act inside one peak core. Work backwards from the peak. When at most one arc is internal to the two cores, each induced reading is at least 3 and one remains 11. To justify the lower bound on the whole height, [Busch's path-merging lemma](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) preserves the relative orders of paths in both disjoint sets, so the product of their induced counts is a lower bound. One edit cannot make a height-11 core transitive, since a single edit out of a transitive tournament has height 1 or \(1+2^s\). If the two internal arcs act in different cores, the inter-core arcs are unchanged and the product is exact, with both factors exceeding one; it cannot equal prime 11. Thus the initial support is precisely the untouched peak core. On \(121\to99\), the active core drops to 9. The next step cannot leave 9 as a factor of 55, so that same core drops to 5. In schedule \(A\) the single next step changes it from 5 to 3, as 5 cannot remain a factor of 33.

**The support at 77 overlaps the peak in at least five vertices.** Height 77 has exactly one nontrivial SCC and it has order at most nine: a split would require a factor 7, which cannot be a tournament height (the strong order-four maximum is 5 and the order-five minimum is 9). In \(A\), the previous height 33 consists of peak-derived cores of heights 11 and 3. Forming the 77 core in two edits requires at least one cross-component edit: otherwise no SCC could grow and the order-five maximum is 15. Hence at most one edit changes arcs inside the original cores. If the 11-core is intact it contributes five vertices to the single 77-core. Otherwise an induced triangle in its edited copy survives and the original 3-core is untouched, giving at least six. For \(B\), the one-arc \(33\to77\) must start from a strong order-seven 33-core plus pads. This is checked at order nine by monotone contraction: four hits from the 29 padded strong-33 shapes, all using order-seven cores, and zero from six reducible shapes. Forming that strong 33 core in two edits from the peak-derived 11 and 5 cores again requires a cross-component edit (a strong 33 needs at least six vertices). An untouched 11-core gives five peak vertices, or an edited 11-core retains a triangle and the untouched 5-core supplies four. This overlap persists in the subsequent increasing step to 77.

**The final 11 core adds at most one outside vertex.** In the last two arcs, at least one must edit the induced 77 core, or its count would stay 77 and bound the entire final height from below. If both do so, the ordinal exterior is unchanged, and no final support lies outside the 77 core. If just one is internal, the edited induced count is an odd value at most 11. The possible one-arc predecessors of a 77-core in that range have count only 5: exclude 7 as above, 1 by the transitive formula, and check 3, 5, 9, 11 on the padded order-nine decks by monotone contraction. Their source counts and arc hits are \((7,0),(6,4),(20,0),(5,0)\). The induced order-four 5-core remains intact through the exterior edit, so it lies inside the final order-five 11-core.

**Compress the whole numerical walk.** All named support lies in the union of the 121 and 77 supports plus at most one vertex, giving at most \(10+9-5+1=15\) vertices. Delete vertices that are singleton SCCs at every named height. This preserves all named counts. If any deleted vertex carried an edit, the restricted numerical walk would have length at most eight, already impossible by the three disjoint detours. Thus all omitted vertices are untouched. This also preserves any endpoint isomorphism: such a vertex has the same outdegree at both endpoints, and singleton SCCs are distinguished from one another by their outdegrees. Pad the resulting walk to order 15.

There are 21 such height-121 peak templates, the placements of five pads around the two 11 cores. Choose the moving core and enumerate the \(99,55\) internal steps. The exact appearances of later heights are:

| Schedule | At 99 | At 55 | At 33 | At 77 | At final 11 |
|---|---:|---:|---:|---:|---:|
| \(A\) | 126 | 252 | 504 | 216 | 0 |
| \(B\) | 126 | 252 | 5460 | 0 | 0 |

The \(B\) entry at 33 counts all two-arc arrivals before the required strong-order-seven filter. At 77 the final two-arc test need only consider pairs with an internal arc. Multiplicities are retained. Ignoring the initial \(11\to121\) in this last scan is safe: the compressed candidate deck is a superset. The zero final hits rule out both nine-edit schedules.

**A rooted veto for schedule \(B\).** Its two edits \(55\to33\) would have to build a strong order-seven 33-core out of a height-11 order-five and height-5 order-four core, with possible pads. At least one edit is cross-component. If neither changes an arc inside these original cores, both strong supports survive, totaling nine vertices. If one edits the 11-core, its induced one-edit count is at least five, as checked on the unique 11-core, and it retains at least four strong vertices alongside the untouched four of the 5-core. So the internal edit has to act in the 5-core, and it has to make it transitive: even a surviving triangle together with the intact five would be too large. These two different arcs commute, so perform that edit first. The remaining jump is \(11\to33\) with one nontrivial initial SCC and gain 22. Two singleton endpoints spanning the 11-core would contribute \(Q_2\) divisible by four, impossible; endpoints both within the core cannot grow it. Thus the new order-seven core uses one endpoint in the 11-core, one singleton endpoint, and exactly one intervening singleton. The rooted gain then reads \(22=2U\) or \(2V\). On the unique 11-core, with rows \((18,4,9,19,6)\), the rooted vectors are
\[
U=(20,12,13,27,17),\qquad V=(17,27,13,12,20).
\]
Neither contains 11. This independently excludes \(B\); the full peak census above checks it as well.


There is also a sharper order-11 cap for \(A\). At \(33\), the original supports are an 11-core and triangle. If the internal arc in a two-edit \(33\to77\) lies in the triangle, commute it first: the other arc would give the forbidden \(11\to77\). If instead it lies in the 11-core, that core's one-edit reading is at least five. A reading above five already keeps its full five-vertex support. In the remaining case commute the internal edit first; the intermediate count is 15, with strong cores of heights 5 and 3. The cross edit must merge both, with gain \(77-15=62\) of valuation one. A singleton endpoint would leave one nontrivial intermediate core, spending at least two valuation units. Thus the endpoints must be in the two cores. If separated, the rooted triangle factor is 3, which cannot divide 62; if adjacent, the gain is at most \(3\cdot Q_2(H=5)/2=48\), even before subtraction. Here \(Q_2(H=5)=32\) for the unique strong four-vertex 5-core, and a fixed vertex cannot occupy the rooted position in both slots. Therefore all eight original strong vertices survive into 77. If the final 11 adds a vertex outside 77, there must be a single arc between the induced 5 and 77 cores. Its gain 72 has valuation three and nonexceptional ratio, so it can absorb at most four singleton pads; the 77 core then has at most eight vertices. In either case the support union has size at most 11. The order-11 replay reduces \(A\) to just three peak templates, with \(18,36,72,72,0\) arrivals at the successive checkpoints after the peak.

## The ten-edit order bounds

The only length allocations to test are the eleven obtained by adding one to a coordinate of \(A\) or \(B\) and deduplicating. With four order-eight padded 11 sources, no schedule even reaches 55 after the proposed 121 and 99 phases. With five order-nine sources and the labelled-return requirement, no schedule closes. The labelled search tracks the set of arcs flipped an odd number of times; if its size exceeds the number of remaining moves, exact closure is already impossible. Only the allocations \((2,1,2,1,2,2)\) and \((2,1,2,2,1,2)\) reach 55 under this pruning, six arrivals each; they reach 33 four and sixteen times, respectively, but never reach 77. A segment can use distinct unordered arcs: cancelling a repeated arc within it preserves named heights and would produce a forbidden walk of at most eight edits. Padding preserves each named height and the chosen closure convention, so the order-eight and order-nine obstructions establish both minimum orders.

The [prefix/peak deck replay](kaprekar-base12-prefix.js), [order-15 suffix replay](kaprekar-base12-peak.js), and [order and witness replay](kaprekar-base12-order.js) use only integer subset path counts and finite orientation enumerations.

## The same ten edits in infinitely many bases

The [exact lift rule](BASE16-SECOND-CYCLE.md#both-base-16-orbits-lift-together) asks for attained factors \(d=11x+1\) and \(e=13x+1\) to transport the base-12 walks to \(B=143x+12\). Their product is the exact ordinal multiplier \(d(B-1)/11\). Set \(x=2^k-2\), \(k\ge5\). The signed shell with a transitive filler of width \(w>0\) reads \(G2^w+D-L\); use the following two shells.

| Factor | Core outgoing rows | Tag mask | \((D,G,L)\) | Filler width |
|---|---|---:|---|---:|
| \(11\cdot2^k-21\) | \((14,28,24,16,33,15)\) | 28 | \((23,176,44)\) | \(k-4\) |
| \(13\cdot2^k-25\) | \((6,4,8,19,39,15)\) | 37 | \((15,104,40)\) | \(k-3\) |

Thus every base \(B=143\cdot2^k-274\), starting at 4302, carries the same ten-arc itinerary, scaled by \(de\). The labelled and relabelled constructions use \(2k+17\) and \(2k+16\) vertices respectively. This gives upper bounds in the lifted bases only.
