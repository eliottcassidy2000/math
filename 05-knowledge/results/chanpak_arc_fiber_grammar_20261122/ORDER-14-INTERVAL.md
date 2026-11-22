# The complete odd prefix through 3,423,417 on fourteen vertices

**FINITE-EXACT.** Every positive odd integer through 3,423,417, except 7 and 21, occurs as the number of Hamiltonian paths of a tournament on exactly 14 vertices. For every such height greater than one, every tournament of minimum order realizing it is strong. In particular every odd height from 1,039,171 through 3,423,417 has minimum order either 13 or 14. The [order-13 interval certificate](ORDER-13-BIRTH-INTERVAL.md) determines the minimum exactly for more than a million of these heights. The endpoint 3,423,417 and the earlier two-apex boundary 2,327,675 both have minimum order exactly 13.

The prior unrestricted prefix reaches 4,294,967,295. The present result supplies a uniform order bound on part of that range. The earlier two-apex word census gives every admissible odd through 2,327,673 on at most 16 vertices in its restricted class. At 14 total vertices its least additional omission is 230,521. A separate single-target construction realizes 11,390,625 on 14 vertices. [Corrected unrestricted prefix](../../../session203/coverage-recertification/README.md) · [Two-apex order census](../session216_overnight_constructions_20261107/two-apex-advances/ADVANCES.md) · [Single high value on fourteen vertices](../../../session208/balanced-order14/README.md)

## Exact cover

There are 1,711,707 requested heights. The order-11 census and singleton tournament already supply 43,911 of them. The other 1,667,796 heights run from 80,453 through 3,423,417. Use precisely the two deterministic parent schedules in the [component-pair cover](BIRTH-CUT-COVER.md). The first schedule examines all apex contact masks on 20,000 parents of orders 12 and 13. It leaves the following 16 heights.

    2856797 2934439 2949955 3112239 3163127 3206281 3241091 3279609
    3289311 3292653 3296165 3299811 3308559 3378997 3398455 3420785

The second schedule examines 1,057 further parents of order 13. It supplies all sixteen. The parent row masks and apex contact mask specify an actual tournament at each equality. Appending an ordinal singleton preserves the count, giving exactly 14 vertices for every requested height.

For the strong conclusion, a smallest representative of height greater than one contains no singleton strong component: deleting such a component preserves the count. By the [two-component size obstruction](BIRTH-CUT-COVER.md), two nontrivial components require combined order at least 16 in any smallest representative. Thus every smallest representative in the interval is strong. The largest order-11 count is 95,095. Choosing the first vertex of a Hamiltonian path bounds every order-12 count by 12 times this, or 1,141,140. The sharper first-two-vertices bound is 66 times the order-10 maximum, or 1,039,170. Padding makes both upper bounds valid at every smaller order as well.

## Two exact births at thirteen vertices

The following outgoing masks, indexed from zero, give strong order-13 tournaments. The direct endpoint recurrence verifies their counts and the upper bound at order 12 proves their minimality.

| Height | Outgoing masks |
|---:|---|
| 2,327,675 | 2372, 6221, 5744, 1877, 7203, 1867, 912, 3391, 1558, 5267, 6211, 4972, 489 |
| 3,423,417 | 4554, 5932, 3161, 5808, 6755, 6341, 3466, 3606, 700, 1125, 6449, 779, 3012 |

Both heights have **at least two** smallest isomorphism classes. Reversing all arcs preserves the path count by reading each path backwards. For the 2,327,675 witness the largest outdegree is 9 and the smallest is 4, so the converse has largest outdegree 8. At 3,423,417 the full score sequence is symmetric. The nine vertices of outdegree six induce a tournament with score sequence \((2,4,4,4,4,4,4,5,5)\). The corresponding sequence in the converse is \((3,3,4,4,4,4,4,4,6)\). An isomorphism must preserve this vertex set, proving the second pair is inequivalent as well.

The first height requires 17 vertices inside the two-apex word class. Its unrestricted minimum is four vertices smaller. The endpoint witness was found by starting from the regular cyclic tournament on 13 vertices, reversing cyclic triangles to retain the score sequence, deleting one vertex and ranging over all ways of restoring its contacts. The displayed rows suffice to reproduce its exact count. [Restricted minimum at 2,327,675](../session216_overnight_constructions_20261107/two-apex-advances/ADVANCES.md)

Run `python3 -B birth_cut_cover.py interval14` from this directory. The program reads the inherited census, verifies both order-13 controls, and compiles the C++17 all-mask oracle in a temporary directory. The second stage must leave zero residual targets. For any height in the residual set, append `--witness H` to return and independently recount one actual graph from the cover. [Driver](birth_cut_cover.py) · [Quadratic endpoint oracle](birth_cut_cover.cpp)
