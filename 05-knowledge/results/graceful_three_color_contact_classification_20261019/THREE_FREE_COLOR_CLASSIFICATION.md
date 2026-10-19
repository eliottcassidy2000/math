# Three free colors: finite contact classification

Let F be a fixed forest with exactly four components, and let three disjoint edge-color classes admit a spanning-tree completion. Contract F, discard loops, and retain for each color its set S_i of unordered pairs of distinct components. Multiple actual edges of the same color and pair are distinct realizations of one contact type.

If any S_i has graphic rank at least three, the earlier three-bad-color rank criterion gives width-two connectivity. In the remaining case all three supports have rank at most two.

**Finite classification.** For three nonempty rank-at-most-two supports on four labeled vertices, the nonempty rainbow spanning-tree type graph is disconnected if and only if there is a spanning tree with edges x,y,z such that, after permuting the colors, the supports are exactly

    {x,y}, {y,z}, {z,x}.

There are exactly two valid contact-type states, the two cyclic assignments. They differ in all three colors. Every other nonempty type graph is connected under changes in at most two coordinates.

## Exact certificate

A simple graph on four fixed vertices has six possible edges. Its nonempty supports of rank at most two consist of six single edges, fifteen two-edge sets, and four triangles: exactly 25 supports. The accompanying standard-library program enumerates all 25^3=15,625 ordered support triples. It lists independent selections by testing connectivity of three selected edges on four vertices, then computes the exact type-state components. For three-coordinate states, width-two adjacency is equivalent to agreement in at least one coordinate.

For every disconnected case, it checks that the union is a three-edge spanning tree, each support has size two, and the supports are pairwise distinct. Those conditions are exactly the displayed cyclic pattern. There are 96 disconnected triples, agreeing with 16 labeled trees times six assignments of their three pair supports to the colors. The explicit cyclic pattern proves the converse without computation.

This is a finite, computer-assisted classification, not a claim that a noncomputational classification proof has been supplied.

## Returning to actual edge realizations

Changing an edge to another realization with the same color and component pair preserves the quotient tree and is a one-color move. Every type-state fiber is therefore connected. Conversely any actual width-two move projects to a type move (or a repeated type state). Any type move can be lifted by choosing arbitrary available realizations in its changed coordinates and preserving all unchanged realizations. Therefore actual and type reconfiguration graphs have the same connected-component structure.

## Integer-difference question reduced exactly

For an arbitrary fixed-color forest on labels 0,...,m leaving three colors free, nonemptiness implies it has four components. A disconnected completion fiber must therefore realize exactly the cyclic three-pair support above; no extra contact type may be available in any of the three colors. Each color can have any number of loop edges inside fixed components and parallel realizations of its two allowed pair types.

The classification itself does not decide arithmetic realizability. The finite interval checks through m=8 found none; subsequent targeted work produced an m=12 realization, staged separately. Empty fibers are excluded throughout.
