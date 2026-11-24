# When two arc edits fit into one Hamilton-path viewpoint

Fix a labelled tournament and a set of unordered pairs whose directions may vary. Write \(D\) for the directed graph obtained by erasing those pairs. A **common viewpoint** is a directed Hamilton path of \(D\). In its positional coordinates every completion is a fixed-path tiling, and the variable pairs are free tile coordinates. Existence of separate paths in the individual completions is automatic. The requirement here is one path retained in every completion.

In the complement-line convention, \(\kappa\) flips all nonpath pairs. Staircase reflection \(\sigma\) sends the tile at \((i,j)\) to \((n-1-j,n-1-i)\). The line joining \(t\) to \(\kappa t\) is blue precisely when \(\sigma t=t\), and is black otherwise. The reflection exchanges converse classes. A single class has \(H/|\operatorname{Aut}|\) tilings. Thus every common-viewpoint cube has a second cube consisting of its line partners, with the same editable positions. A single tile edit and the complement-line involution are different operations. [THM-643, blue and black line parity](../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md) · [THM-4731, the staircase reflection and black twins](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)

## One erased pair

**PROVED.** An arc lies on every Hamilton path of a tournament precisely when its endpoints are consecutive singleton strongly connected components (SCCs) in the ordinal component order.

Inside a nonsingleton SCC take a directed Hamilton cycle. If the arc belongs to that cycle, cut the cycle at the arc. If it does not, cut anywhere. This gives a spanning path of the component avoiding the arc. A Hamilton path of the full tournament concatenates spanning paths of its SCCs. It uses no arc between nonconsecutive components. At consecutive components an arc is avoidable whenever one component has at least two vertices, since rotating its cycle gives any ending or starting vertex. Adjacent singleton components necessarily occur consecutively on every path. Their transposition identifies the tournaments before and after the reversal.

Consequently the Hamming-distance-one quotient of the fixed-path tiling cube has exactly the nonloop adjacencies of the one-arc isomorphism metagraph at every order. Every edge between distinct classes has a common viewpoint. This proves the all-order implication left conditional in section 5 of the dictionary. [THM-1430, tiling/class/metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md)

## Two erased pairs

**PROVED.** Suppose \(D\) has at most two missing unordered pairs and all other pairs carry one arc. The following conditions are equivalent.

1. All its completions share a Hamilton path.
2. The SCCs of \(D\) are totally ordered by reachability.
3. For each missing pair, in \(D\) one endpoint can reach the other.
4. No completion has a missing pair's endpoints as consecutive singleton SCCs of that completion.

In the fourth condition the designated pair has of course been given its direction. Checking the third condition uses at most four ordinary graph searches.

A spanning directed path compares all its vertices by reachability. All nonmissing pairs are already comparable. This gives the necessity of (2) and its equivalence with (3). Suppose (2). In the unique component order there is a direct connection between consecutive components. Let \(b\) be the number of missing pairs between different components.

* For \(b=0\), every component has a spanning path. For an incomplete strong component use the Chen–Manalastas theorem, since its stability number is at most two. Complete cross-arcs concatenate the paths.
* For \(b=1\), at most one pair is missing internally. Every nonsingleton tournament component has a spanning path with any chosen start or end. An incomplete strong component with one missing pair has at least two available starts and at least two available ends, by Havet's endpoint theorem. If the absent cross-arc is at a consecutive boundary, at least one of its endpoint components is nonsingleton, since otherwise total reachability fails. Choose a different endpoint there. All other boundaries are complete.
* For \(b=2\), all components are tournaments. Independently choose a uniformly random cut of one fixed Hamilton cycle in each nonsingleton component. At a specified boundary the chance that the resulting end and start use a particular missing pair is at most \(1/3\), since some endpoint component has at least three vertices. A missing pair away from a boundary has probability zero. The sum of the two failure probabilities is less than one, giving a successful choice.

This proves (2) implies (1). The cited stability theorem applies to any finite strong digraph of stability two. Havet's stronger result supplies two nonadjacent initial vertices and two nonadjacent terminal vertices. [Chen and Manalastas, 1983](https://www.sciencedirect.com/science/article/pii/0012365X83901887) · [Havet, 2004](https://doi.org/10.1016/j.disc.2004.07.013)

If (2) fails, incomparability of two SCCs requires every pair between them to be missing. An SCC here has order one or at least three, so the two incomparable components must be singletons. There is a linear extension of the component poset with these two adjacent: list all strict predecessors of either one, then the two, then extend the remainder. Fill cross-component pairs forward, and fill internal pairs arbitrarily. The two singled-out vertices are consecutive singleton SCCs of this completion. Conversely such a completion cannot have a spanning path avoiding its forced arc by the one-pair theorem. This proves (4).

**Geometric consequences.** A two-edit square with no common viewpoint has a whole side whose two completions differ only by a transposition of consecutive singleton components. In particular any two-edit square with different Hamilton counts across all four sides has a common viewpoint. So does a square with strongly connected tournaments at two diagonally opposite corners: the forced side would have two reducible corners and must meet that diagonal. The conclusions apply to the given labelled edits, retaining the actual incidence at the middle class.

**Chan–Pak unit steps.** Suppose two distinct edits from the same tournament each change its Hamilton count by exactly two in absolute value, and all three counts are at least five. Then the edits occupy two tiles in one common viewpoint. Indeed the inherited strong-core theorem for each edge says that both its orientations consist of a strong active core and singleton exteriors, with the same active vertex set for that edge. The middle tournament has only one nonsingleton SCC. Both edits are therefore in its same core, and reversing either one keeps that core strong. These two neighbors occupy opposite corners of the core square. The preceding result supplies a common path there. Append the singleton exteriors. This reduces the search for composable unit steps to two coordinate flips at one fixed path, together with their antipodal line partners. Constructing the middle tournament is still necessary. [Strong active core for unit flips](../chanpak_arc_fiber_grammar_20261122/README.md)

The small-height exception is real. Take outgoing masks \((0,5,9,3)\), mark the pairs \(\{1,2\}\) and \(\{0,1\}\), and index vertices from zero. The four counts in unedited, first, second, both order are \((3,1,5,1)\). The two steps form a coherent \(1\to3\to5\) itinerary on these labels, but no viewpoint survives both. At order three all counts are at most three, so this is a smallest such monotone unit example.

## The boundary at three holes

**PROVED.** With an arbitrary number of missing pairs that form a triangle-free undirected graph, strong connectivity of the fixed directed part guarantees a common viewpoint, by Chen–Manalastas. The hypothesis about triangles is sharp already for three holes. Take three independent vertices \(X\), a vertex \(a\) and a vertex \(b\), with \(X\to a\to b\to X\). This digraph is strong. Going from one member of \(X\) to another consumes both \(a,b\), so a path visits at most two of the three. All eight fillings are strong tournaments but have no common viewpoint.

Total reachability alone also fails for three triangle-free holes. In the ordinal sum \(s\to C_3\to t\), name the directed cycle \(a\to b\to c\to a\). Erase \(sa,at,bt\). All components remain totally ordered, yet a common path would have to start at \(s\), enter the cycle at \(b\) or \(c\), and exit it at \(c\). The spanning paths entering at \(b,c\) end at \(a,b\), respectively. This exhibits the unremembered condition: a component may permit the required entry and required exit individually while forbidding them jointly.

## Blue and black placement on the same two-edit square

For a fixed common path, let \(F\) be the free tile coordinates. A cube has zero blue corners if some reflection orbit of coordinates has contradictory pinned bits. Otherwise it has exactly \(2^q\) blue corners, where \(q\) counts the reflection orbits entirely contained in \(F\). Indeed each blue filling must be constant on each reflection orbit. A partly pinned compatible orbit forces the free bits on it. Thus, for example, a free mirror pair gives two opposite blue corners when the pinned part is symmetric. All corner colors are preserved under \(\kappa\).

Use the five-vertex height-13 class from the earlier quadratic-character construction with masks \((28,17,10,18,4)\). Mark \(\{0,3\},\{1,4\}\). Its four counts are \((13,11,15,13)\). There are precisely seven common path words for this one square. Changing the path leaves those four counts fixed and can change all the complementary targets. In the table, bit positions in the last column mark the blue corners in the same unedited, first, second, both order. [Inherited character fork](../square36_signed_collatz_20261123/CHARACTER.md)

| Common path | Four complement-partner counts | Blue indicators |
|---|---|---|
| 04231 | 13, 9, 15, 13 | 0110 |
| 10234 | 13, 11, 5, 5 | 0000 |
| 10423 | 9, 5, 5, 3 | 0001 |
| 23104 | 13, 11, 9, 9 | 0000 |
| 31042 | 9, 11, 9, 13 | 0000 |
| 34210 | 5, 11, 5, 13 | 0000 |
| 42310 | 3, 5, 5, 9 | 1000 |

The common-path test exactly classifies which two-edit squares can be drawn. Numerical transport along the blue or black lines further depends on the chosen common path. Run `python3 -B lines.py` for the exact examples and the small exhaustive comparator.
