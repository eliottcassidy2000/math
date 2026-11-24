# Three ten-vertex parents suffice for every strong odd residue at eleven

**PROVED using the inherited exact lower bound and explicit finite certificates.** Eleven is the smallest order at which **strong** tournaments with transitive odd Hamilton-path precedence attain every odd residue modulo 24. One can keep both properties before and after adjoining a vertex. Exactly **three** ten-vertex parents are necessary and sufficient for this one-vertex construction, even if reducible parents with transitive precedence are also allowed.

The new vertex is the only one whose incidences vary on each of the three boards. Thus the residue chooses both a concrete strong parent and a single row of contacts. For a ten-vertex tournament \(T\), \(D(T)\) directs \(i\to j\) when an odd number of Hamilton words of \(T\) put \(i\) before \(j\). This theorem restricts the parents to those with transitive \(D(T)\); the twelve output heights below are ordinary Hamiltonian-path counts.

## A three-parent physical selector

The parent name gives its exact height. In each row the integers are outgoing bitmasks with vertices numbered from zero.

| Parent | Outgoing rows |
|---:|---|
|153|990, 1020, 1016, 1008, 992, 705, 896, 768, 32, 256|
|1563|542, 1020, 1016, 1008, 864, 705, 385, 785, 33, 320|
|1243|654, 1020, 1016, 1008, 737, 321, 897, 800, 529, 32|

Put vertex 10 above precisely those old vertices whose bit occurs in its contact mask. Each resulting tournament is strong and has transitive odd precedence. The direct replay enumerates its Hamilton words, counts odd precedence for every pair, and checks reachability.

| Requested residue | Parent height | Contact mask | Exact child height |
|---:|---:|---:|---:|
|1|153|256|289|
|3|1243|1008|2547|
|5|153|8|1349|
|7|1243|8|7783|
|9|153|258|2721|
|11|1563|1000|5363|
|13|153|510|2797|
|15|1243|1000|3807|
|17|153|760|929|
|19|1563|1008|3187|
|21|153|260|1821|
|23|1563|1010|6191|

Every one of these parents and children is rigid: any automorphism of a tournament preserves its precedence tournament, and a transitive tournament has trivial automorphism group. At every order through ten some residue is missing even without the strong condition. Indeed the [ordinary transitive-precedence census](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md#eleven-is-the-first-order-with-every-odd-residue-modulo-24-in-this-fiber) excludes residue 11 at ten; a smaller example would pad to ten in the same ordinary fiber. This proves the lower bound for strong witnesses as well.

## Why three parents are necessary

**FINITE-EXACT.** For each of the 911 ordinary ten-vertex classes with transitive precedence, vary every one of the \(2^{10}\) possible contact rows and keep just the strong children with transitive precedence. There are 11,203 qualifying rooted choices across this domain. Of the 911 parents, 468 are themselves strong. Here is the histogram of the number of distinct residues among the qualifying children.

| Number of residues | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|
| Number of ordinary parents | 112 | 156 | 167 | 175 | 164 | 137 |

Every six-residue support is the same set \(\{1,5,9,13,17,21\}\). Consequently two parent supports cannot union to size twelve. The three displayed parents have full qualifying supports \(\{1,5,9,13,17,21\}\), \(\{11,15,19,23\}\), and \(\{3,7,15\}\) respectively. Two parents in fact cover at most ten residues, achieved by the first two displayed supports. These counts refer to distinct marks on the parent extension; an unmarked child may have more than one qualifying parent.

The optimal triples are completely classified. Grouping the parents by their actual residue support gives 70 different supports. A cover by three must take one of each of the following rows, choosing either alternative in the last row:

| Role | Full residue support | Ordinary parent choices | Strong parent choices | Forced height |
|---|---|---:|---:|---:|
| Split | 1, 5, 9, 13, 17, 21 | 137 | 70 | varies |
| Inert | 11, 15, 19, 23 | 2 | 2 | 1563 |
| Inert, first alternative | 3, 7, 15 | 2 | 2 | 1243 |
| Inert, second alternative | 3, 7, 11 | 2 | 0 | 591 |

The proof enumerates all triples of distinct support sets among those 70. A repeated support cannot occur in a minimum cover since a pair never suffices. This gives exactly \(137\cdot2\cdot(2+2)=1096\) unordered covering triples of ordinary classes and \(70\cdot2\cdot2=280\) consisting entirely of strong parents. The two rare strong families are precisely the [two first-born converse pairs at residues 3 and 19](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md#eleven-is-the-first-order-with-every-odd-residue-modulo-24-in-this-fiber). Consequently **every** optimal strong-parent selector at eleven uses both of those doubled nodes at ten. The node masses are 3126 and 2486, and their ordinary heights are 1563 and 1243.

For reproduction, run `python3 -B strong11.py` for the explicit witnesses or `python3 -B strong11.py --census` for all extension masks. The [911-parent catalog](strong11_parents.txt) lists the upper triangular matrix in canonical precedence order, read by rows as a 45-bit binary integer and displayed in hexadecimal. The replay checks that the 911 encodings differ, that each displayed labeling has the ascending precedence order, and that the counts by residue agree with the complete ten-vertex census. This certifies that the catalog exhausts the required parent classes; the unique precedence order rules out isomorphic duplicates.

For each parent the replay evaluates 1024 lanes together using integer bit operations. Directed cycles are evaluated on each vertex subset, followed by the superset parity sum on each distinguished pair. The [pair-order theorem](../chanpak_tiling_viewpoint_20261123/ORDER-PARITY.md) gives the child precedence orientation by adding this parity to its arc indicator. The absence of cyclic triples recognizes transitivity. To recover the integer heights, split a word before and after the apex. A prefix ending at \(i\) and a suffix starting at \(j\) are compatible precisely when \(i\) beats the apex and the apex beats \(j\). Summing the independent prefix and suffix counts over complementary subsets gives the [quadratic contact table](../chanpak_arc_fiber_grammar_20261122/BIRTH-CUT-COVER.md#one-parent-supplies-a-quadratic-table-of-counts). Thus both the group of six split residues at the most prolific parents and the need for additional parents to supply the inert residues reflect restrictions of actual contact placements.

## The obstruction to repeatedly using the same extension rule

**FINITE-EXACT.** Suppose the final tournament has eleven vertices. Insist that it be obtained by adding one vertex at a time and that every intermediate tournament, including the initial one, be strong with transitive odd precedence. Allow every initial class at the given order. Beginning at nine or ten, all twelve residues occur. Beginning at any order from four through eight, the entire attainable list at eleven is
\[
                         \{1,5,9,13,17,19,21\}\pmod{24}.
\]
Thus the earliest starting size for complete residue selection by such a hereditary chain is nine, even though smaller transitive-precedence components can be joined in other constructions. At orders two and three there is no strong tournament with transitive precedence; a directed triangle has cyclic precedence.

To see the finite reduction, keep the set of attainable final residues attached to each of the 468 strong ten-vertex parents. For each vertex deletion, retain the smaller tournament precisely when it is strong with transitive precedence, canonicalize it in that order, and transfer the whole set of final residues downward by union. An induction on the deletion depth proves that this records exactly the admissible chains. At starting orders 9, 8, 7, 6, 5 and 4 the numbers of distinct ancestor classes with a surviving chain are respectively 126, 29, 12, 3, 1 and 1. Thus the lost residues at eight occur even when every qualifying eight-vertex starting class can be used. The same replay computes these sets using independent direct odd-word recurrence on the smaller parents.


At the earliest hereditary starting size nine, exactly three initial classes are again necessary. The supports obtained after two qualifying extensions have only two types of covering triples. One initial class must attain each of the sets \(\{11,15,19,23\}\) and \(\{3,7,15\}\). There are two choices for each, forced to have respective heights **751** and **591**. The third class can be any of 78 whose support is the six split residues, or any of 8 whose support additionally includes 19. Thus precisely \((78+8)\cdot2\cdot2=344\) unordered triples work. No two supports cover all residues. The replay verifies these claims by grouping the 126 nonempty nine-vertex supports and checking their unions.

This gives two independent reasons for the same nine-vertex height 751. Every optimal three-source bank in this hereditary residue problem uses it, and the [shifted-square intersection](BANK.md) requires its prime factor in the only six-block solution \(113^2-2=17\cdot751\). The first conclusion depends on a path through two actual apex extensions; the second depends on multiplicative prime support and the norm-two Pell equation.
