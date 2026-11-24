# Exact failure of translation lifts and odd two-sided symmetries

Let G be the degree-56 Gewirtz automorphism group, of order 80,640, on an orbit of hyperovals of PG(2,4). This is the group in the [standard Gewirtz description](https://aeb.win.tue.nl/graphs/Sims-Gewirtz.html); the replay builds the projective plane, elementary transvections, Frobenius and exterior-line duality and enumerates its action. The independent [PrimGrp](https://gap-packages.github.io/primgrp/doc/chap1.html) generator lists give a calibration for C2 and C7. The results below are finite-exact for G, and therefore apply also to any subgroup under the appropriate symmetry hypothesis.

A set S of 56 permutations is **sharp** when for every ordered cell (x,y) exactly one g in S satisfies g(x)=y. If it exists, it can be translated to contain the identity, so its other 55 permutations are derangements. Any two distinct members disagree everywhere. Left and right multiplication preserve sharpness.

## One-sided symmetries: no common coset construction

Suppose S is invariant by left translation by a prime-order h. Then the order divides 56 and h must act freely on labels, since g and hg agree at a label mapping to any fixed point of h. Conversely, if h has orbits of equal size p and c is their coloring, the filling by the entire left coset \(<h>g\) uses every color class selected by c(g(x)) exactly once from x. A sharp coset union is therefore 56/p colorings with pairwise distinct colors at *every* coordinate. Right translation is equivalent by inversion, and any nontrivial multiplier subgroup has order dividing 56 because its action on S is free. Thus p=2,7 suffice.

For p=7, normalize one coloring to c. There are 201 distinct colorings compatible with c; the compatibility graph among them has 252 edges and no triangle. Its largest clique has size two, giving max three cosets with c, rather than eight. The [replay](gewirtz-sharp-c7.py) also checks the four proper catalogue actions.

For p=2, there are 280 fixed-point-free involutions, one conjugacy class (centralizer size 288). Fix any one h and normalize one coloring to c. Among the 40,320 distinct colorings, 10,152 are compatible with c. The centralizer permutes the compatible colorings by permuting inputs and relabeling colors; the data replay verifies 41 orbits. Any hypothetical clique of 27 additional colorings can move one of its vertices to an orbit representative. A bitset clique search starts from each representative and branches with a proper-coloring upper bound: color classes are independent, so a clique uses at most one vertex from each class. It visits 90,184 nodes and fails; the [data generator](gewirtz-sharp-c2-data.py) and [clique replay](gewirtz-sharp-c2.cpp) contain the complete test. This is a finite proof, not a timed heuristic search.

Thus a hypothetical S has trivial left and right multipliers individually. In particular there cannot be a sharp *subgroup* in G, but failure of all such coset constructions is stronger.

## Two-sided symmetries: no order 3, 5, or 7

An order-p element (a,b) acts by g ↦ a g b⁻¹. Each of p=3,5,7 has a single class of nonidentity order-p elements in G. The [replay](odd-stabilizers.py) checks this by counting order-p elements and the centralizer of one element (the counts multiply to |G|). When both a,b are nonidentity, left and right multiplication separately reduce to conjugation by one h, up to translating S. If only one is nonidentity, the one-sided results and the order dividing 56 rule apply.

For h acting diagonally on cells, each nonfixed orbit of permutations must cover whole orbits of cells. Its support is legal iff its representative hits 56 distinct nonfixed cell orbits: otherwise the expanded supports overlap. A fixed cover member lies in the centralizer; if one exists, translating by it preserves invariance and normalizes that member to identity. If h fixes f>0 labels, an invariant sharp set has exactly f fixed permutations: at any fixed input, the unique maps to fixed outputs are fixed cover members, and a commuting map must send that input to a fixed output. If f=0, the number of fixed cover members is still congruent to 56 mod p.

| p | Number of quotient cells | Legal nonfixed orbit footprints in normalized case | Fixed choices after normalization |
|---:|---:|---:|---:|
| 3 | 1048 | 5250 | 15 choices for the second fixed member |
| 5 | 628 | 1106 | identity only |
| 7 | 448 | 507 | 64 compatible choices of six further fixed members |

Every case fails by branching on an uncovered cell with the fewest admissible footprints. In the extra p=7 case with no fixed permutation, allow every group element initially (no identity normalization): there are 1102 legal nonfixed footprints and that cover search fails too. The algorithm has no time cutoff and excludes all covers in each listed case. In particular the entire setwise stabilizer in G×G, if nontrivial, must be a **2-group** projecting injectively to each factor. This does not exclude all diagonal involution models or all arbitrary sharp sets.

## Consequence for a hypothetical global Moore symmetry

Suppose the holonomy at a root r is a degree-56 rank-three Gewirtz action with subdegrees 1,10,45, as in the surviving catalogue cases. An automorphism of odd prime order p fixing r and at least two of its neighbors cannot exist in this chart. In the matching gauge anchored at one fixed neighbor, preservation of the identity matchings forces the same label permutation h on each fiber. It normalizes the generated holonomy. Its unequal rank-three suborbits imply the normalizer preserves the Gewirtz adjacency relation (the unique nontrivial orbital of valency 10), so h∈G. For a second fixed neighbor the sharp set indexed by the intermediate branches is conjugation-invariant under h. If h=1 the branch permutation also fixes every branch, as the maps in a sharp set are distinct, and the automorphism is trivial. Otherwise p divides |G| and hence p=3,5,7, precisely the cases just excluded.

For example, this rules out the Gewirtz transport **at a root inside** a fixed Petersen graph (order 3), fixed pentagon or Hoffman–Singleton graph (order 5), or the fixed center of a star with at least eight leaves (order 7). The possible fixed-subgraph types are described in [Ishida, §2.2](https://arxiv.org/html/2606.29183v2#S2.SS2). It does not exclude the same abstract automorphisms in a chart with alternating or symmetric transport, or cases where there is only one fixed branch.
