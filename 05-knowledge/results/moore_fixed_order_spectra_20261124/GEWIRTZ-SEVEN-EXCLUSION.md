# The seven- and 49-duality Gewirtz cases are excluded

This is a finite exclusion in the **full 80640-element Gewirtz action** for a sharply transitive set of 56 label transports, not an exclusion of the whole action at all swap counts, nor a claim about existence of a degree-57 Moore graph. It completes the previous [seven-duality pencil/Fano branch](FANO-AND-FOUR-SEAM-DECODER.md).

## 1. Necessary geometry and exhaustive seven-set reduction

If a sharply transitive set has seven type swaps, their projective point images at *every* source point, and their line images at *every* source line, satisfy the pencil-or-Fano rule proved in the preceding pass. Translate on the left by the inverse of one of the seven swaps. The seven become point-preserving maps including identity, the other 49 become type-swapping maps, and the local pencil/Fano conditions are transported by a correlation. In particular the seven normalized maps must be pairwise deranging on the hyperovals. Each of the six other nonidentity maps is among the 19080 derangements of the point-preserving group.

Enumerate the six point-preserving conjugacy classes, of sizes 2520, 5040, and four of size 2880. Conjugation preserves both projective source sets and the geometric rules. If any pair of maps agrees on a projective point, translate within the point-preserving group so that pair becomes identity and one fixed-projective-point representative in one of the first two classes. This translation preserves the target geometry. After testing those two cases, one may restrict the remaining cases to no projective point coincidences. A line coincidence would give a point coincidence after conjugating by a duality (the point-preserving group is normal), so the remaining cases need **no** coincidence on either type. In those remaining cases all 42 sources must therefore be Fano.

The exact replay takes one of these representatives as the first nonidentity map, and adds other maps in increasing index order. It intersects the pairwise hyperoval-derangement conditions and stores the current image multiset at all 42 projective sources. A map is kept only if each source can still extend to a permissible final multiset. The image masks are computed over the explicit \(\mathrm{PG}(2,4)\) from the [parent replay](GEWIRTZ-AND-DUAL-GUARDS.md).
 
- A current distinct support can be extended to a Fano set containing it, or to a line of five points with one point doubled outside or tripled inside. Exact masks for the next image are the union of these geometric choices.
- Once a duplicate appears, all other images must lie on one line, and a third copy of the duplicate is allowed exactly when it can also be on such a line. Once it is tripled, use only lines containing the full current support and add no repetition.
- In the coincidence-free case, use only subsets of the 360 Fano subplanes. At any completed seven-set, independently check odd incidence with all 21 projective lines from all 42 sources. This rechecks the final geometry, apart from the pruning.

Membership in the candidate lists is checked by masks, not by sampling or a numerical criterion. The 19080-vertex relative derangement graph has 85,142,880 edges. A bitset filter on eight of the source coordinates is just an acceleration: all surviving vertices are tested on all 42.

The two fixed-one representatives produce **no** seven-set: the combined search reached 35,140,682 states with three selected nonidentity maps (plus identity), and no extendible state at the next level. The four coincidence-free representatives produce respectively **one, one, zero, zero** normalized leaves. Both leaves, together with identity, are cyclic subgroups of order 7, closed under composition. Both have all 21 points and all 21 lines as Fano sources. The earlier geometric condition alone therefore **does not** exclude seven dualities; the last exact-cover step is essential.

## 2. Exclude both exact-cover completions

For each surviving normalized seven-set, complete the \(56\times56\) all-ones matrix by **49 type-swapping permutations**. Enumerate the entire type-swapping coset, retaining only permutations that do not occupy any of the seven old entries of a row. Both leaves leave 2744 cells and 819 candidate permutations; at the initial state every remaining cell has at least eight candidates.

Use exact cover: pick any remaining cell with the fewest surviving placements, branch on those placements, remove their 56 cells and every intersecting placement, and continue. Empty available list fails, empty uncovered cell set succeeds. This covers all potential completions because every missing cell has to be covered exactly once. The full branching search fails after 36 and 37 visited states respectively. Hence neither seven-set extends. Translating by a duality also excludes 49 swaps. As a result the remaining necessary count of swaps per branch pair is \(10,13,\ldots,46\). A valid deck, or even a sharply transitive set in this action at one of those other counts, remains unknown here.

The stronger **unrestricted completion** test also fails: after allowing either type for each of the remaining permutations there are 1407 candidates, each missing cell having at least 24 initially. The exact search fails in 461 and 502 nodes on the two seeds. Since all order-seven subgroups are Sylow conjugate, **no sharply transitive set in the full group contains even one complete left coset of a subgroup of order seven** (translate a coset to its subgroup). This extends the older obstruction to unions of cosets. It does not rule out sets distributed among partial cosets.

Reproduce in this directory:

    python3 gewirtz-seven-export.py gewirtz-seven.bin
    c++ -O3 -std=c++17 gewirtz-seven.cpp -o gewirtz-seven
    ./gewirtz-seven gewirtz-seven.bin > seven.out
    python3 gewirtz-seven-cover.py seven.out

The optional class argument of the C++ replay reruns only one of the six cases; the full conclusion uses all six. The two tail-variable tournament results from this session are recorded [separately](VARIABLE-TAIL-SELECTOR.md).
