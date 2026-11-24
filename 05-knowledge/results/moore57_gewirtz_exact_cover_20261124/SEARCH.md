# The unrestricted 56-cover reduces to two seed classes

**PROVED reduction with finite class inventory.** This is not an infeasibility certificate. Both seeds together cover every possible sharply transitive set up to the specified moves.

For a derangement g of the Gewirtz sheet graph define e(g) to be the number of sheet labels x with g(x) adjacent to x. Class enumeration in [seeded-cover.py](seeded-cover.py) gives:

| Geometric type | e=8 | e=14 | e=20 |
|---|---:|---:|---:|
| Collineation | 7560 | 11520 | 0 |
| Correlation | 16800 | 0 | 5040 |

For a sharp set S and each g in S, the normalized set g⁻¹S is again sharp. Since the sheet graph is 10-regular, its 55 nonidentity maps together move labels onto 56·10=560 ordered adjacent cells. Define n14(g), n20(g) as the counts of the indicated e values in g⁻¹S. Subtracting 55·8 yields

    n14(g) + 2n20(g) = 20.

This is a pairwise identity on S. In particular, each g has at most twenty same-geometric-type partners whose relative movement is 14. One of the collineation and correlation parts of S has at least 28 members, regardless of the congruence on the split. Any g in that part must have some distinct same-type partner h with e(g⁻¹h)=8. Translate by g⁻¹, then conjugate to a representative of the appropriate full G conjugacy class. The replay finds exactly two such classes, of sizes 2520 and 5040. Thus it is enough to force identity and one of two fixed seeds, and delete every candidate that overlaps either. The two candidate universes are 21,172 and 20,292, each covering the same number (3024) of remaining cells. To recover the full set, seek 54 maps, not 55, one cover of every cell.

The default command just generates and verifies the reduced counts. Optional bounded search:

    python3 seeded-cover.py --class-index 0 --seconds 100 --witness seed0.npy
    python3 seeded-cover.py --class-index 1 --seconds 100 --witness seed1.npy

The initial 40,920-candidate unrestricted HiGHS search reached its 240-second limit without an incumbent or infeasibility certificate. Both 100-second seeded runs also reached their time limits without an incumbent or infeasibility certificate. None of the three HiGHS runs got past the initial LP node. Exact structure, not the duration of these runs, motivates the two-seed split; it also suggests moving the next unrestricted attempt to a direct exact-cover or SAT solver. Successful solver output is independently checked for sharpness before saving. To claim nonexistence, both seeds need an exhaustive check or independently verifiable certificate.

## A positive sign-only test, and its sharp boundary

The prior [ternary sign equation](../moore57_switching_flags_20261124/GEWIRTZ.md) is not contradictory by itself. Put the actual Gewirtz adjacency K on the 56 nonanchor branches. Its parameters (56,10,0,2) give

    K² = 8I − 2K + 2J ≡ K − J − I (mod 3),

and its row degrees are 10. In the anchored chart the correlation counts in the relative matching families are then exactly 10 for an anchor pair, 16 for a nonedge, and 37 for an edge: the two row supports off an edge have 18 differing bits (reverse the sign at that edge), and off a nonedge they have 16. All pass the numerical bounds and congruences. This is a sign template only, not a permutation witness.

In fact this is the **smallest maximum anchor degree** for any sign template obeying the known bounds and ternary law. If every degree is ≤10, admissible degrees are 7 or 10. The ternary law forces every nonadjacent pair to have a common-neighbor count congruent to 2, hence at least 2. A degree-7 vertex would need at least 2·48=96 two-step walks to nonneighbors but can have at most 7·9=63. Thus all degrees are 10. Now the 45 nonneighbors already use all 90 nonbacktracking two-step walks, forcing no triangles and exactly two common neighbors at nonedges. So the extremal template has the Gewirtz parameters. A search restricted to this template must support sharp sets of the three displayed correlation sizes; the unrestricted search need not fall in this minimum-degree sign chart.
