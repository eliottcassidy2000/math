# Selective common second returns, including an exact strong times-three move

The strict cover order of the [two S320 returns](FOUR_VISIT.md) rules out a common substitution **of the entire returned graph** into a strong quotient. It does not rule out opening the prior quotient and changing one side module. The following are literal examples.

Label each returned graph by putting the original source vertices first, then the nonmarked quotient vertices in their original order. Add one vertex with a common out-cut U on these labels. Among all nonempty proper common cuts, the equalities are exactly:

| old equal count | common cut mask | new equal count | refined quotient visit counts |
|---:|---:|---:|---|
| 532623 | 7167 (all except vertex 10) | 1068387 | (29,24,4) |
| 4591771 | 2048 (vertex 11) | 7655505 | (43,40,10) |
| 4591771 | 5120 (vertices 10,12) | 13775313 | (75,72,18) |
| 4591771 | 6144 (vertices 11,12) | 9189287 | (53,48,12) |

The first three visit counts refer to the old source module, not to the entire returned graph. Each listed new graph is strong because the old graph is strong and the added vertex has at least one in-neighbor and at least one out-neighbor. The equality can be certified directly by the refined quotient visit counts: for the first source pair the differences in (F1,F2,F3) are (0,240,−1440); for the prime pair they are (0,720,−2880). No refined quotient has a four-run itinerary. The standalone layout replay enumerates the quotient words.

For the exhaustive assertion, write h_i^+(S), h_i^−(S) for the counts of paths starting or ending at vertex i on S, set to zero when i is absent, and define
$$
C_{ij}=\sum_{S\subsetneq V,\,S\ne\varnothing}h_i^-(S)h_j^+(V-S).
$$
The extension count is
$$
\sum_{j\in U}h_j^+(V)+\sum_{i\notin U}h_i^-(V)+\sum_{i\notin U,j\in U} C_{ij}.
$$
This partitions by whether the new apex starts the path, ends it, or separates two nonempty slots. The cut replay evaluates this exact quadratic form for every shared nontrivial cut: 8,190 at order 13, 32,766 at order 15. It also finds the smallest strong extension increments 471288 and 472008 for the two order-13 graphs, and 3063734 for each order-15 graph. In each case the minimum is at the unique mask 2048. Hence none admits a strong count-preserving one-vertex extension individually. This does not restrict other label alignments, different cuts on the two graphs, or simultaneous edits of existing arcs.

## Why the exact times-three move works for any old source

In the prime-return quotient (44,33,2,6,15,28), vertices 2 and 3 form a transitive pair that is a homogeneous module, with sole exterior out-neighbor 1. The cut 5120 on the full graph corresponds to quotient cut {1,3}: the new vertex beats 1 and 3 and loses to every other old vertex. It joins the transitive pair to form a homogeneous directed triangle 3→2→new→3. For every choice of the old marked source, a Hamilton path can visit that side module at most twice, since every departure must use vertex 1 and it occurs once. Replacing the transitive pair by the triangle scales each permitted internal factor by three:
$$
(F_1,F_2)(TT_2)=(1,2),\qquad (F_1,F_2)(C_3)=(3,6).
$$
Consequently every itinerary count triples exactly. This gives a physical strong selective times-three step for the prime pair, and is not in conflict with the all-visit obstruction (the entire old graph was not kept as a module).

## Sharp exhaustion of that side clock

Collapse quotient vertices 2,3 into a side marker M and retain the old source marker. Ordered itineraries may visit the source one, two or three times, and M once or twice. Their two column vectors, indexed by source run count, are
$$
(13,10,2),\qquad (6,7,2).
$$
The first column multiplies H(M), the second F_2(M). Applied to the two original prime sources, the output Hamilton difference is therefore
$$
720(2H(M)-F_2(M))=720(4H(M)-Q_M(2)).
$$
Thus equality holds exactly when M is a transitive pair or a directed triangle. Indeed F_2(M) is at least (|M|−1)H(M) by splitting each Hamilton path, excluding size at least four. The orders one through three give only the two stated cases. The exact times-three move is real, but cannot be iterated by growing this same side module. It leaves a genuinely changing module or selective repair as the next requirement.

## The selective endpoints still have one-sided full cover order

All four outputs are again equal in Hamilton count but strictly ordered in every nonempty-cover coordinate from 2 through n−2. Each exact difference is 2z²(z²−1)R(z²), with coefficients of R listed in ascending powers:

| common cut | coefficients |
|---:|---|
| 7167 | 916, 1292, 472, 42, 1 |
| 2048 | 3648, 11906, 6054, 940, 50, 1 |
| 5120 | 8784, 21466, 9818, 1398, 68, 1 |
| 6144 | 4368, 14722, 7012, 1034, 52, 1 |

The companion `replay_second_profiles.cpp` computes the exact cover counts by subset convolution and compares the polynomial at the n+1 integer nodes 0 through n. The nonnegative shifted-coefficient proof of the [first obstruction](FOUR_VISIT.md) applies verbatim, with one extra positive factor z. Thus further equality cannot be retained by wrapping any of these *complete* outputs in a common strong homogeneous context; this does not preclude further selective side-module refinements.
