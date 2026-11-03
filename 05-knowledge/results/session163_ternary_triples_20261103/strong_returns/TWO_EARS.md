# Two-vertex strong returns and exact offset obstruction

**Status.** The two graphs below continue both strong-source paths from [the first-return note](NOTE.md) by a further literal 0 branch, with the full marked observer scaled by 3. The exhaustive finite calculation excludes the two offset counts only for extensions of those two specified source graphs by at most two vertices without changing their internal edges. It does not prove a closed ternary bank or a tournament-spectrum classification. No novelty claim against the literature is made.

Recall `R2=(285,123,225)` and `R4=(291,120,225)`, with 10 vertices and marked port 9, from [NOTE.md](NOTE.md). Append vertices 10 and 11, orient 10→11, and keep 9 as the marked port. For each new vertex, the given 10-bit mask has bit i set precisely when that new vertex beats old vertex i. Take masks `(1017,1017)` for R2, and `(1007,1023)` for R4. The complete resulting adjacency row masks are

    S2=(1008,4089,4083,997,776,912,816,592,128,256,3065,1017),
    S4=( 996, 993, 994,999,1999,912,800,576,128,256,3055,1023).

Each contains its named R source literally as the first ten vertices, is strongly connected, and keeps its original port 9. Direct endpoint/subset Hamilton-path dynamic programming gives

    S2: (H,H_a,H_aa)=(855,1593,2943), observer=(855,369,675)=3R2;
    S4: (H,H_a,H_aa)=(873,1593,2943), observer=(873,360,675)=3R4.

Here `a` adds a singleton that beats only the marked port and becomes the new port. Thus, including the earlier sources, the explicit strong count paths are

    31 → 95 → 285 → 855,    and    31 → 97 → 291 → 873.

Because both new full states equal three times their predecessors, the inherited b/q2 source transfer carries each equality through any common word W in `{b,q2}`. The previous +2 and +4 relations from the 31-source remain valid when W has an even number of b's. Every source in this construction is strong: the two seed extensions are checked directly, and the b and q2 joins preserve strong connectivity as shown in [NOTE.md](NOTE.md). The new marked port is an old vertex, not the newest of the two; for later joins attach relative to that marked port and mark the last vertex of the new block. This is a source-state return, not merely equality of scalar counts.

## Exhaustive two-vertex obstruction to the offsets

No graph obtained by adding one or two vertices to R2, while leaving its ten internal vertices and their orientations unchanged, has count 857 or 859. Likewise no such extension of R4 has count 875 or 877. The statement is independent of strongness and marked-port choice.

For the two-vertex search, swapping the new labels makes 10→11 without loss of generality. Each new vertex has one of 2^10 possible orientations toward the old source. Crucially, the odd-cycle-packing identity makes Hamilton-path count monotone under taking induced subtournaments: deleting vertices can only remove nonnegative packing terms. Hence a candidate with total h≤859 (R2) must have each single-new-vertex induced extension at count at most 859; similarly the R4 bound is 877. Exact subset DP checks all 1,024 one-vertex masks per source. The surviving sets have sizes 32 and 46, so checking their Cartesian squares requires only 1,024 and 2,116 pairs, respectively. This pruning is logically complete, with no assumption about strongness or about ports.

The exact pair counts at the three targets are

    source  target counts          all two-vertex extensions
    R2      855, 857, 859          7, 0, 0
    R4      873, 875, 877          6, 0, 0.

These counts fix the 10→11 orientation and retain labelled old vertices; they do not count isomorphism classes. All 13 graphs at the zero-branch scalar counts are strong. Testing all twelve possible marked ports gives five R2 and three R4 graph-port pairs with the entire state exactly `3R2` or `3R4`. The displayed S2 and S4 are one example each. A single-vertex extension achieving one of the excluded counts would yield a two-vertex extension with the same count by adding a source or sink, so the same conclusion covers both one and two added vertices.

## Reproduction

Run `g++ -std=c++17 -O3 -Wall -Wextra two_ears_check.cpp -o two_ears_check && ./two_ears_check`, optionally also compiling with `-DNDEBUG`. The checker uses explicit failure tests, enumerates all 1,024 single masks, performs the pruned exhaustive pair census, counts Hamilton paths independently by endpoint/subset DP, tests strongness, and checks both observer queries for all ports of every scalar zero-branch candidate. The output prints all eight marked matches. The all-suffix consequences use the inherited general source-transfer identity, not a finite test of words.
