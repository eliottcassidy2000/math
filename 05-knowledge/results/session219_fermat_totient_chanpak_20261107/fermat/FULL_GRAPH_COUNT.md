# Independent exact count of the complete 36-vertex escape graph

**PROVED, exact subset-DP certificate on16 outside vertices.** This is an independent path-count validation of the full36-vertex adjacency witness, not merely a small-prefix regression or a replay of the source-profile recurrence.

Let A be the15-vertex literal graph for prefix OOZJOZJJOJJOZ. Insert a transitive block B of m>=1 O vertices, then add the closing vertex c. Every vertex of A points to every vertex of B, and every vertex of B points to c. Edges within B point in increasing birth order. The remaining16-vertex tournament on A union{c} is read directly from the literal graph rule.

In a Hamilton path, each maximal nonempty segment of B vertices is increasing. A segment can end only at c or at the end of the entire path. It can begin only after a vertex of A or at the start of the entire path. Hence there are at most two such segments. Removing them leaves an ordered traversal of all16 outside vertices; consecutive outside vertices need not be connected by a native edge when a removed B segment supplied the bridge.

Define B_k as the number of outside traversals with k marked nonempty B slots, for k=0,1,2. These counts depend only on the outside graph:

- Native transitions a->b have weight1.
- A transition from any a in A to c has the additional possible weight z for a B slot. It may be allowed this way even if the native edge is absent.
- A traversal beginning at c can have an initial slot, weight z.
- A traversal ending at a vertex of A can have a final slot, weight z.

A standard Hamilton-path subset DP over these16 vertices, carrying polynomial coefficients through degree2, counts B_0,B_1,B_2 exactly. `check_full_36_dp.cpp` constructs the actual outside adjacency and implements this weighted DP, without calling O/Z/J profile formulas.

For a fixed outside traversal with one nonempty slot, all m B vertices go into that slot in their unique increasing order: one assignment. With two slots, every B vertex chooses slot1 or slot2, with both choices used. There are2^m−2 assignments, and each chosen slot is then traversed in increasing order. These assignments are in bijection with the full paths and never impose an incorrect order between distinct slots.

Therefore

    H_m = B_1 + (2^m−2) B_2,    m>=1.

The independently computed coefficients are

    B_0=10063,   B_1=16935,   B_2=6013.

Thus H_m=4909+6013*2^m, including H_20=6305092397 for the full36-vertex graph and all m=20+64t for the infinite escape family. This recovers the profile-derived formula by an independent combinatorial method.

The separate72-vertex five-run witness remains certified by the inherited physical recurrence plus literal adjacency generation; no full72-vertex DP is claimed. The16-vertex weighted DP is small enough to rerun quickly.

    g++ -O2 check_full_36_dp.cpp -o /tmp/session219_full_36
    /tmp/session219_full_36
