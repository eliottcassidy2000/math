# Unused edges as fundamental-cycle coordinates

Let `G` be any finite simple graph with `n>1` vertices and `m` edges, and let `P` be a fixed undirected Hamilton path. Then `P` is a spanning tree with `n−1` edges. The set `D=E(G)\E(P)` contains exactly `β=m−n+1` edges. These are the chords **unused by this particular path**; they need not be unused by every Hamilton path.

For each `e={u,v}∈D`, adding `e` to `P` creates the unique simple cycle `C_e`, consisting of `e` and the unique `u`–`v` subpath of `P`. Order the graph edges with chords first and path edges second. The F₂ incidence matrix with one column for each `C_e` has block form

    [ I_β ]
    [  B  ].

Thus its columns are linearly independent. They span the cycle space as well: for any Eulerian edge set, add the fundamental cycles corresponding to its present chords. This eliminates all chords. An Eulerian subset of the remaining path has no edges, since a nonempty subgraph of a forest has a degree-one vertex. Hence the eliminated Eulerian set was a sum of the `C_e`. This proves directly that the `β` cycles form a basis. The chord bits are independent coordinates for the **cycle space**, not a claim that arbitrary chord orientations yield a Hamilton path other than the one already fixed.

Once one of the two directed orientations of `P` is fixed, each of its `n−1` edges has a forced direction, and the `β` chords each have two free orientations. Exactly `2^β` actual-edge orientations extend that directed path. The two directions of `P` have disjoint extension sets because `n>1`. Thus one undirected path contributes exactly `2^(β+1)` incidences between edge orientations and oriented Hamilton paths. If the undirected Hamilton path is unique, this is also the number of orientations with at least one directed Hamilton path. If there are several, overlaps must be accounted for; their contribution cannot just be added as a count of distinct orientations.

For G15–G17, S147 proves the unique undirected Hamilton path and `β=1`. Its single chord is `{1,3}`, and its fundamental cycle is `1–3–6–10–15–1`. Therefore exactly `2²=4` actual-edge orientations admit a directed Hamilton path at each of these three orders. The extra chord is a cycle coordinate, not a discarded edge. This uniqueness statement is not asserted at other orders.

For G23, `m=28, β=6, U=3`, so there are `3·2⁷=384` path-orientation incidences. The independent compatibility census finds no overlaps, so there are 384 successful orientations. For G25, `m=32, β=8, U=10`, giving `10·2⁹=5120` incidences. The census instead finds 4,608 orientations with exactly one directed Hamilton path and 256 with exactly two, yielding 4,864 successful orientations and `4608+2·256=5120` incidences. These histograms and the guarded G25 polynomial were already established in the pinned Session141 source identified in [the orientation report](../oriented/REPORT.md).

Reproduce with `python fundamental_cycles.py` and `python -O fundamental_cycles.py`. The check enumerates every path in G15–G17, G23, G25–G27, verifies the unique chord projection and Eulerian parity of every constructed cycle, and enumerates all `2^β` generated vectors to confirm independence. This finite check supports the general proof, which does not rely on it. Path banks and orientation counts have the independent checks described in [the orientation report](../oriented/REPORT.md).
