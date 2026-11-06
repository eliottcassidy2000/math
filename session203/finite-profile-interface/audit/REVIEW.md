# Independent audit: finite weighted presentation of bounded-start profiles

**PASS in the stated scope.** For every positive integer K, the proof gives a finite nonnegative weighted presentation of exactly the actual rooted-profile pairs (B,s) with 1<=s<=K. In particular K=26 needs at most 15 active registers in the existence construction. It does not prove regularity or semilinearity of the numerical B-support, F(8,137), or an eventual Hamilton-path spectrum threshold.

Reviewed proof SHA256: `774e5b4fc3b4a263871e3559af2b7145a7db0e888137c391c5bf795f95ca892e`. This is the stronger TOTAL marked-budget proof, not the earlier 129-register draft. The original frozen candidate and reproducible controls are preserved under `reviewed-inputs/`.

## 1. The register-width construction is complete

Write m(0)=0 and m(t)=floor(sqrt(4t-1)) for positive integers t. The published endpoint theorem implies that a marked set with total start count at most t has at most m(t) vertices. Therefore

    W(K)=1+max_(a=0,...,K) [m(a)+m(K-a)]

is a finite integer bound. The stated estimate W(K)<=1+floor(2*sqrt(2K)) is valid. For K=26, a=13 gives 1+7+7=15, and the upper estimate is also 15. The small boundary W(1)=2 is consistent with the construction.

The induction maintains a prescribed nonempty marked set S with total start count at most K, not the set of all vertices whose individual start is at most K. This distinction is essential to the improved width.

Choose r in S and let C,D1,...,Dt be the ordered SCCs of T-r. The entry set A is nonempty, and

    start_r(T)=q*sum_(a in A) start_a(C),
    q=product_i H(Di).

The required child marks are S_C=(S intersect C) union A. The budget inequality is correct even when these two sets overlap: summing the two sets separately is an upper bound, not a presumed disjoint decomposition. Rooted insertion gives start_x(C)<=start_x(T) for x in C, while the entry sum is start_r(T)/q<=start_r(T). Hence the marked child budget is at most the original budget K.

The insertion fact used here has no hidden order assumption. For an induced strong C in strong T, choose a fixed reachability order of the outside vertices, each with an incoming edge from the current set. In any root-starting Hamilton order, insert the next vertex immediately after the last path vertex that beats it. If that vertex is not last, the new vertex beats its successor. The first root never moves, and deletion recovers the old path. This gives the required injection for each retained old root.

The recursive call constructs C before r or any tail block has been introduced. Thus no uncounted external live registers are carried into that recursive call. At return, precisely S_C is retained. Introducing r can access every member of A; every forgotten C vertex necessarily beats r. After r is marked, old C vertices outside S may be forgotten because C dominates all later tail blocks.

At a tail stage, put a=start_r(T). All retained marks other than r have total start at most K-a, so there are at most m(K-a) of them. A current tail block D satisfies H(D)<=q<=a. Applying the endpoint theorem to all vertices of D gives |D|<=m(H(D))<=m(a), including a singleton. Keeping the whole current block therefore needs at most 1+m(K-a)+m(a) registers. Double-counting a marked vertex of that current block only makes this an upper bound.

Every possible backward edge of a newly introduced tail vertex goes to r or to an already introduced vertex in the same block. Those vertices remain active. Earlier blocks dominate it, so forgotten vertices never need to be selected again. Exactly S remains at completion. The procedure may use knowledge of the target graph and its final marked set; no unsupported online decision rule is asserted.

For weak rooted sources with s>0, the root lies in the initial SCC. Deleting later singleton SCCs preserves B and s. If the remaining tail multiplier is q, then s=q*s0, hence q<=K. Construct the root SCC with its single root mark and append the nontrivial tail blocks, each of order at most m(K). The root plus a whole block needs at most 1+m(K)<=W(K) registers. This proves representation of every actual bounded-start pair, not merely of a restricted source family.

## 2. Endpoint patterns are sufficient and do not overcount paths

For fixed W and p, a retained forest has at most p bulk-start components. Every other component starts at a different active vertex, so there are at most W+p components. An unordered multiset of their endpoint types, together with the active-label control, belongs to a finite set. Active vertices absent from the endpoint types are internal vertices. No hidden internal assignment is needed for future transitions.

A new vertex can be a singleton, appended, prepended, or internal between two old components. The required arcs depend only on the selected endpoints: every bulk endpoint beats the new vertex, while a new vertex can precede only an active start in its chosen out-neighbor set. All four cases are exact.

For the inverse, delete the new vertex and cut its component at that position. If the vertex was internal, this ALWAYS produces two old components, even when the old predecessor happens to beat the old successor. The transfer does not additionally insert into an unsplit old path. This convention is what makes the inverse unique and avoids a second count of the same child path.

Selecting the same old component on both sides would close a directed cycle and is correctly forbidden. Selecting two different components of the same endpoint type is allowed, with the ordered multiplicity m(m-1). Selecting one of m identical endpoint types has multiplicity m. The actual paths are distinct vertex-disjoint objects even when their endpoint types agree, so these multiplicities are correct and independent of hidden vertex identities.

An active vertex that was internal in one old forest can still become the start of the right component in another old forest. Carrying all spanning path forests, rather than only Hamilton paths, supplies exactly the required split states. There is no lost insertion at an internal active vertex.

Forgetting changes endpoint labels to bulk and regroups counts. Every actual forest had a unique old pattern, so regrouping cannot duplicate it. Reusing the label for a later new vertex is safe: all occurrences of the former endpoint label have already become bulk, and the old vertex is never selected again.

The bulk-start count cannot decrease. An introduction may merge only with a right component starting at an active port; it cannot remove a bulk start. Forgetting can only increase that count. Therefore a discarded pattern can never reenter the retained class, and truncation loses no cover that could contribute to a final cover of at most p paths.

Different legal construction words may describe the same graph or the same pair (B,s). That is harmless: the theorem is a presentation of the value set, not a unique graph census. For each individual word, the forest counts and final profile are exact.

## 3. Readouts and the finite capped control

With p=2 and the root retained, the one-component sums give h and s exactly. A one- or two-component forest with a component ending at the root gives exactly one L-cover, since the root-ending component is uniquely the designated first path. The second path may be empty.

The nonnegative B readout is correct:

- Add every one-component forest not starting at the root.
- Add every one- or two-component forest with a component ending at the root.

A Hamilton path counted in both classes correctly has coefficient two. The singleton boundary also gives B=1. This uses no subtraction in the finite weighted system.

For fixed labeled-bin cover evaluations, assigning q actual unordered nonempty paths to distinct bins gives the stated falling-factorial weight. Taking p at least the required total number of paths is always sufficient; the analogous designated-root readouts remain finite as well.

All transition and readout coefficients are nonnegative integers. For c=K+1, saturation min(c,x) commutes with nonnegative integer addition and multiplication, including zero coefficients. Thus the capped forest vector and finite register/root control form a finite deterministic recognizer for legal words whose final s is any prescribed value in 1,...,K. Keeping the uncapped vector in parallel retains the exact B value.

The initial empty graph uses an unmarked control. A root mark may be placed only on an active vertex; a marked vertex cannot be forgotten until the mark has moved to another active vertex. Final acceptance requires a marked active root. These finite-control conventions can be imposed without changing the construction: a recursive temporary root is replaced by r before old marks are forgotten. They make the empty-state and marker boundary explicit.

The capped state does not encode exact B. A regular construction language with nonnegative weighted outputs need not have regular binary output support or semilinear numerical support. The proof correctly stops at an effective weighted presentation and does not supply an all-target inverse or coverage theorem.

## 4. Dependencies and frozen scope

The endpoint and bounded-entry theorem is now independently audited and published at [f51bbf5905807957a44057a68ecb7eba3bbd8f9d](https://github.com/eliottcassidy2000/math/blob/f51bbf5905807957a44057a68ecb7eba3bbd8f9d/session203/bounded-entry-grammar/PROOF.md). Its frozen proof SHA256 is `651a45b8d3a4b3bb307b1366d2a9b8199925788e99f42d5375abec4f265c3b76`.

The weak-tail profile preservation is frozen from [the published rooted-profile order bound](https://github.com/eliottcassidy2000/math/blob/61f06512d0eafdb1c177181d4556dda2e12f2fcb/session203/rooted-profile-order-bound/PROOF.md), SHA256 `839e6b74ba3344a93fa32240b939f57bd1dcfb1fd0114bf578568085a7e6aa7b`.

The candidate's opening description of the first dependency as still under audit is historical, not its current status. Exact published copies and provenance are retained. The independent review checked the dependency arguments needed here; no reserved statement is treated as a proved premise.

The reproducible register fixture contains 281 literal graphs. Its provenance explains deterministic selection from research data; completeness of that research bank is neither assumed nor needed. The earlier 65,125-source exploratory receipt is not included in this audit's certified replay claims.

## 5. Independent controls

The independent forest checker constructs actual unordered path forests by enumerating successor assignments, rejecting indegree collisions and cycles. It then compares the entire endpoint-pattern dictionary, not just scalar final counts, against an independently written four-case transfer.

Normal and optimized Python receipts agree on:

- All 1,099 labeled tournaments through order 5.
- 44,074 exact pattern comparisons, including intermediate introductions, every active subset, repeated bulk types, internal active vertices, forgetting and label reuse.
- 96,391 saturation-commutation checks.
- 90 reused-label construction words of widths 2,3,4 through order 6, with 869 introduction/forget operations.
- 5,616 root readouts, including a separate literal Hamilton-path count of the physical edge extension to check B.

The independent register compiler uses Tarjan SCCs, chooses the largest retained root rather than the producer's smallest, introduces tail vertices in reverse label order, and uses actual reusable register labels and root-mark operations. It checks every constructed incidence and every backward-edge access. Optimized and debug C++ receipts agree on:

- 571 strong labeled inputs through order 5 and all 17,239 eligible nonempty marked sets with total start at most 26.
- 3,807 positive-rooted weak-profile normalizations through order 5, including direct preservation of both B and s.
- All 281 frozen larger fixtures through order 13, with their complete B,s profiles independently recounted.
- 63,343 budget-contraction checks and 226,722 literal register events. The largest observed use is eight registers; this finite observation does not replace the proved all-order bound of fifteen.

The producer's forest normal/optimized replays, register-fixture optimized/debug replays, and optimized marked-set replay also match their frozen receipts. The producer marked-set include renames an unused main function, which triggers a missing-return warning in that unused function. It is not called and does not affect these runs; an explicit return-zero cleanup may be separately reviewed for release. The independent checks have no such dependency.

## 6. Reproduction

From this audit directory:

    python3 check_forest.py > FOREST_NORMAL.json
    python3 -O check_forest.py > FOREST_OPTIMIZED.json
    cmp FOREST_NORMAL.json FOREST_OPTIMIZED.json
    g++ -O3 -std=c++17 check_register.cpp -o register
    ./register > REGISTER_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check_register.cpp -o register_debug
    ./register_debug > REGISTER_DEBUG.json
    cmp REGISTER_NORMAL.json REGISTER_DEBUG.json
    python3 reviewed-inputs/check_forest.py > PRODUCER_FOREST_NORMAL.json
    python3 -O reviewed-inputs/check_forest.py > PRODUCER_FOREST_OPTIMIZED.json
    cmp PRODUCER_FOREST_NORMAL.json PRODUCER_FOREST_OPTIMIZED.json
    g++ -O3 -std=c++17 reviewed-inputs/check_register.cpp -o producer_register
    ./producer_register reviewed-inputs/REGISTER_FIXTURES.tsv > PRODUCER_REGISTER_NORMAL.json
    cmp PRODUCER_REGISTER_NORMAL.json reviewed-inputs/REGISTER_FIXTURES_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 reviewed-inputs/check_register.cpp -o producer_register_debug
    ./producer_register_debug reviewed-inputs/REGISTER_FIXTURES.tsv > PRODUCER_REGISTER_DEBUG.json
    cmp PRODUCER_REGISTER_NORMAL.json PRODUCER_REGISTER_DEBUG.json
    g++ -O3 -std=c++17 reviewed-inputs/check_marked.cpp -o producer_marked
    ./producer_marked > PRODUCER_MARKED_NORMAL.json
    cmp PRODUCER_MARKED_NORMAL.json reviewed-inputs/MARKED_NORMAL.json

Independent forest source SHA256: `ef523de8a1b80480ad452f95a814024f5230c0998f9e4c5204335e17d7a67a9c`.

Identical forest receipts SHA256: `1adc17201162b5c713f13c5438406ae1a5c040cf42a8e8e9365ef0cb4ba22e3c`.

Independent register source SHA256: `111236c027c245e1564d32ab446bd7c3d2008b13b219b450919663a4e011c429`.

Identical register receipts SHA256: `47c5b9c4fd63418230096d87bbbcdcb80df9f855cef02e3dc961c8a2c3b71fce`.

The outer manifest freezes all reviewed inputs, independent sources and receipts, excluding binaries. No original proof, source fixture or shared checkout was edited during the audit.
