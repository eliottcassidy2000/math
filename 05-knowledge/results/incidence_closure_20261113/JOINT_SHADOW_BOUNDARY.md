# Full slots plus weighted arc usage still lose exponentially many classes

**Status: PROVED IN NOTES after independent proof review; see [review](audit/JOINT_SHADOW_REVIEW.md).**
The concrete 21-vertex core and weighted-shadow collision are incoming Session319
results, not new discoveries here. The extra observation is that the *same*
collision preserves the entire slot polynomial. No minimal-order claim and no
strong-output claim are made.

## The source and the jointly retained observer

Use the vertex-transitive, strong, non-self-converse tournament U on
C7 semidirect C3 from [ARC_USAGE_SHADOW.md](../session319_arc_usage_and_signed_cycles_20261113/ARC_USAGE_SHADOW.md).
Its actual Cayley connection set, tournament/translation checks, strongness and
non-self-converse neighborhood fingerprints are verified by
[cayley_certificate.py](../session319_arc_usage_and_signed_cycles_20261113/cayley_certificate.py).
Both normal and optimized root replays return the two distinct fingerprints

    [3,4,4,4,5,5,5],   [3,4,4,5,5,5,5].

Let h=H(U), and let W(u,v) count Hamilton paths using the unordered pair {u,v}
as a consecutive pair. Retain the full labelled W and the separate labelled
start and end vectors, as well as the global full slot polynomial Q. Converse
preserves W pointwise by reversing every path. Vertex transitivity gives
s_U(v)=e_U(v)=h/21 at every vertex, so it preserves both endpoint vectors
pointwise too.

## Exact corollary at every k>=1

For a binary word b of length k, let U_0=U and U_1=U^op, and take the ordinal
join, on the same labelled blocks,

    T_b = U_b1 -> U_b2 -> ... -> U_bk -> e.

There are **2^k distinct merged tournament classes at order 21k+1** with
identical joint observer `(Q,W,start,end)`.

**Proof of polynomial equality.** For a positive integer t, Q_T(t) counts
assignments of all vertices to t ordered directed-path slots, allowing empty
slots. In an ordinal join A->B, each slot consists of its A path followed by
its B path; it can never return from B to A. Restriction to the two blocks
and concatenation are inverse, so Q_(A->B)(t)=Q_A(t)Q_B(t). Equality at all
positive integers gives polynomial equality. Reversing paths and reversing
the slot order gives Q_U=Q_(U^op). Since Q_e(t)=t, every word therefore has

    Q_(T_b)(t)=t Q_U(t)^k.

**Proof for W and endpoints.** All Hamilton paths concatenate the SCC paths
in their unique ordinal order. Internal block weights are h^(k-1)W_U.
An adjacent 21-by-21 seam has constant weight h^k/441; the final seam into e
has weight h^k/21. Nonadjacent blocks have zero weight. These formulas are
independent of b, with the 21-by-21 formula used only for k>=2. Start counts
are h^k/21 on the first block and zero elsewhere; only e has nonzero end
count h^k. This proves labelled equality, not merely isomorphism of shadows.

**Proof of class separation.** Isomorphisms preserve the ordered SCC sequence.
Since U is not isomorphic to U^op, different words give different ordinary
classes. Converse puts the unique singleton SCC first, whereas every member
of this family has it last. Thus no member is isomorphic to the converse of
another member. All 2^k classes remain distinct after merging by converse.
The outputs are not strong; their strong components are exactly the displayed
blocks and terminal singleton. This proves the claimed scope. QED.

At k=1 this is an order-22 collision. It does not contradict the current complete
through-seven observer calculations. Nor does it refute a classifier using
actual *directed path sets*, native partner incidence, or an additional marked
response: the joint observer here has forgotten those coordinates. The new
finite anchor minima and this all-scale hostile address different operator
palettes. The core's Hamilton paths are not enumerated in this proof; flat
endpoints follow from an actual transitive group action.
