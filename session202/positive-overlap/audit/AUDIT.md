# Independent audit: free edge-order sector and positive two-overlap laws

**PASS for the revised frozen proof, with scope clarified before approval.** The actual graph construction, exact sector and scalar formulas, graph-extension identities, nonproduct example and conditional arithmetic implication are correct. The proposed 6u-1 and6u+1 constructors remain **OPEN**. No eventual-coverage threshold or unsupported target existence is approved.

## Frozen inputs and clarification lineage

Approved proof: `EDGE_SECTOR.md`, SHA-256
`56fd5b839160090a6bb0dd19117c9735872655fba7c5549636a1c165193068b3`.
Approved verifier: `verify_sector.py`, SHA-256
`19d18a017e5fb3bb6187cfafbc7249d07900edcc2d292d62410359cb8df2e8f1`.
Both frozen author receipts have hash4128e6949ca71d2398691d4af15405246ab19df01b6acb1576e86071a752adfc.

The initial proof hash10c0c69b58ca42ff650047bd8dea16c894f6fd3a0a63e4bfd5bc63108f4cca24 is preserved under the author's originals/. Three prose clarifications precede this approval: the large variable input is E(u), while the second input is the fixed triangle; the extraction implication is explicitly positive existence for odd counts; and global exclusion of7/21 is a separate inherited result. No formula, graph prescription, proof argument or verifier changed. `SCOPE_CLARIFICATION.diff` preserves the exact edits. Other exploratory probe files are outside this audit.

## Physical construction and free sector

Every unordered vertex pair receives exactly one prescribed orientation. For u>=2 the initial binary digit is one, so r,c followed by the increasing private chain and then r gives a Hamilton cycle. Since r has exactly one out-neighbor, every cycle must use r->c; the remaining private path is uniquely increasing. This proves exactly one cycle modulo rotation and strongness at every u>=2. The exceptional u=1 graph is the two-vertex edge, with zero cycles and one path. The order formula2+bit_length(u-1) covers it correctly.

In a forward-order path r must immediately precede c. The empty suffix contributes one; a nonempty suffix starting at a private vertex selected by a one bit contributes exactly the indicated power of two. These disjoint cases sum to u, with no free choice of the reverse sector.

In the reverse sector, r is last. Every deletion path extends by r except the path ending at c. That exceptional path is unique precisely when the last private vertex beats c, giving the parity u mod2; the singleton deletion at u=1 is handled separately and gives the same parity value.

For the deletion recurrence, a new private vertex beaten by c is a universal sink. A new private vertex beating c either ends the path, excluding the old c-ending deletion path, or immediately precedes c. The latter case is precisely the previously counted forward sector with its marked r replaced by the new vertex. Thus D(2u)=D(u) and D(2u-1)=D(u)+u-(u mod2), with D(1)=1, proving D=g. The resulting physical sectors are exactly(u,g(u)-(u mod2)) and total u+g(u)-(u mod2).

## Literal two-overlap operations

Appending a one or zero to the canonical binary word gives exactly the specified adjacency extensions. Move A maps E(u) to E(2u) and gives total2u+g(u). Move B maps E(u) to E(2u-1) for u>=2 and gives total3u+g(u)-(u mod2)-2. Its u>=2 domain is essential: appending a zero to the empty private word at u=1 does not give the canonical E(1) graph, even though a scalar coincidence occurs. The author correctly excludes that graph identity.

The two added three-vertex pieces really are a triangle and a transitive tournament sharing exactly the original marked edge. All private cross edges are stated. These are exact physical count identities, not an application of the overlap lower bound as an equality.

The recurrence proves g(2^k+1)=2^k+1 by induction. Hence the claimed ratio to the triangle-product count is exactly u/(2u-1), tending to one half. The input E(u) grows while the triangle remains fixed. No universal half-product inequality is asserted.

## Scalar-hostile and conditional sufficiency

Since H(E(u))>=u, the list through u=13 is a complete inverse check for target13. The table is correct and omits13. A separate literal five-vertex graph has H=13 and two cycles, so this is a genuine family obstruction, not a global missing count. The physical E(39) has H=63.

The proposed graph constructors for6u-1 and6u+1, for every u>=2, are a sufficient but unproved target. Under that hypothesis their actual outputs cover exactly all odd numbers coprime to three except1,5,7. Counts1 and5 are supplied by the displayed edge and E(3). For any odd target write n=3^a b with b coprime to three. If b is covered, ordered joins with directed triangles implement multiplication by3^a. If b=7, only a=0,1 remain excluded; for a>=2, the actual E(39) graph and a-2 triangle factors give63*3^(a-2)=n.

This proves the stated conditional positive-existence implication for every odd n except7/21. It neither constructs the two extraction operations, proves their existence by finite data, nor establishes the separate global nonexistence of7/21. No bound on a graft, graph order, or operation class has been silently added to the open target.

## Independent exact controls

The new `check.py` imports no author code. It labels private vertices first, reconstructs every edge by a separate pair predicate, and counts paths by a memoized outgoing-path recursion with an explicit first-mark flag. This differs from the author's incoming subset DP. Deletion counts and c-ending deletion paths are counted independently as well.

It verifies129 indices1..129, all256 literal A/B extensions, cycles and exact deletion/sector data through order11. It separately checks the u=1 move boundary, the full13-hostile table, the unrestricted H13 witness, E(39)=63 and the nonproduct arithmetic family through k=100. A conditional factorization stress test covers49,999 allowed odd targets through100,001; these are checks of the arithmetic implication, not claimed graph constructions from the still-open extraction hypothesis.

Independent normal/optimized outputs agree. Both author modes replay exactly and match the two frozen receipts. All validation remains active under optimization. The all-order conclusions rest on the audited arguments, not finite survival.

Reproduce from this audit directory:

    python3 check.py > check-normal.json
    python3 -O check.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json NORMAL.json

The revised proof is eligible for status-only promotion in this exact scope. The extraction hypothesis and eventual odd-count coverage must remain OPEN.
