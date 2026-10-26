# Cross-review of the producer's final guard theorem and decoder

**Verdict: PASS.** This review follows the independent proof/source inspection and independent verification recorded in `REPORT.md`; it is not a substitute for that independent attempt.

Reviewed producer files:

- `session70-guards/REPORT.md`, SHA-256 `3159646021bf7b2c2ad6025f636274a1e6f4f60df87f1728e3f2f4b095923563`.
- `session70-guards/check.py`, SHA-256 `62332b585c5ef31dd83bcce56449683d2a509f02a05bef76a4a77ee68c6ed68d`.

## Mathematical review

1. G1's equivalence between Hamilton-word language and directed arc support is exact, including general digraphs.
2. The cited BJR hypothesis is correctly restricted to strong tournaments and nontraceable arcs. The four-component lower bound and ordered endpoint locations match the inspected primary statement.
3. G2's interior-vertex argument does not assume the damaged component remains strong. The two witnesses are distinct, and each of their incident arcs is traceable by a spanning cycle after deleting that witness.
4. G3 replaces each missing arc by a detour avoiding the deleted vertex; repeated vertices only create a walk, which is correctly sufficient for reachability.
5. G4's decoder applies the preservation lemma to actual tournament cores, not arbitrary completions. Every eligible separator gives the true orientation; contradictory directions cannot occur.
6. G5 restores the exact ordered SCC decomposition, including singleton blocks and missing long forward cross-component arcs. All Hamilton words of a block extend globally, which justifies using the local core.
7. The exact `2^M` state bound concerns named guards with observable failure. It does not rule out smaller unnamed or aggregate interfaces.
8. G6 correctly invokes labelled-word relocation connectivity, rather than only connectivity after quotienting by automorphisms. The same sequence of current-position operations has the same accumulated original-label permutation in both tournaments. Both leftward and rightward arbitrary-span relocations are available. The stated `H(T)−1` bound is finite, not asserted polynomial.
9. G7's Boolean point-indicator product works over every field, including characteristic two: at the selected state each factor is one, while any other state has a zero factor. The full **pointwise algebra** claim is valid. The raw linear-span warning is necessary and the order-six bound `6!+1=721<1024` is correct.
10. The compact Hamilton-cycle proof supplied in the report is sound: after exclusion of insertable mixed vertices, strong connectivity yields an outside arc from the side dominated by the cycle to the side dominating it, inserting two vertices.

## Independent hostile and implementation review

The order-four masks 0 and 5 have identical one-step elementary guard signatures. Their displayed two-step trace composes to `(1,3,0,2)` and has the stated failure/success outcome. This was checked independently; every smaller fixed-path universe has an injective one-step signature.

The producer's subset-endpoint dynamic program is independent of path DFS. Its transition records possible endpoints of a Hamilton path on each subset. Applying the same construction to the converse yields possible starting vertices; joining a prefix endpoint to a disjoint complementary suffix start characterizes traceable arcs exactly.

The bit-mask conventions, induced-block remapping, SCC order, and decoder orientation were reviewed. SCC ordering uses a snapshot rather than reading the in-place sorting list. All correctness gates use explicit exceptions, so optimization does not remove them. The intentional decoded-arc corruption reverses one actual pair and must fail exact reconstruction.

`cross-review-receipt.json` records **176** additional cross-implementation cases: all `X_n` through order five, and 20 seeded dense tournaments at each order six through ten. For each case, the auditor's DFS and the producer's DFS agreed on the exact core and Hamilton count; the producer's subset DP agreed with the auditor's core; and both decoders restored the same original tournament. No mismatch occurred.

The producer envelope verifier was inspected: it checks the listed file/hash universe, exact deterministic normal/optimized receipts, full subprocess replays, and intentional decoded-arc corruption in both modes. The publisher/root owns the final envelope replay and Git publication check; this review makes no claim to have performed their later push.
