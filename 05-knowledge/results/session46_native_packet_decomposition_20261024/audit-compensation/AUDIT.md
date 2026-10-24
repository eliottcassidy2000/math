# Independent audit: native compensation foundation

**FINITE-EXACT / PASS.** This audit reconstructs the supplied native domain directly from the original order-eight canonical matrices, owner-map export, colored edge export and electrical vectors. It imports neither the Session 46 producer checker nor its compressed domain report. No Git operation or producer-file change was made.

## Foundation result

The full census tests all 40,320 owner permutations of each of the 28 supplied canonical source matrices: **1,128,960 actual permutations**. Every source automorphism group is trivial, and the complete Hamiltonian path count agrees with its original two-class source fibre. Unrestricted U/V capacities also agree with the inventory.

Imposing every clock at which the two advertised endpoint paths agree gives exactly **57 actual Hamiltonian paths**, not merely 56. Exactly 56 have targets U=1965 or V=2234, one of each per source. The sole extra path is at source 1923 and target 1559. Its presence is explicitly retained in the receipt. Because the baseline has zero target-1559 margin and target occupancies are nonnegative, no state preserving every target margin can use it. Thus the complete 28-source common-clock universe with the advertised margins really is binary, even though the raw pinned domain is not.

All 57 paths use different normalized black resources: 228 distinct closed ports. For every path the audit reconstructs the actual complement partner, converse source and converse partner matrices and all four literal words; it checks Hamiltonicity, the K involution, the commuting K/converse square, the actual pure-black edge and every fixed physical owner clock. Four shared owner roles and private source-owner labels are retained explicitly. Trivial actual automorphism groups mean hidden owner/kernel lifts do not create additional capacity.

The actual original electrical vectors are checked against all 438,216 edge rows. These tests verify the five Y vectors modulo four and the primitive X modulo 256, plus every selected order-256 separator. Completeness of the entire nodewise dyadic bank is the independently audited Session 45 dyadic-filtration theorem, not an inference that five randomly chosen characters suffice. A fresh normal/optimized re-run of that separate audit is recorded alongside this audit: each run passes 530,454 checks with byte-identical receipts.

## Full-supply obstruction

An independent collision census of all two-element source subsets finds exactly **27 disjoint balanced four-slot charge trades** among the 28 labels. All source labels are distinct, so there is no balanced two-slot trade. Every four-slot trade preserves b1638+b1974. Exactly two such trades touch either critical source:

- (1638,1841) ↔ (1974,2176);
- (1638,2176) ↔ (1974,2945).

The original 16-source before and after masks have critical sums zero and two. Therefore no route whose moves change at most four slots can connect them, even with the other twelve supplied sources acting as restored buffers and regardless of their starting bits. Since every margin-preserving move changes an even number of slots, **at least six simultaneous slots are necessary in this complete 28-source common-clock universe**. This is stronger than a lower bound in a single selected compensator menu, but not a lower bound after adding new native sources, loosening clocks, changing target margins, or adding non-unit capacity.

## Receipts and reproduction

Normal and optimized execution each pass **841 explicit, optimization-safe checks** and produce byte-identical JSON receipts. Full raw domains, source matrices, owner frames, capacities, normalized resources, target labels and the complete four-trade list are in normal.json.

From the session workspace:

    python session46-audit-compensation/audit.py --repo session35-publication/clean-repo --supply session45-multisource-trades/target_pair_inventory.json --out session46-audit-compensation/normal.json
    python -O session46-audit-compensation/audit.py --repo session35-publication/clean-repo --supply session45-multisource-trades/target_pair_inventory.json --out session46-audit-compensation/optimized.json

## Scope boundaries

These are selected partial-packet domains, not a replacement of the fully occupied canonical atlas. All-dyadic means the full original **nodewise** dyadic harmonic bank; it does not mean every ordinary per-port character, odd-primary preservation, or evaluation of one unchanged integer X at arbitrarily large moduli. The chosen ordinary per-slot X phase modulo eight is separately retained. Shared physical owner records do not imply agreement of every induced arc in a single common-host tournament. Buffer resources may be revisited at different times, but cannot be used twice in the same packet state.

Route-specific physical-current, image-witness and exact bottleneck audits are separate from this foundation result and are recorded in ROUTE_AUDIT.md.
