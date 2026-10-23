# Independent audit: converse repair lowers the alphabet minimum from twelve to nine

23 October 2026. **PASS.** This addendum does not alter the earlier frozen rule audit.

## Exact scope of the new distinction

At order eight, for the two-coordinate alphabet

    (neighbor cyclic-triangle count mod p, edge cyclic-triangle count mod q),

there are two different optimal alphabet sizes:

- **Twelve symbols** are necessary and sufficient for ordinary class separation, and for uniqueness of the stated ordinary endpoint/color/signature transportation table without reflection.
- **Nine symbols** are necessary and sufficient for uniqueness of that entire native ordinary-class table when its authenticated converse-reflection equalities are imposed.

Both minima are within this rectangular modular family, with positive integer p and q. Neither is a minimum over all possible encodings. Neither says all numerically larger moduli succeed. They do not imply uniqueness of the matching of individual ports, all-order signature completeness, or that a feasible alternate matching is the literal tournament-complement matching.

## 1. All-size converse-repair proof

The proof in `CONVERSE_MARGINAL_REPAIR.md` is sound. For a fixed source/endpoint/color category, let x_c be the unknown count routed to ordinary target class c. The retained S-row margins give sums grouped by S(c). Native reflection identifies these cells with the reflected source category, whose retained margins give sums grouped by S(c*).

One must not silently replace these two families of margins by joint (S(c),S(c*)) margins. The correct incidence object is the bipartite multigraph with one labelled edge for each class c, from S(c) to S(c*). Parallel edges remain distinct.

If this graph is a forest, its incidence columns are independent after changing the signs on one bipartition, so the two marginal families force every edge count. Equivalently, leaf elimination forces every zero-margin difference to vanish. Deleting incompatible target edges preserves this sufficient condition. Capacities and nonnegativity may also force uniqueness when the graph has cycles, so the forest test is not necessary in general.

If each class or its converse belongs to a singleton S-fiber, each edge already has a degree-one endpoint and is read directly from one margin. This stronger simple criterion applies to the actual E_(3,3) data.

The paired-label injectivity warning is also correct: the four different labels (0,0),(0,1),(1,0),(1,1) admit the nonzero zero-margin change (+1,-1,-1,+1). Distinct pairs alone do not force their separate margins to determine the counts. The independent checker contains this control.

## 2. Independent complete E_(3,3) certificate

The addendum checker reconstructs all eight-vertex tournaments from the canonical codes and recounts cyclic triples directly, without importing the producer's signature functions. It computes E_(3,3), E_(4,3), and the exact-edge neighbor-modulo-three refinement on all 6,880 classes.

Findings:

- E_(4,3) separates all 6,880 ordinary classes.
- E_(3,3) has exactly 6,876 signature fibers and uses all nine possible symbols.
- Its partition is identical to that of the neighbor-modulo-three signature with exact edge counts, not merely equal in size.
- Its four collision pairs have canonical codes (4328706,4330626), (25694272,25694466), (27467792,27467812), and (50676048,50676504).
- The converse of each of these eight classes belongs to a globally singleton E_(3,3) fiber. The receipt contains each class ID, canonical code, converse ID, converse code, and verified converse fiber size.
- The full signature multigraph has 13,752 vertices, 6,880 edges, 6,872 components, maximum degree two, and no cycles. Thus it consists of eight two-edge stars and 6,864 isolated single-edge components.

Consequently the 6,872 direct singleton classes and eight classes recovered through their converses account for every target class. All native ordinary table cells are forced under reflection, despite ordinary signature collisions.

## 3. Independent verification of the counterfeits

The checker reads the producer's explicit edits, but recomputes their acceptance conditions rather than trusting reported booleans. Against the full native ports and endpoint capacities, it verifies:

- Every reported before/after count and canonical class identifier.
- Nonnegativity, existence of source signature rows and open target columns.
- Every retained source-signature margin and every target endpoint/color margin.
- Exact directed and undirected merged edge changes, including nonzero projection.
- Unchanged blue table, literal merged degrees, and authenticated cyclic-triangle current.
- The endpoint triangle identity for every changed cell.
- Distinct literal closed-mask witnesses for every removed unit, with exact source/target classes, endpoint scores, and mirror classes.
- Full native reflection for each symmetric example, even directed black multiplicities, and even black counts in reflection-fixed native cells.

The symmetric J4, neighbor-modulo-two/exact-edge, and exact-neighbor/edge-parity witnesses each change eight native cells by one unit, and change four directed merged black multiplicities by two. Their undirected projections also change. The neighbor-modulo-three/exact-edge example is a valid four-cell ordinary counterfeit, but its naive mirror extension violates retained signature margins. Its absence from the symmetric examples is not treated as empirical evidence of uniqueness; the complete singleton-converse proof supplies that conclusion.

The stronger literal matching certificates are independently checked as well. Each symmetric example replaces the partners of four distinct actual closed ports, permutes exactly their four original open partners, respects endpoint and target-signature compatibility, commutes with reflection, and induces exactly the reported eight-cell table change. They therefore describe equivariant alternate matchings of real ports, not only formal degree changes.

These alternate matchings are not asserted to be realizable by literal tournament complementation. The retained bucket constraints deliberately omit the complete authenticated complement decoder.

## 4. Why the two minima follow

The reflected neighbor-modulo-two/exact-edge counterfeit remains accepted after reducing the exact edge coordinate modulo any q. It rules out p=2 for every q, and therefore also p=1. The reflected exact-neighbor/edge-parity counterfeit remains accepted after reducing the exact neighbor count modulo any p. It rules out q=2 for every p, and therefore also q=1. Reflected uniqueness thus requires p>=3 and q>=3; E_(3,3) attains the resulting minimum pq=9.

Without reflection, the additional neighbor-modulo-three/exact-edge counterfeit rules out p=3 for every q. Thus p>=4 and q>=3 are necessary; the independently verified class separation of E_(4,3) forces its target class in every row and attains pq=12.

The twelve-symbol class-separation statement remains correct. It was never, by itself, a proof of a twelve-symbol minimum after adding reflection constraints.

## Reproduction and provenance

The checker is `check_reflection_counterfeits.py`; its frozen output is `reflection_audit_receipt.json`. It imports no producer checker, uses only the Python standard library, and compares the complete frozen result by default. It reads:

- `../wiring/classes.tsv.gz`
- `../wiring/ports.tsv.gz`
- `../wiring/endpoint_capacities.tsv.gz`
- `../wiring/class_of_mask.u16le.gz`
- `../theory/residual_counterfeits.json`

Run in the published layout:

    python rule-audit/check_reflection_counterfeits.py
    python -O rule-audit/check_reflection_counterfeits.py

For the scratch layout, add `--wiring session38-wiring --theory session38-root`. The receipt hashes all five inputs by basename, avoiding publication-path normalization. All tests use explicit exceptions and survive Python optimization. `--write` intentionally regenerates the receipt.

The large native ledger and canonical class universe have their own independent full reconstruction audit; this addendum does not duplicate that full enumeration. Its fresh signature computation, direct native-table checks, literal mask checks, and port-bijection checks independently validate the new finite conclusion against those authenticated inputs.
