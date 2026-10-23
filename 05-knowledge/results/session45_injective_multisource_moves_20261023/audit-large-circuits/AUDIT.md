# Independent audit: large injective native source circuits

23 October 2026. Computational audit PASS. The supplied `THEOREM.md` proof has been reviewed: its rank-two criterion, complete pinned-menu atomicity, all-dyadic dependency, conditional compatible-supply theorem, binomial thresholds and DP complexity are sound within the stated scope. The revised wording has normal word/number spacing and retains the necessary scope limits.

## Independent checks

`audit.py` imports neither the producer checker nor its witness generator. It reads the original pinned atlas, owner map, colored edges and marked electrical lifts, and separately authenticates the supplied native paths and image witnesses.

Normal and optimized Python executions pass 42,981 explicit, optimization-stable checks, with byte-identical JSON receipts. Receipt SHA256:

    f4d14bd1f0cb9551d2cda390a8c81e8287e9e72f453ad6e0b5709e6852ff14ba

The audit:

- Rebuilds all twenty source tournaments from their original canonical codes and tests all 40,320 owner permutations for each source, 806,400 source-permutation pairs in total.
- Finds trivial source automorphism groups, nonselfconverse ordinary source carriers, and complete Hamiltonian inventory sizes consistent with each source's original merged fibre.
- Recounts the complete native capacities at targets 1965 and 2234 independently of the producer's selected paths.
- Filters the complete Hamiltonian inventories by all six fixed owner clocks and finds exactly the two supplied paths for every slot.
- Tests the actual indexed endpoint and internal-arc pins, injective owner occurrence maps, and retention of every fixed owner clock.
- Constructs every actual K partner and converse matrix, checks the retained directed paths, K involution and converse square, and derives all four forward words directly from those matrices and paths.
- Checks actual source/target atlas ownership, original paired-black edge availability, and all literal F/R tile identities.
- Proves the stronger pool-wide resource fact: all forty before/after tiles are pairwise disjoint, totaling 160 closed ports. Thus there is no hidden cross-state resource reuse in any of the six circuits.
- Independently multiplies all seven full 3,528-coordinate image antecedents by every original colored edge. Every supplied equation Lz = charge + 256 correction holds over the integers, including the complete correction ledger, grounded potential and zero-sum conditions.
- Constructs each actual signed edge-conductance change from before/after target owners and checks its row action against the rank-two formula, rather than assuming that formula as a native certificate.
- Exhausts all 87,360 joint menu states across the six cases, including 17,568 balanced states. Every proper balanced rejected state has a nonzero current at an actual owned target coordinate for a directly authenticated original harmonic character.

## Finite results

| Slots | All states | Balanced states | Full H8 states | Full H256 states | H8 minimax slot budget | H256 minimax slot budget | Ordinary clock cost | Closed ports per endpoint |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|6|64|20|2|2|6|6|12|24|
|8|256|70|2|2|8|8|16|32|
|10|1,024|252|8|2|4|10|20|40|
|12|4,096|924|16|2|4|12|24|48|
|14|16,384|3,432|68|2|4|14|28|56|
|16|65,536|12,870|196|2|4|16|32|64|

The H256 two-state claim is independently sufficient and necessary: the two endpoints have complete original-L image witnesses, and every other balanced state has an actual harmonic separator. It does not depend on assuming that a saved list exhausts the harmonic bank.

The exact H8 state count uses the complete paired-bank coordinate theorem independently audited in the sibling `../audit-dyadic` packet: the same original edge and electrical-receipt hashes occur in its 530,454-check receipt. The H8 minimax budgets exhaust reachability at every smaller even Hamming-distance threshold and save explicit attaining paths. All their vertices have independently checked literal current equations. Thus the sizes 10 through 16 separate genuine higher-bank atomicity from mod-eight-only connectivity: four-slot H8 moves suffice, while the full higher bank forces the entire large switch.

## Necessary scope boundaries

- The shared owner maps authenticate the declared shared endpoint and oriented-arc coordinates, plus all additional per-slot fixed clocks. They do not assert a common-host tournament amalgam on all source vertices.
- This is a complete pinned binary-menu / closed selected-packet result. There is no claim that the new resources avoid every port of an unchanged complete canonical atlas, nor that an outside compensating packet is forbidden globally.
- The original selected X phase modulo eight is retained per slot. This is not preservation of every per-port phase for every ordinary harmonic character.
- All-dyadic preservation uses the separately audited exponent-256 theorem for the original 2-primary critical group and zero-sum image classes. It cannot be obtained by treating the raw saved X vector as harmonic at arbitrary larger moduli. The eight-slot and sixteen-slot cases have exact integer X target debt -256; the other four have debt zero.
- “Primitive” means no nonempty proper balanced subset of the supplied binary switch passes the full bank. It does not assert an unconstrained matroid or integer-kernel circuit theorem, or an unconditional native construction at every size.
- Ordinary owner-clock cost is twice the number of changed slots. Counting converse or full K/converse occurrences uses different conventions.
- The conclusion preserves the original harmonic banks by inclusion. It does not assert equality of old and new kernels or critical groups, or preservation at odd moduli.

## Reproduction

Run these commands from the published `audit-large-circuits/` directory. The repository root is `../../../..` and the native packet is the sibling `../balanced-circuits/` directory.

    python3 -B audit.py --repo ../../../.. --packet ../balanced-circuits --out normal.json
    python3 -O -B audit.py --repo ../../../.. --packet ../balanced-circuits --out optimized.json
    cmp normal.json optimized.json

No Git action, checkout modification, producer-file modification, or source-atlas modification was performed.

## Independent sharp finite support boundary

`audit_maximum.py` uses a different exhaustive strategy from `max_support.py`: it recomputes all labels from the original electrical receipt, enumerates all 167,960 nine-element subsets of the twenty sources, groups them into all 2,048 charge classes, then tests all 6,853,280 unordered same-charge pairs for disjointness. It does not use the producer's omitted-pair or five-bit prune.

Exactly 1,129 disjoint nine-versus-nine charge relations exist. Their omitted source pairs and counts are (176,2945):379, (1457,3000):378 and (1480,2176):372. Every relation has a proper balanced subrelation: minimum side-size counts are 288 of size two, 832 of size three and nine of size four. Complementation justifies testing side sizes at most four for a complete nine-versus-nine relation.

The total five-bit parity of all twenty source labels is six, so no ten-versus-ten relation exists. Odd support cannot be balanced. An exhibited sixteen-support relation is separately verified primitive by exhaustive smaller balanced subset sums. Consequently sixteen is the sharp maximum primitive support within this exact one-slot-per-source twenty-option menu, not a global order-eight or all-order bound.

Both normal and optimized executions of the new upper-bound audit pass 3,437 checks with identical receipts. The generated proper-subrelation witness digest is recorded in those receipts. The main audit authenticates all paths and resources; this second checker authenticates the exhaustive finite combinatorial upper bound.

A direct full-edge check also verifies that every nonloop incidence at U=1965 has even weight. Thus the theorem's hostile ordinary character 4e_U is genuinely harmonic modulo eight, as asserted.

    python3 -B audit_maximum.py --repo ../../../.. --packet ../balanced-circuits --out maximum.normal.json
    python3 -O -B audit_maximum.py --repo ../../../.. --packet ../balanced-circuits --out maximum.optimized.json
    cmp maximum.normal.json maximum.optimized.json

Final reviewed proof-document SHA256, including the sharp support bound, Boolean boundary, margin-preserving majority hostile, precise complementary-component-pair wording and concurrent provenance: `130ad9f105e174b1689a4879177ce479fed924c075837e9cb06b832bed1ab8e2`.

## Boolean boundary and final producer replay

The revised theorem's fixed-menu maximum, binary-assembly boundary, and concurrent-carrier distinction were reviewed. The independently read incoming `assembly/BINARY_ASSEMBLY.md` explicitly excludes arbitrary higher-arity electrical guards until their reduction is authenticated, so the comparison accurately describes its interface rather than contradicting its theorem. The incoming `search/LARGER_CIRCUITS.md` supplies a different sixteen-source list with owner-clock cost 36, consistent with the revised provenance paragraph; the commit identifiers themselves are the parent's verified provenance, not independently Git-queried here.

`audit_boolean.py` recomputes the complete accepted relations from the original electrical receipt and directly tests every literal pair projection. Its normal and optimized runs pass 893 checks with identical receipts. The H8 actual/pairwise-closure sizes are 2/2, 2/2, 8/64, 16/2048, 68/16384 and 196/65536. Every H256 relation is exactly the two all-equal states. All supplied majority counterexamples are verified.

A stronger independently found sixteen-slot hostile is majority(0,1806,6500)=260. All three input masks are H8 accepted. Mask 260 switches positive source 1477 and negative source 176, so it retains target margins, but the actual harmonic character 2Y0 has nonzero target current 4 modulo eight. Thus the obstruction is not only a failure of target cardinality. Other current-specific majority hostiles are saved in the Boolean audit receipt.

For any signed binary clause, coordinatewise majority preserves satisfaction: if its majority made both literals false, at least two of the three assignments would make each literal false, forcing one input assignment to make both false. Conjunction preserves this property. Existential projection also preserves it, by taking the majority of complete satisfying witnesses including auxiliary coordinates. Hence the displayed majority failures exclude every 2-CNF representation on these visible literal bits, even with existential auxiliary variables. The audit additionally checks the binary-clause fact on all signed two-bit truth tables. Nonlinear recodings are outside this conclusion.

All six current producer executions were replayed into this audit directory, not into the producer directory: main, maximum-support, and Boolean checkers in normal and optimized Python. Their paired outputs agree byte-for-byte and agree with the producer's saved receipts. The main producer now reports 18,321 checks, including its new original-label assertions.

    python3 -B audit_boolean.py --repo ../../../.. --packet ../balanced-circuits --out boolean.normal.json
    python3 -O -B audit_boolean.py --repo ../../../.. --packet ../balanced-circuits --out boolean.optimized.json
    cmp boolean.normal.json boolean.optimized.json

The final `SOURCE_PINS.json` was also verified against the checkout: all seven original file hashes and both inherited local-artifact hashes match. Its source commit identifier is parent-supplied provenance; no Git action was used in this audit. The producer Boolean checker was replayed again after its addition of the margin-preserving majority hostile, and the new normal/optimized receipts again agree with the saved producer receipt. `MANIFEST.json` freezes the final proof, all three independent audits, all three producer checker pairs and their inputs.

## Publication manifests

`PUBLISH_MANIFEST.json` lists only this audit packet's deliverable files. It includes `MANIFEST.json` as a cross-pin record, not as a claim that the external files it authenticates belong to this audit directory. `MANIFEST.json` resolves its audit, native-packet, dyadic-audit, inherited-supply and repository references using the published relative paths. Its separately named development-path map records the original workspace locations for provenance only. The workspace-only manifest generator, standard-output logs and redundant exploratory H8 path file are not publication dependencies. Neither manifest hashes itself.
