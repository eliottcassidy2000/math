# Boundary hubs, double descent, and bridge cuts

Status: **PROVED**, independently audited and accepted by the coordinating reviewer on 2026-10-20. All-order nonleaf-bridge classification remains **OPEN**.

Read the [main proof](BOUNDARY_HUB_AND_DOUBLE_DESCENT.md), [independent audit](AUDIT.md), and [exact dependency pins](DEPENDENCIES.json).

## Results and graph conventions

- A shared literal black hub, with nonzero triangle change, certifies Q cycles of length at most four.
- Every incident Q edge at the stated rigidly guarded nonstrong nodes lies on such a cycle.
- Every strong non-self-converse node has at least two distinct strictly triangle-descending Q edges.
- A terminal-free side of a Q bridge must hang from a self-converse single-descent endpoint, its unique minimum of triangle potential.

Q is the **paired-black multigraph**, retaining loops and parallel edges. K is the literal black graph and doubles every edge of Q. F is K plus the actual blue lines. The actual blue-only graph and the quotient obtained by contracting black components are different graphs. Their bridge and leaf conventions must not be interchanged.

The proof does not establish black-core connectivity, universal short-cycle coverage, or the absence of nonleaf Q bridges. The newer singleton-padded blue-bypass results do not settle the corresponding Q problem. Their finite coverage through total order eleven concerns nonstrong mixed sources, not a complete order-eleven census; see the audit's separate freshness comparison.

## Portable, read-only controls

From this directory, run:

    python -B independent_controls.py
    python -B -O independent_controls.py
    python -B independent_common_hubs.py
    python -B -O independent_common_hubs.py
    python -B replay_guard_certificates.py
    python -B -O replay_guard_certificates.py

No third-party packages are required. All scripts compare their computed result with bundled expected JSON and do not rewrite the packet. `-B` prevents Python bytecode cache files.

- [Independent quotient, two-orbit and block controls](independent_controls.py): 150 strong non-self-converse classes through n=7, 18,448 cycle-cut pair comparisons, 2,184 two-prime-orbit systems, and 279,488 signed prime-block certificates.
- [Independent shared-edge hub checks](independent_common_hubs.py): 3,777 literal alternate routes through n=7.
- [Saved guard-certificate replay](replay_guard_certificates.py): all 65 literal certificates and every Q edge at the nine listed sources through n=8. This does not independently enumerate all eligible n=8 sources.

The small input files contain precisely the original node ID, fixed-path representative and self-converse flag for every merged class through n=7. Their full-source commit and SHA-256 pins are included. This is a projection of an authenticated prior class table, not a new class enumeration.

Results: [independent controls](results/independent_controls.json), [common hubs](results/independent_common_hubs.json), [guard replay](results/guard_certificate_replay.json). The [original 65 certificates](controls/boundary_hub_certificates.json) and [historical full-atlas replay summary](controls/boundary_hub_checks.json) are retained byte-for-byte.

## Optional historical full-atlas replay

[replay_full_atlas.py](replay_full_atlas.py) preserves the original producer checker's mathematical logic, with explicit input-directory and read-only output adaptations:

    python -B replay_full_atlas.py --atlas-dir /path/to/full/atlas
    python -B -O replay_full_atlas.py --atlas-dir /path/to/full/atlas

It requires census_n3.json through census_n7.json at the supplied root and n8/census_n8.json beneath it. The full n8 atlas is omitted for compactness; its exact SHA-256 is in [DEPENDENCIES.json](DEPENDENCIES.json) and the prior [full-atlas hash inventory](../merged_metagraph_black_quotient_bridges_20261020/FULL_ATLAS_SHA256SUMS). Regenerate it using the [prior atlas packet](../merged_metagraph_black_quotient_bridges_20261020/README.md). This optional replay checks all eligible nodes and the full reported descent statistics through n=8, then compares the original receipts exactly.

## Proof provenance and packaging

The saved candidate was based on commit 7d21b9fdcfaa30e723b92a51ebe02d35ff57b250. Its original common-hub operation is credited to [the earlier blue-detour proof](../defect_transitions_20261020/blue.md). Independent review checked the exact [prime-orbit block and all-order descent proof](../black_descent_20261020c/cycles.md), including its required dependencies.

The only mathematical proof-text completion is an explicit sentence excluding ell=1 in the nonregular double-descent case. Status and completed-obligation wording were updated, and links normalized. The original saved candidate is unchanged and its hash is retained in DEPENDENCIES.json. The newer 0fccdd7bc00ab288a7dd182b0c8fc50fd72bd4ae comparison is not retroactive proof provenance.

Publication packaging changes verifier input paths to bundled representative projections and replaces result-writing with exact comparisons. The full-atlas producer replay has only path/output adaptations. The standalone 65-certificate replay is an additional packaging check of existing certificates, not a new mathematical claim. [PACKET_VALIDATION.json](PACKET_VALIDATION.json) records the portable and optional replay checks. All packet files are plaintext; the add-only publication manifest is maintained outside this payload to avoid a self-hash cycle.
