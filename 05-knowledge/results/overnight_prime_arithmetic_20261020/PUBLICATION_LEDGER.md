# Publication ledger: ongoing overnight checkpoint

**Not a final close-out.** This map reconciles reviewed session-13 work through commit `0ad127453c8ad4b2e31b96bb29a9b7e3d5ae56e7`. The session baseline is `3a6d7f5702feffdfcf979f6425f221f9a3746e0d`. There are 17 verified research packets in this checkpoint; review/publication of further continuations remains active.

The publication queue and each original packet manifest were checked against the exact committed bytes at the checkpoint head: **259 payload files**, with matching SHA-256 and byte counts. That verifies this checkpoint's local publication mirror; remote push verification remains with the publisher. Source proof status and scoped independent audits control mathematical claims; a frozen author's old review-pending label is not a new retraction.

| Packet | Published commit | Files | Entry |
|---|---|---:|---|
| affine | `a593dc550aab6a9f96a5d3afc64e870a81d81ed8` | 8 | [proof / route](../affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md) |
| affine25 | `1e797475cf3f98344e42f6727402cd6b0774a47b` | 5 | [proof / route](../affine_square_field_census_20261020/AFFINE_SQUARE_FIELD_CENSUS.md) |
| collatz | `4a7883035673b83149b68c3e833be043996d98b4` | 10 | [proof / route](../collatz_local_prime_power_gates_20261020/README.md) |
| bridges | `42f30a3475dfd17f541d88dcd8c0004dd49120cd` | 10 | [proof / route](../prime_local_branch_code_lifts_20261020/README.md) |
| census | `6955aab88dbf6e92943669c72d8e8110ad441179` | 11 | [proof / route](../odd_mixed_prime_census_20261020/README.md) |
| twin-defects | `d9ddd4c83a45be62b05f395112d60fec85a9e8d6` | 16 | [proof / route](../twin_quotient_defect_weights_20261020/README.md) |
| collatz-resultants | `bb1482536ac983b8381d5bbb2b4b0610eb5dc0f0` | 9 | [proof / route](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md) |
| census-extension | `5609dd53519d59767eb12937deb820ef68e98215` | 20 | [proof / route](../mixed_prime_census_extension_20261020/README.md) |
| degree | `1b1bbf3453a59de3e72e370b939afd2d7ff21ce7` | 37 | [proof / route](../circulant_boolean_degree_20261020/BOOLEAN_DEGREE_BOUNDARY.md) |
| collatz-sqrt7 | `3ada8876bfac9875a275bf18ec6bba1c8ea59985` | 7 | [proof / route](../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md) |
| closure | `be69d6d7acba0e610c7fbf53cac5170eaa731cb8` | 19 | [proof / route](../collatz_one_swap_exclusion_20261020/ONE_SWAP_EXCLUSION.md) |
| trace-lattice | `e34b69dcad1eff43cfef2a809e06a4c4bda84a21` | 21 | [proof / route](../source_conductor_collision_multiplicity_20261020/README.md) |
| collatz-negative | `3d4d19b75654c40cbd7061f811218ba0f1f1be3f` | 17 | [proof / route](../collatz_one_swap_negative_exclusion_20261020/NEGATIVE_ONE_SWAP_EXCLUSION.md) |
| nodal-exclusion | `abb4fb7ca2e7b7a40edb9ab4b601c1b38a44448c` | 23 | [proof / route](../nodal_cylinder_period_obstruction_20261020/NODAL_CYLINDER_PERIOD_OBSTRUCTION.md) |
| double-triple | `c2e9c86d6a5b2e5d99d07b78f8185859633061ae` | 25 | [proof / route](../double_triple_symmetry_boundary_20261020/README.md) |
| collatz-two-edit | `715d0c8ce4ebcaaac3ef7ef7a4609757ba3a5d4a` | 16 | [proof / route](../collatz_two_swap_norm_barrier_20261020/TWO_EDIT_FRONTIER.md) |
| collatz-two-edit-addendum | `0a4194f3263f15b75b42c85a6678e9efe4853664` | 5 | [proof / route](../collatz_two_swap_stripped_norm_barrier_20261020/NEAR_CRITICAL_STRIPPED_NORM_BARRIER.md) |

## Status distinctions

- Proofs with all-parameter statements are labelled PROVED only within their stated scope. Their finite tests are controls, not replacements for the proofs.
- The mixed-prime subgroup classification is proved through twelve; numerical orbit counts are FINITE-EXACT. Selected degree-thirteen/fourteen probes are not a complete larger census.
- Classical Frobenius, Matveev, code and resultant inputs remain credited as CITED or classical mechanisms in their packets. No external historical-priority claim is introduced here.
- The one-swap positive-family theorem is fully published and reviewed. Its predecessors retain their initial open-question text as lineage; newer scoped pointers control current navigation. The source-conductor packet likewise has a completed separate review despite its author's preserved initial header.
- The negative-clock one-swap classification and nodal-cylinder exclusion are now published and independently reviewed, with their exact exceptions and hypotheses retained in the map. The two specified degree-thirteen double-triple deficits are now proved and published, but do not extend the complete subgroup census beyond twelve. Both two-swap norm packets are now published: the six sparse forms and finite-weight exclusion are separate from the explicit nearby-positive-clock stripped-norm root-limit theorem. Exact-critical-only and actual-content upgrades are not asserted. Later cylinder, nonprimitive mechanical two-edit and low-delta continuations lie outside this baseline's reviewed packet set; subsequent approvals do not automatically enlarge it. Concurrent THM-4724 now has author status PROVED (pre-audit), not RESERVED; its incoming claims are preserved but not independently audited or used as proved dependencies here.

## Concurrent incoming work

The new baseline also includes incoming commits `8d83026c634fb088d6526ab160fd9cfa4bb5c261`, `ac09647acff799596bdce71ebf2cb03b4a6ce43b` and `0ad127453c8ad4b2e31b96bb29a9b7e3d5ae56e7`. These are not counted among the seventeen reviewed session-13 packets or their 259 payload hashes. The new CURRENT-FRONTIER link and the full HYP-9341 incoming status/update are retained. [THM-4724](../../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md) is untouched: its author status is PROVED (pre-audit), and this navigation checkpoint does not independently audit or promote its new claims.

## Navigation-only reconciliation

The START-HERE and CURRENT-FRONTIER compact pointer is replaced in place, with no added startup lines. HYP-9340 and HYP-9341 receive scoped current statuses; the latter preserves its original proposals explicitly as historical questions. Earlier maps link forward rather than silently rewriting their original research narrative. The local-carry README's source inventory is updated only for its routing text, with the prior README fingerprint and prior inventory fingerprint retained; all mathematical source, scripts and receipts are unchanged. A new sparse-norm README routes its historical OPEN claim to the later closure without editing that frozen proof or audit. The source-conductor README points to the later generic-conductor-degree-two exclusion without altering the reviewed source or its checksums.

No theorem-ID reservation, canon theorem modification, research replay claim or final-session completion claim is made by this navigation packet.

## Prior session maps

- [Prime/symmetry map](../prime_power_bijections_session_20261020/README.md).
- [Geometry/arithmetic map](../geometric_arithmetic_session_20261020/README.md).
- [Earlier map chain](../prime_power_bijections_session_20261020/PUBLICATION_LEDGER.md#prior-session-maps).
