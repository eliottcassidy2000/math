# THM-4281 curated promotion manifest

This manifest records the promoted proof packet for the theorem
“Global frozen 421-mask joint-deck collapse, carrier closure through 664, and
residual top 663,” whose status is **PROVED RELATIVE TO THM-4276/4277 +
FINITE-EXACT + DETACHED NEWCOMER-LITERAL-WALL AUDIT PASS**.  The packet has 15
canonical source files under
`04-computation/lrc14_joint421_global_common_carrier_thm4281/` and 38 canonical
result/input artifacts in this directory, together with this manifest, the
reproduction guide, and `SHA256SUMS`.  The canonical theorem lives separately
at the path below.  No binary is canonical.

## Canonical destinations

Replace the existing reserved THM-4281 stub at

```text
01-canon/theorems/THM-4281-rectangle-common-joint-deck-endpoint-670-bridge.md
```

with the audited theorem text.  Copy every `src/` file byte-for-byte to

```text
04-computation/lrc14_joint421_global_common_carrier_thm4281/
```

and every `results/` file, plus `MANIFEST.md`, `REPRODUCTION.md`, and
`SHA256SUMS`, byte-for-byte to

```text
05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/
```

The staging packet's `THEOREM.md` was only an assembly draft.  The canonical
theorem file named above is authoritative and is not duplicated here.

## Included sources

```text
endpoint664_full_response_atlas.cpp
endpoint665_single_response.cpp
endpoint667_full_response_atlas.cpp
endpoint669_full_response_atlas.cpp
joint420_global_literal_activity_scout.cpp
joint421_redundant_role.cpp
joint421_response_fibre_census.cpp
proof_graph_consequence.py
response_pattern_atlas.cpp
scan_augmented_carrier_endpoint668.cpp
scan_augmented_carrier_endpoint669.cpp
scout_joint421_post4277_primary.cpp
verify_common_complement_carrier_overlap.py
verify_joint421_literal_post4277.cpp
verify_joint421_literal_r670.cpp
```

These respectively freeze the four exact response certificates, shared
response enumerator, primary endpoint layers, primary global atom scan,
detached global literal partition, detached endpoint literal controls, the
420-core global replay and mask-role census, and two independent
proof-graph/set replays.

## Included results and exact CSV inputs

```text
endpoint664_full_response_atlas_O3.out
endpoint665_single_response_O3.out
endpoint667_full_response_atlas_O3.out
endpoint669_full_response_atlas_joint421_O3.out
joint420_core_without_003c900c_masks.txt
joint420_global_literal_activity_scout.err
joint420_global_literal_activity_scout.out
joint421_common_carrier_overlap_layers.csv
joint421_redundant_role_O3_NDEBUG.out
joint421_redundant_role_UBSan.err
joint421_response_fibre_census_O3.out
joint421_response_fibre_census_UBSan.err
joint_iter39_shrunk_forward_masks.txt
post_thm4271_residual.csv
post_thm4277_joint420_common_literal.csv
post_thm4277_joint421_common_primary.csv
post_thm4277_joint421_common_primary_layers.csv
post_thm4277_joint421_literal_hostile_controls_O3_NDEBUG.csv
post_thm4277_joint421_not_common_primary.csv
post_thm4277_joint421_not_common_primary_layers.csv
post_thm4277_joint421_primary_hostile_controls.csv
post_thm4281_residual.csv
proof_graph_consequence.out
scan_augmented_carrier_endpoint664_final_O3.out
scan_augmented_carrier_endpoint664_final_UBSan.err
scan_augmented_carrier_endpoint664_joint421_O3.out
scan_augmented_carrier_endpoint665_joint421_O3.out
scan_augmented_carrier_endpoint666_joint421_O3.out
scan_augmented_carrier_endpoint667_joint421_O3.out
scan_augmented_carrier_endpoint668_joint421_O3_NDEBUG.out
scan_augmented_carrier_endpoint669_joint421_O3.out
scout_joint421_post4277_primary_O3.out
thm4281_carrier_not_common.csv
verify_common_complement_carrier_overlap.out
verify_joint421_literal_post4277_O3_NDEBUG.out
verify_joint421_literal_post4277_ubsan.err
verify_joint421_literal_r670_O3_NDEBUG.out
verify_joint421_literal_r670_UBSan.err
```

The large CSVs are retained only where they are an exact input, a claimed
universe/set, a hostile control matrix, a genuinely useful layer ledger, the
36-edge proof-graph novelty, or the final residual.  Regenerable intermediate
set copies are omitted.

## External canonical source dependencies

The packet deliberately imports four already promoted source dependencies:

| Canonical path | SHA-256 |
|---|---|
| `04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp` | `d2388ad0901ea23d9ba2360c1dd4fd8ec2d1055b8cbd3d5a86125add5b9eff14` |
| `04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp` | `c4acbbc18eb0d9b3bb3105efe090ae4bb08ccc15ac8c230b4af14dec0db00627` |
| `04-computation/lrc14_three_round_learned_carrier_thm4266/dependency_lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp` | `25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d` |
| `04-computation/lrc14_three_round_learned_carrier_thm4266/cleanroom_resistance_controls.cpp` | `1be1f341dc8ffc1cb9bc0c55f597ac71220b46aa4b7c347242adacb39b5cb21e` |

The clean-room source itself includes the promoted direct-wall implementation

```text
04-computation/lrc14_endpoint_cascade_direct_wall_body_audit.cpp
SHA-256 bbfc55a3181882cf6456900951658f95634680fb30f68ad763ffa830812c66e7
```

and the cascade source includes the listed THM-4188 dependency.  The required
prefix/replay transcripts are inherited through THM-4266 and THM-4271 and are
named exactly in `REPRODUCTION.md`.

## Byte-identical build families

```text
global literal O0/O3/O3+NDEBUG/O1+UBSan stdout:
  3081efc81b43cffe912d38f22e8b2e09719980cd0291f42551ee386a2b3100b2
global literal hostile CSV in all four modes:
  7d2e0de7ed1650edf0ca7746bc870e82f83ab0168302fd22eb1e650450d230ef
endpoint literal O0/O3/O3+NDEBUG/UBSan stdout:
  88fa48d3130dbe1a1948a22c79fa5f6c7fa3003a20cbd1c4a365b99d71b0a530
final endpoint-664 O3/O3+NDEBUG/UBSan stdout:
  07244f273ff2cc026d8038fb7cee58f6a138aa2ffe69ef757447f47864baa34c
421 response-fibre census O0/O3/O1+UBSan stdout:
  fe38cb900033085c0ad0584b730d6baaa840bcba21fc7e7462c3a4b41dad3f36
420-core body/endpoint-role O0/O3+NDEBUG/O1+UBSan stdout:
  1b463f1e7f0a5eae462f92ef87de0e8716b9aa3efead1e1ab3cdc3f253330efd
sanitizer diagnostic stderr (all named controls):
  e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
normal/Python-optimized overlap replay stdout:
  8ed72b405f3dbcbfa0355544647ea3c8f15f93eb30bc8dc872096a8df5dea20c
```

Only one semantic representative per byte-identical family is retained,
along with a separately named empty sanitizer control for each proof role.

## Serialization

FNV-1a-64 mask ledgers serialize each ordered mask as one little-endian u64.
Edge ledgers serialize `q` then `r` as little-endian u64 in inherited order.
Edge SHA-256 values use raw ASCII `q,r\n` bytes.  Typed hostile/layer ledgers
serialize their documented integer columns as little-endian u64 words.

## Explicit exclusions

The following are not in this promotion packet:

- binaries, dSYM directories, timings, progress stderr, and redundant build
  copies;
- CP-SAT/MILP/greedy/CEGAR discovery traces and temporary decks;
- the 515-mask fallback deck and the obsolete 515-hardcoded
  `verify_joint_deck_certificate.py`;
- the failed order-sensitive primary 420 scout and its zero-byte transcript;
  the packet instead freezes the successful detached-literal 420 replay and
  still makes no global deck-minimality claim;
- the obsolete temporary import shim used while making the endpoint-role
  checker portable;
- the interrupted endpoint-663 scan and every statement that would depend on
  it;
- the simpler preliminary 421 global literal scout superseded by the complete
  common-matrix plus hostile-partition verifier.

The authoritative byte ledger is `SHA256SUMS`.  Regenerate it only after all
packet files are final; do not edit a file after that freeze.
