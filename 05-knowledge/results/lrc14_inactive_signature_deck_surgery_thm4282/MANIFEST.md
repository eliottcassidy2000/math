# THM-4282 inactive-signature surgery and carrier-descent packet

Status: **PROVED RELATIVE TO THM-4281 + FINITE-EXACT + DETACHED
LITERAL-WALL AUDITS PASS.** This packet proves an exact 850-row deletion from
the post-THM-4281 fixed-pool residual. It proves neither physical entry nor
LRC(14).

## Consequence

The proof graph is the typed union of:

```text
index-367 fibre plus index-520 surgery C     586 rows
index-256 explicit surgery S                 188 rows
augmented-carrier band K                      90 rows
C intersect S                                 9 rows
C intersect K                                 0 rows
S intersect K                                 5 rows
triple intersection                           0 rows
union                                        850 rows
```

The union has FNV `8f595510210a5785` and SHA-256
`7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b`.
Its exact complement in the 24,223-row post-THM-4281 residual has 23,373 rows,
FNV `c6ab0ae49ee32273`, SHA-256
`c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3`,
maximum endpoint 644, and seven top rows.

## Component packets

The nested packets retain their own byte ledgers:

| component | role | nested `SHA256SUMS` SHA-256 |
|---|---|---|
| `components/index367` | complete inactive-signature atlas and exact 26-row singleton-fibre surgery | `569610a0ff6b7b867baeeb958a91bba7fabdcadc97e79464fb29a6c77156dae5` |
| `components/surgery520` | exact minimum-six five-deletion surgery and detached complete response universe | `6115b35d1c60b54812cb7bbb3b5d5afa0806d0d15968bccb61c843cddfd7059c` |
| `components/carrier663` | independent reconstruction of the inherited 8,951-mask carrier and first top-row closure | `b6e6bf879739585c4430880f7f5fdc9448ff6867848a8b4c943b665bb1aea390` |
| `components/carrier645` | detached derivation and literal audit of the exact 90-row band under the 45-mask augmentation | `28eb127fe0611cf7a2e6b5ce2f7bd56c682b500a6bbe337a40bc509de895e294` |
| `components/surgery256` | solver-free explicit 422-mask deck, 188-row literal audit, and independent proof-graph replay | `f0b5bb0dc5b94d0a955c164cad4e7bde746b263c66bc4956c423339881476201` |

Nested `REPRODUCTION.md` files preserve their original staging paths as
provenance. `REPRODUCTION.md` at this root supplies the canonical repository
paths and supersedes only those path assignments, not the commands or claims.

## Top-level primary evidence

- `results/final8996_primary_audit.out` is the source-pinned primary replay on
  all 90 carrier-band rows. It contains 90 pair records, total failures zero,
  and SHA-256
  `73dbe44577fd0850e7fec453308f726e9e50bc73f03b89d677cbcf2d161e21db`.
- `results/final8996_primary_audit.err` is empty.
- `results/final8996_primary_failures.csv` contains only its header.
- `controls/reject_short_pair_ledger.err` and
  `controls/reject_short_addition_ledger.err` freeze the two hostile identity
  failures.
- `results/endpoint650_failures358.csv` is the exact stage-D obligation ledger
  consumed by the complete active-response greedy program.
- `results/endpoint650_greedy31.out` and
  `results/endpoint650_greedy31_masks.txt` freeze its explicit upper bound and
  disjoint-obligation packing lower bound.
- `results/proof_graph_consequence.out`, `deletion_union850.csv`,
  `carrier90.csv`, and `final_residual23373.csv` are canonical convenience
  copies of independently replayed nested ledgers.

## Canonical sources and ledgers

```text
04-computation/lrc14_inactive_signature_deck_surgery_thm4282/
  final8951_joint_exposure_scan.cpp
  endpoint650_large_response_greedy.cpp
  proof_graph_consequence.py
  carrier_closed_663_645_pairs.csv
  carrier_continuation_masks.txt
  carrier_through650_additions.txt
  carrier_through648_additions.txt
```

The primary band scanner pins, inside the source, the 90-row pair count/FNV,
45-addition count/FNV, and resulting 8,996-mask carrier FNV. It therefore
rejects omissions and truncated stage ledgers before scanning bodies.

## Scope and exclusions

- The exact minimum claims are only those named in THM-4282. The index-256
  eight-mask append and endpoint-650 31-mask append are explicit upper bounds.
  Exploratory SciPy/HiGHS minimum claims are excluded.
- The 216-mask inactive classifier is not a body cover or proof deck.
- A common deck and a carrier closure are not merged into a single
  certificate. Their set overlap is accounted for only in the proof graph.
- Binaries, timing logs, nondeterministic progress streams, and discarded
  greedy orderings are excluded.

The root `SHA256SUMS` is the authoritative raw-LF ledger for the complete
promoted tree and canonical computation paths.
