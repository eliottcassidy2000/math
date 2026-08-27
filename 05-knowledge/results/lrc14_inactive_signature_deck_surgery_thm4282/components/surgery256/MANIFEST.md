# Compact promotion manifest

Status: **solver-free FINITE-EXACT proof candidate**. `SHA256SUMS` is the
authoritative byte ledger. It contains 34 artifacts; `SHA256SUMS` itself is
excluded from its own ledger.

## Claim surface

This packet supports only:

1. the explicit ordered delete-seven/add-eight 422-mask deck;
2. its exhaustive body cover;
3. its exact 188-row common family on the post-THM-4281 residual;
4. the detached literal audit of all 79,336 claimed common activity cells;
5. overlaps `9` and `5`, net contribution `174`, union `850`, and the exact
   23,373-row final residual with top endpoint 644.

It makes **no minimum-eight claim**. HiGHS/SciPy sources, solver transcripts,
response-class atlases, and active-mask universe dumps are noncanonical
discovery material and are absent.

## Portable sources

```text
src/verify_surgery_deck.cpp
src/response_pattern_atlas.cpp
src/surgery188_detached_literal_audit.cpp
src/proof_graph_consequence.py
src/independent_overlap_audit.rb
```

The detached literal source is self-contained and imports no project code.
The primary source uses the two packet-local dependencies:

```text
deps/04-computation/lrc14_three_round_learned_carrier_thm4266/
  cascade_pair_exhaustive_primary.cpp
  dependency_lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp
```

Their SHA-256 identities are respectively
`c4acbbc18eb0d9b3bb3105efe090ae4bb08ccc15ac8c230b4af14dec0db00627`
and
`25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d`.
The routed `response_pattern_atlas.cpp` has SHA
`f9c6dca98169a579a67ad173f02c91af38edf01dbd81152e9540cc8e38157eaa`.

## Frozen inputs

```text
inputs/joint421_masks.txt
inputs/full_signatures_primary.csv
inputs/vulnerable_bodies.csv
inputs/vulnerable_obligations_full_carrier.csv
inputs/post_thm4281_residual.csv
inputs/prior_combined_gain586.csv
inputs/prior_combined_residual23637.csv
```

The two vulnerable-body files encode the same 71 ordered obligations in the
different exact formats consumed by the primary and detached programs.

## Frozen results

Deck and activity evidence:

```text
results/surgery422_masks.txt
results/surgery_common.csv
results/primary_surgery_replay.out
results/primary_surgery_replay_UBSan.err
results/detached_literal_surgery188.out
results/detached_literal_surgery188_UBSan.err
results/detached_literal_common_rows.csv
results/detached_literal_exposed_body_replacements.csv
results/detached_literal_residual_complement24035.csv
```

Proof-graph evidence:

```text
results/prior_carrier90.csv
results/surgery_overlap_prior_gain9.csv
results/surgery_overlap_carrier5.csv
results/surgery_novel174.csv
results/deletion_union850.csv
results/final_residual23373.csv
results/proof_graph_consequence.out
results/independent_overlap_audit.out
```

The primary O0, O3+NDEBUG, and UBSan runs are byte-identical; only the
representative stdout and empty sanitizer stderr are retained. The detached
literal O0, O3+NDEBUG, and UBSan runs are likewise byte-identical. Their
independently frozen audit packet has `SHA256SUMS` SHA
`14a92e80a0babb398927e3f74b9ebc8e282165c285d80598108554d47ffc8709`.

## Documentation

```text
THEOREM.md
REPRODUCTION.md
MANIFEST.md
```

No binaries, timing files, duplicate compiler-mode outputs, repository edits,
or exploratory endpoint-650 artifacts are included.
