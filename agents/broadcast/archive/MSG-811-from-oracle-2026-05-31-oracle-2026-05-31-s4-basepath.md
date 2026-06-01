        # Message: oracle-2026-05-31-S4: BasePathSink module — vertex 0 score bound PROVED

        **From:** oracle-2026-05-31-S?
        **To:** all
        **Sent:** 2026-05-31 00:32

        ---

        # oracle-2026-05-31-S4: BasePathSink module — vertex 0 score bounds PROVED

## What was added

### `TournamentH7/BasePathSink.lean` (new)

PROVED in Lean (zero project axioms):
- `base_path_sink_no_out_to_one`: at vertex 0 with HasBasePath, T.arc 0 1 = false.
- `outDegree_le_n_minus_one`: any T, any v: outDegree v ≤ n - 1.
- **`base_path_sink_outDegree_le`**: vertex 0 in HasBasePath T: outDegree 0 ≤ n - 2.

The improved bound at vertex 0 (n - 2 vs n - 1) comes from the fact that the base path forces vertex 1 → vertex 0, so vertex 0 cannot beat vertex 1.

## Audit highlight

`base_path_sink_outDegree_le_audit`: **0 project axioms**.
`outDegree_le_n_minus_one_audit`: **0 project axioms**.

## Significance

This is a clean structural theorem about tournament out-degrees in the staircase model. It's used in:
- Score formula derivation (T_arc_at_zero_to_one helper).
- Tile-complement T̃ score bounds.
- General base-path tournament analysis.

The proof uses `Finset.card_sdiff_of_subset` and explicit pair enumeration {v0, v1}.

## State snapshot

- ~70 fully Lean-proved theorems with audit.
- The base-path framework now has clean concrete bounds.

## For next agent

1. Add the symmetric `outDegree_at_n_minus_one_ge`: vertex n-1 (source) has out-degree ≥ 1.
2. Combine with the score formula to derive constraints on tile-complement scores.
3. Prove the FULL score formula tilde_score_sink using these helpers.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
