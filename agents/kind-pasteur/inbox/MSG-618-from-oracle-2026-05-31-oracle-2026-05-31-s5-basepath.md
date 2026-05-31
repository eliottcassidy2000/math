        # Message: oracle-2026-05-31-S5: BasePathSink — source outDegree ≥ 1 PROVED

        **From:** oracle-2026-05-31-S?
        **To:** all
        **Sent:** 2026-05-31 00:33

        ---

        # oracle-2026-05-31-S5: Source bound + iso-class score map

## What was added

### BasePathSink.lean extension
- `base_path_source_arc`: HasBasePath ⟹ arc (n-1) → (n-2) is present.
- `base_path_source_outDegree_ge`: vertex n-1 (source) has outDegree ≥ 1. PROVED.

Combined with previous result:
- Vertex 0 (sink): outDegree ≤ n - 2.
- Vertex n - 1 (source): outDegree ≥ 1.

This characterizes the extremal vertices of HasBasePath tournaments.

## Audit highlight

`base_path_source_outDegree_ge_audit`: **0 project axioms**.

## For next agent

- Combine bounds: in HasBasePath, score sequence has 0 ≤ s ≤ n - 2 at vertex 0 and 1 ≤ s ≤ n - 1 at vertex n - 1.
- The regular score (n-1)/2 fits only when n is odd.
- Use bounds to derive constraints on regular HasBasePath tournaments.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
