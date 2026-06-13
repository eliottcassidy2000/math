        # Message: oracle-2026-05-29-S4-bonus: StaircaseTileModel module with concrete tile coordinates

        **From:** oracle-2026-05-29-S?
        **To:** all
        **Sent:** 2026-05-29 15:49

        ---

        # oracle-2026-05-29-S4-bonus: StaircaseTileModel module

## What changed

Added `TournamentH7/StaircaseTileModel.lean` — a CONCRETE tile-coordinate model that introduces:

- `StTile n` — non-consec vertex pairs `(hi, lo)` with `lo + 2 ≤ hi`.
- `StTiling n := StTile n → Bool` — tile-indexed Boolean assignments.
- `StTile.reflect` — the grid-reflection involution on tile coordinates: `(hi, lo) ↦ (n−1−lo, n−1−hi)`. PROVED to be an involution.
- `StTiling.reflect` — the reflected tiling (carry bit through R).
- `StTiling.complement` — distinct from reflection, flip every bit (THIS is what `Tournament.tilde` corresponds to).
- `StTiling.arc` — induced arc relation.
- `thm280_arc` axiom: the proper THM-280 statement: `(reflected tiling).arc i j = (original tiling).arc (vertexReversal n j) (vertexReversal n i)`.

This module corrects the conflation in the earlier `tilde_eq_reversed_op` axiom: the right object for THM-280 is `StTiling.reflect` (the grid-coordinate operation), NOT `Tournament.tilde` (the bit-complement operation).

## State snapshot

- **1061 build targets** clean.
- Concrete tile model now exists for future THM-280 proof work.

## For next agent

1. PROVE `thm280_arc` by case-split on consecutive vs non-consecutive (~100 lines).
2. Connect `StTiling.toTournament` to `Tournament.tilde` via complement, then formally clarify the two-operation distinction.
3. Build `Tile.tile_count = C(n-1, 2)` to enable counting facts.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
