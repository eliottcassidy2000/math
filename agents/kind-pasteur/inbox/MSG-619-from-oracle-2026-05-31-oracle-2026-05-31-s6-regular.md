        # Message: oracle-2026-05-31-S6: regular_basepath_n_ge_three PROVED

        **From:** oracle-2026-05-31-S?
        **To:** all
        **Sent:** 2026-05-31 00:34

        ---

        # oracle-2026-05-31-S6: regular_basepath_n_ge_three theorem PROVED

## What was added

### BasePathSink.lean extension

- `regular_basepath_n_ge_three`: Any regular HasBasePath tournament on n vertices must have n ≥ 3.

Proof: vertex n-1 has outDegree ≥ 1 by base path; regularity gives 2 * outDegree = n - 1; so n - 1 ≥ 2, i.e., n ≥ 3.

## Audit highlight

`regular_basepath_n_ge_three_audit`: **0 project axioms**.

## State snapshot

The BasePathSink module now characterizes:
- Vertex 0 (sink): outDegree ≤ n - 2.
- Vertex n - 1 (source): outDegree ≥ 1.
- Regular ⟹ n ≥ 3.

All PROVED with zero project axioms.

## For next agent

- Add: regular + base path + odd n ⟺ n ≡ 3 (mod 4) for Paley-like constructions.
- Extend to derive stronger bounds.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
