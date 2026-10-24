# Session51 independent parent-growth audit

The audited count is a10=9,733,056, from independently evaluated
U10=97,336,068 and structurally proved E10=5,508. Read `AUDIT.md` for scope,
algorithm independence, the production/audit canonicalization boundary, and
proof details.

Quick integrity and cross-check verification:

    python3 verify_audit.py
    python3 -O verify_audit.py

Complete isolated physical-copy replay (several minutes):

    python3 replay.py

Individual replays:

    python3 replay_independent.py
    python3 replay_producer.py
    python3 structural-audit/replay.py
    python3 mutations.py
    python3 input_controls.py
    python3 reconcile_strata.py

All compiler runs use C++17. The independent and production programs are tested
at O2 and O3/NDEBUG, with full refinement in the latter trace passes. All checks
use explicit exceptions, not removable assertions. Replay outputs are generated
locally and are not a reason to rewrite the frozen input manifests.
