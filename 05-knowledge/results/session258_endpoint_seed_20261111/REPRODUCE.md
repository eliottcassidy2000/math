# Reproduce Session258

Python3 with the standard library is sufficient. From this packet directory:

    python3 -B reproduce_all.py
    python3 -B reproduce_all.py --optimized

The nine engines test new endpoint identities, two independent fixed-q methods, the physical target grammar, positive seed probes, native Zeckendorf inverse/carry interfaces, and the arithmetic phase/chamber claims. They do not re-audit the inherited corpus or classical strong-order theorem.

The strongest seed engine is `seed_decision/core_check.py`, using literal permutations and split locations. The independent nine-state engine uses subset endpoint DP and is separately checked by `audit_independent.py`. Scope and universes are in the adjacent proof notes and JSON receipts. All correctness checks use explicit exceptions and remain active under optimized Python.

`REPLAY.json` and `REPLAY_OPTIMIZED.json` record engine hashes. Three independently generated normal/optimized receipt pairs are required to be byte-identical; the other deterministic receipts are checked during session close-out. `MANIFEST.json` records the final artifact hashes.

The repository documentation checker initially found a one-line/32-byte inherited startup-budget overrun. The session's compact tournament route repairs that local issue; no historical claims were rewritten. The final required checker and all new relative links pass as recorded in `MANIFEST.json`.

## Emit a seed and its actual tournament

    python3 construct_seed.py --q 20 --tail-length 1
    python3 construct_seed.py --q 29 --tail-length 2
    python3 construct_seed.py --q 7

The first two return literal base/core/tail adjacency rows and independently counted totals. The third returns an all-order absence result for exact multiplier7 in the selected-set tail mechanism. General fixed-q input is supported by the proved finite strong-core search, but large inputs may be computationally prohibitive. Failure for one q does not exclude its odd slope at higher exponents.
