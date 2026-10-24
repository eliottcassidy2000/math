# Independent sandwich audit

**PASS.** See `PROOF_AUDIT.md` for the proof and scope verdict.
`SOURCE_AUTHENTICATION.json` authenticates all seven source inputs against the
immutable repository commit. `inputs/frozen-sandwich/` contains exactly the
closed producer packet accepted for review, including its own manifest.

Reproduction, requiring only Python's standard library:

    python3 verify_manifest.py
    python3 check_independent.py
    python3 -O check_independent.py
    python3 inputs/frozen-sandwich/verify_manifest.py
    python3 inputs/frozen-sandwich/check.py
    python3 -O inputs/frozen-sandwich/check.py

Run from any working directory using absolute script paths if desired.
All assertions are explicit runtime checks and retain their force under `-O`.
The outputs `normal.json` and `optimized.json` match byte-for-byte. The original
producer outputs are reproduced separately as `theorem-normal-replay.json` and
`theorem-optimized-replay.json`. `PORTABILITY_REPLAY.json` records closed-manifest
copies replayed in fresh `/tmp` directories with no repository dependency.

The independent implementation exhausts the native pinned menu for r=1..4 and
all75 labelled R of orders1..4, checks the literal resource/capacity/clock maps,
and independently recomputes the famous H/Aut controls through Paley11. Finite
controls support, but do not replace, the all-size mathematical proof.
