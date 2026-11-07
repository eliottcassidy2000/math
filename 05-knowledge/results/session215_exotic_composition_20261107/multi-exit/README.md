# Multi-exit and arbitrary-apex positive compiler

Read PROOF.md. New structural result: for any fixed apex tournament and
chain-incidence letter, slot operators commute and are idempotent, giving
an exact nonnegative Stirling expansion of arbitrary homogeneous blocks.
This proves integer clocks 0,...,k+1 for k apexes and explicit matrix
composition for arbitrary masks. Existing two-apex clocks are credited.

The cyclic two-block descent EGF is explicitly recognized as a presentation
of inherited THM-3121/3134 substitution, not a new general composition law.
It yields an explicit three-apex fourth-clock family and all-k leading
coefficient k!^2. No Chan–Pak completeness claim.

Reproduce: python3 verify.py > RESULT.replay.json; cmp RESULT.json RESULT.replay.json.
Source pins: SOURCES.json. All artifact paths are local to this packet.
