# Reproduce Session240

From the repository root, run:

    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/replay.py
    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/replay.py --optimized
    python3 agents/check_docs.py

The replay runs the producer and independent verifier for the prime gap,
matching actions, and recovered arithmetic/coding controls; checks the mixed
finite kernel; and builds/runs a separately written C++ vertex-subset DP for
the 19-vertex witness. No network or third-party Python package is needed.
C++ requires a C++17 compiler. Outputs are deterministic after removing the
independent design verifier's timing fields. Explicit failures survive -O.

The gap producer includes uncached direct path enumeration, separate from
subset DP. The independent gap audit uses integer inclusion–exclusion of
walks and a different subset recurrence. The design auditor uses a matching
constraint solver and tournament-isomorphism lifts rather than the producer's
incidence-cell refinement. Complete matching search is checked, not inferred
merely by dividing two group orders. All32 order11 circulants are controls.

Evidence boundaries: these are ordinary exact integer computations and
human-readable proofs, not Lean formalization. No global tournament
extremality, arbitrary-prime matching classification, or Chan–Pak coverage
is verified. The 19-vertex DP certifies the displayed witness; completeness
of its min-order stratum uses the orbit/module proof and the 8*32 factor
census, rather than enumerating every labelled 19-vertex tournament.
