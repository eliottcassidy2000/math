# Independent packaging signoff

Status: PASS for the documented proof-packet entry points.

Copied only the five manifested files and SHA256SUMS into a fresh temporary directory, without bytecode or repository neighbors, then executed them from /tmp. All five SHA-256 hashes pass. Both `python -B verify.py` and `python -B -O verify.py` pass: formal dimension 210, reachable dimension 8, observable dimension 4, seven literal controls, root-start count 1. Both normal and optimized independent.py runs pass all 15 endpoint-DP controls through word length three. Python and SymPy are the only dependencies for these entry points.

REVIEW.md is included and hashed; RESULT.md explicitly lists it and independent.py. The initializer's corrected docstring and independent checker's file-relative import are present. Claims remain limited to the exact physical family and the fixed-candidate bounded-additive-correction obstruction, expressly excluding a universal compiler, tail occupancy, and F(8,137).

Nonblocking packaging caveats: literal_forest.py retains an undocumented legacy standalone __main__ test that imports a neighboring session; the documented verifier and independent checker do not execute that branch. REVIEW.md's historical references to three manifest hashes and the former docstring are contextualized by its final packaging-change section. The original release contains unmanifested generated __pycache__ files; these are not needed, and were excluded from the portability test. Prefer distributing only manifested files plus SHA256SUMS.

No release source was modified. This is packaging signoff, not a replacement mathematical audit.
