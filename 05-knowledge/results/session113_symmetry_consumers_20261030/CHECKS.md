# Publication verification and reproduction

All Python verification gates remain active under `python -O`. Dependencies: Python 3; SymPy for the Jacobian checker; NumPy for the Collatz vector checker; a C++ compiler for regenerating the Collatz vectors. No vector binaries are stored in the repository.

From this packet directory:

```sh
python lrc/check_boundary_germ.py
python composite/check_composite.py
python jacobian/check_symmetry.py
python audit/check_interfaces.py .
(cd prime_ranking && python check.py && python check_larger.py)
python audit_prime/independent.py
python audit_prime/thresholds.py
python collatz/word_controls.py
python collatz/envelope_obstruction.py
python uniform_content/check_obstruction.py
python uniform_content/check_phase_refinement.py
python audit_collatz/check_new_boundaries.py
```

The word and boundary scripts regenerate their local JSON receipts. The independent prime audit enumerates all 692,224 named-group orbits at p=5,k=2. The larger prime tests check total counts, not large rank/unrank queries.

For each fixed precision E=6,7,8:

```sh
g++ -O2 collatz/budget.cpp -o /tmp/collatz_budget
/tmp/collatz_budget 6 64 /tmp/collatz_V6.bin
python collatz/audit.py 6 /tmp/collatz_V6.bin
python audit_collatz/check.py 6 /tmp/collatz_V6.bin
```

Replace both occurrences of 6 with 7 or 8. Regeneration is deterministic for these inputs; both auditors compare the same SHA-256 certificate hash. The scalar independent verifier checks 62 million edges in total, plus entry constraints, both infinite-letter tails, and 262,143 exact positive words at each precision. Its `seconds` field is timing metadata, not a mathematical receipt field.

Publication reruns use both normal Python and `python -O`; deterministic fields and certificate hashes must agree with the frozen receipts. Documentation, relative links, and whitespace are checked separately. The mathematical proofs have the scopes in each note; finite controls are not evidence of a major-conjecture closure.

## Repaired E9 certificate

```sh
g++ -O2 uniform_content/repaired.cpp -o /tmp/repaired
/tmp/repaired 9 /tmp/repairedV9.bin
python audit_collatz/audit_repaired.py 9 /tmp/repairedV9.bin --output /tmp/repaired9_normal.json
python -O audit_collatz/audit_repaired.py 9 /tmp/repairedV9.bin --output /tmp/repaired9_optimized.json
for E in 1 2 3; do /tmp/repaired "$E" "/tmp/repairedV$E.bin"; done
python audit_collatz/check_repaired_scalar.py
python -O audit_collatz/check_repaired_scalar.py
```

The E9 vector is 78,125,000 bytes and is deliberately not committed. Its SHA-256 must equal `83e016db40ab86a1ae8d7e6aef605a4710c0d1279f23f1741e11f617a187670f`. Compare all receipt fields except elapsed `seconds` with `audit_collatz/repaired9_independent.json`. The producer is restricted to E=1,...,9 and writes explicit little-endian signed64 slots. Normal and optimized publication audits independently rechecked all 525,000,000 edges and 126 entries; the scalar controls include rejection of a corrupted vector under optimization. E10 and all-E feasibility remain open.
