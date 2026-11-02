# Validation

The primary verifier completed normally and under `python3 -O`; their JSON
receipts are byte-identical. It uses no `assert` statements. All load-bearing
small tables are recomputed from actual adjacency bit rows, not imported
from an earlier atlas.

A separate Session 142 auditor (`audit_session142_gap_proofs`) independently
recomputed the 5- and 7-circulant tables, all 192 marked strong4 substitutions,
and the threshold witness with its cyclic automorphism. Its normal and
optimized receipts agree. The auditor then inspected the entire proof and
reported PASS for:

* singleton-SCC handling and the prime/p-divisor exclusions;
* the p=5 and p=7 boundary cases;
* the deterministic path-merge injection;
* the elementary marked strong-four lemma (supplied during that audit);
* lifting the full actual seven-vertex module from an induced subquotient;
* sharpness and exact automorphism order seven of the 4151 witness.

Independent files are in the sibling `../independent/` directory:
`verify_symmetry.py`, `symmetry.normal.json`, `symmetry.optimized.json`.
The parent should preserve those with this packet. The manifest hashes the
frozen local proof, primary code, receipts, README, this validation note,
and source pins. It does not claim a full repository-wide or literature
novelty audit, a global spectrum proof, or validation of unrelated packets.
