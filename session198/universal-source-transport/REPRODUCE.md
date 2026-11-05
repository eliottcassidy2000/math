# Reproduce this packet

Run `python3 check_transport.py` and `python3 -O check_transport.py`. Both use
only Python's standard library and write their deterministic certificate to
standard output. Compare each result with TRANSPORT_NORMAL.json. The two frozen
receipts agree byte-for-byte.

PROOF.md is byte-identical to the independently audited author file named
NO_UNIVERSAL_TRANSPORT.md in audit/AUDIT.md. Its all-word arguments use the finite
kernel certificate, not a bounded-depth search. The auditor checked all matrices
against the separately audited positive-source grammar.

The three named proof snapshots are in dependencies/. POSITIVE_SOURCE_GRAMMAR.md
supplies the actual O/Z/J source operations and their physical completeness;
BILINEAR_C2.md and C3_TWO_SOURCE.md give the original scalar forms and physical
constructors. These are frozen supporting inputs, not new proofs in this packet.

Scope: no universal whole-word transfer and no single fixed universal affine
radix identity across every pair of independent physical sources in the eight
specified forms. Target-dependent source replacements, restricted physical
fibers, and exact coverage of the odd integers remain open.
