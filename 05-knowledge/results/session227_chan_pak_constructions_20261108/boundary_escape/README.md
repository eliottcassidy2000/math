# Rooted grafts escaping the source faces

The central theorem characterizes exactly which arbitrary rooted grafts lie outside the three possible terminal faces of nonempty O/Z/J words: for a gadget of at least two vertices, both marked input vertices must be nonsinks. This separates full rooted states, not Hamilton-path-count spectra.

The smallest escaping gadget is a directed edge marked at its source. It has rank three, an explicit image plane and an integer inverse on that plane. The packet also proves the root-start multiplication law, provides literal strong witnesses and a cubic recurrence family, and gives a self-contained proof of the inherited defect inequality needed for the classification.

Files:

- `PROOF.md`: precise statements, proofs, attribution, and scope.
- `verify.py`: self-contained standard-library verifier with optimization-safe checks.
- `CHECK.json`, `CHECK.optimized.json`: deterministic normal and optimized receipts.
- `SHA256SUMS`: hashes of the deliverables in this directory.

Reproduce with:

    python3 verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

No external data, package installation, or repository checkout is needed for the verifier. Independent audit status is recorded in the containing Session227 packet. Chan–Pak completeness remains open. Experimental probes outside this listed deliverable set are not part of the theorem or its verification receipt.
