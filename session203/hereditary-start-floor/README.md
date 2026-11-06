# A hereditary endpoint floor for the four-choice interface

Status: **PROVED, independently audited in the stated scope.**

For the fixed eight-vertex core (102,0,10,83,135,30,182,47), replace vertex 5 by any tournament body of order at least three. Allow arbitrary incidences between old vertex 3 and the body; inherit every other core/body incidence uniformly. Any tournament containing this interface as an induced subtournament has, at every root, either zero Hamilton paths starting there or at least 27.

The body-order-at-least-three condition is essential. Two-vertex bodies can have minimum positive start 15. The construction fixes the other core incidences and internal core arcs. Altering them is not a guaranteed escape: any retained certified kernel still enforces its floor.

The proof selects three body vertices and six fixed nonsink core vertices. All 64 possible colored nine-vertex kernels are strong and have minimum start at least 27. A root-preserving path-insertion injection carries that floor to every induced extension, including roots outside the kernel. The all-order conclusion follows from this literal induced selection and injection, not extrapolation from larger bodies.

For the four special wrapped words 1^k c 0^k, where c is 01010, 01110, 10001 or 10101, retained eleven-vertex kernels give the sharper floors 119, 127, 133 and 141 respectively, at every scale k.

F(8,137) permits only positive start counts at most 26. Consequently these interfaces cannot supply its required sources while retaining the relevant induced kernel. This is an obstruction to one specified construction route, not a counterexample to F(8,137), an unrestricted Hamilton-count exclusion, or an eventual-threshold theorem. Other interfaces and non-induced rewrites remain possible.

## Proof and exact evidence

- PROOF.md is the audited proof with only status promotion and the verified immutable citation replacement. REVIEWED_PROOF.md preserves the original, and STATUS_AND_CITATION.diff is the complete change.
- verify.py and KERNELS.json check the four special kernels and all their root counts. NORMAL.json and OPTIMIZED.json are identical receipts.
- verify_colored_floor.py and COLORED_KERNELS.json check all 64 labeled three-vertex bodies and colorings. Their two receipts agree, and a permutation count separately checks one sharp kernel.
- audit/ preserves all 15 entries of the independent audit manifest plus its unchanged SHA256SUMS. Its separately labeled construction and direct path enumeration check every kernel, all displayed minima, the order-two boundary and 104 special embeddings. The actual insertion maps are tested on all labeled host tournaments through order five: 115,211 root/kernel maps, including 55,702 roots outside the kernel, with inverse deletion checked for every injected path.
- dependencies/ freezes the published four-choice construction and its immutable commit, Git blob and SHA-256 provenance. Its large-scale count theorem is not needed: this proof fully specifies the interface and counts its own finite kernels.

## Reproduction

From this directory:

    sha256sum -c SHA256SUMS
    python3 verify.py > replay.four.json
    cmp replay.four.json NORMAL.json
    python3 -O verify.py > replay.four.optimized.json
    cmp replay.four.json replay.four.optimized.json
    python3 verify_colored_floor.py > replay.colored.json
    cmp replay.colored.json COLORED_FLOOR_NORMAL.json
    python3 -O verify_colored_floor.py > replay.colored.optimized.json
    cmp replay.colored.json replay.colored.optimized.json

From audit/:

    sha256sum -c SHA256SUMS
    python3 check.py > replay.json
    cmp replay.json NORMAL.json
    python3 -O check.py > replay.optimized.json
    cmp replay.json replay.optimized.json

All checks use explicit failures and remain active under optimization. The independent audit explains the precise graph interface, arbitrary-root quantifier and body-size boundary.
