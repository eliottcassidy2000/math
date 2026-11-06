# Eight short tests for all-word synchronization

**Proved and independently audited.** Two initially synchronized nontrivial positive-start rooted sources, assigned to residues 1 and 2 with fixed admissible shifts, remain synchronized under every common O/Z/J word exactly when their eight non-start coordinates agree and their fixed starts differ by two in the specified direction. Eight words of length at most three decide this all-order condition; the exact determinant is 32.

This is a criterion for a possible paired construction. It does not assert existence of a suitable physical seed pair, occupancy of its orbit, or an eventual Hamilton-path threshold. Author-only finite discovery probes are not part of this release.

- [Proof and scope](PROOF.md)
- [Independent audit](audit/REVIEW.md)
- [Published physical transfer provenance](dependencies/PROVENANCE.json)
- [Frozen four-core count readout](audit/READOUT_DEPENDENCY.md), published at [commit e5570ebddf67175234d0504af71da481747380ec](https://github.com/eliottcassidy2000/math/blob/e5570ebddf67175234d0504af71da481747380ec/session203/staggered-profile-criterion/dependencies/FREE_ADDITIVE_PROOF.md)
- [Original reviewed proof](audit/reviewed-inputs/PROOF.md)
- [Status-only promotion diff](STATUS_PROMOTION.diff)

## Reproduction

From this directory, with SymPy installed:

    python3 verify.py > /tmp/synchronized-producer.json
    cmp /tmp/synchronized-producer.json NORMAL.json

From the audit directory:

    python3 check.py > /tmp/synchronized-independent.json
    cmp /tmp/synchronized-independent.json NORMAL.json
    g++ -O3 -std=c++17 check_physical.cpp -o /tmp/synchronized-physical
    /tmp/synchronized-physical > /tmp/synchronized-physical.json
    cmp /tmp/synchronized-physical.json PHYSICAL_NORMAL.json

All 18 audit-manifest entries and the audit manifest are preserved byte-for-byte. The package manifest also covers the promoted proof and release documentation; no local binaries are included.
