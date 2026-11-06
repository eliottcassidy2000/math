# Reviewed release and provenance

The independently reviewed input proofs are preserved unchanged under `originals/`:

- CONTENT.md: SHA-256 947b483f77b80013429ad3db2eb27fc20c19653fa16eccbd40f0710cc09fc64e.
- TRIANGLE_67.md: SHA-256 ed6caf6d6c0563732b32eac301a3d017c46ab14f4b6d77e53c691ea1b5e40707.

Only their status sentences changed for this reviewed release. `SCOPE.md` spells out the exact factorization independently approved in `audit/AUDIT.md`, including its corrected expanded-bank scope. The original mathematical claims required no correction.

The audit report is pinned at SHA-256 a647a13773cb6299e62fca96a8cbc5c9dc020443292c731145b828c883028064. The transitive executable dependency `verify.py` is included and pinned, addressing the audit's reproducibility recommendation. Both author scripts and the independent verifier are replayed in the release directory.

This packet is separate from the fixed-block entropy theorem. It proves no all-order target occupancy and no universal odd-count threshold.
