# Complete reviewed C3 certificate packet

The source proof, both independent reviews, and the complete producer,
canonical-census and physical-response data are included in the full-data
checkpoint. All large files are preserved losslessly through
`../ARCHIVES_CP21.json`. From the parent directory restore the original
filenames before following the producer and audit reproduction instructions:

    python3 restore_receipts.py ARCHIVES_CP21.json

- [Producer proof](PROOF.md)
- [Independent complete-census review](../c3_census_audit/REVIEW.md)
- [Independent response and minimum-size review](../minimum27_audit/REVIEW.md)

The frozen producer's pending-review wording is historical. The approved
conclusion is a necessary minimum of 28 old vertices for a nontrivial full
named-response collision obtained from simultaneous XOR of original source
modules. The finite C3-child comparison is not arbitrary-child or
arbitrary-fixed-ambient prime-nine closure, unrestricted reconstruction, or
a collision construction. The earlier proof-only checkpoint disclosed that
its complete archives were still uploading; this full-data checkpoint
supplies them without dropping any retained records.
