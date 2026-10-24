# Validation record

The independent checker and both unchanged inherited authentication scripts pass in normal and optimized Python with byte-identical output in each pair. See VALIDATION_RECEIPT.json for the recorded executable results and timestamp.

- Independent threshold/domain audit:1,521,278 explicit checks, including all474,896 native endpoint pairs and393,210 arbitrary-domain joint-budget checks.
- Full native menu/kernel/owner/capacity authentication:115,250 checks.
- Full inherited dyadic-filtration replay:441,763 checks.

FIBRE_SUMMARY.tsv has23,444 exact nonempty-fibre rows. AUDIT_RESULT.json contains census, grade histograms and exact lower controls. Complete pair results are streamed into the deterministic SHA256 recorded there, rather than duplicated in a bulky pair archive.

Historical source snapshots are text archives. All native inputs/authentication files are byte-preserved and independently hash-pinned. The threshold checker does not execute or import producer code. It derives the output itself and only then compares the complete primitive list, complete pair digest and grade endpoint histograms with the frozen producer certificate.

The reproducible portability command is:

    python -B verify_manifest.py --relocate

This verifies the full closed manifest, creates a fresh temporary directory by copying only the listed files, and executes all three scripts in both Python modes there. Its terminal result is emitted without modifying the frozen packet. A parent-level portability receipt may retain that output and the manifest digest externally, avoiding a self-referential receipt/manifest cycle.
