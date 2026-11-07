# Chunked transport release review

Verdict: packaging-only variant of the independently audited Session216 advance release. Mathematical source programs, numerical receipts, independent checker programs, and certificate contents are unchanged. The complete original gzip is represented by31 ordered binary parts; it is intentionally absent as a standalone file.

Original release manifest SHA256:
7e0fc606963ebc2f7a7305b5f24d14fa5cc33711b1127243be4715b06095413c

Original gzip SHA256:
b0786da688454f1cad753b13ac1ee055c0210a06c67c700fa699e04679537394

Parts000..029 each contain131,072 bytes; part030 contains70,491 bytes. Reconstruction totals4,002,651 bytes. The ordered certificate/MANIFEST.json records hashes for each part and for the original gzip. Its own hash and all other release-file hashes are recorded by the current SHA256SUMS.

Verification used independent shell concatenation, SHA256 comparison, and byte comparison against the original unchunked gzip. The reconstruction was decompressed and the independent certificate checker rerun: all1,163,835 saved witnesses agree with the compressed36-coordinate transfer;297 sampled certificate witnesses agree with independent subset DP; the boundary witness agrees with both. The exact checker receipt is provenance/CHUNKED_CERTIFICATE_CHECK.out.

The updated replay.sh reconstructs only in temporary space, verifies hashes before decompression, compares the complete certificate against the full census, and runs the mode and probe controls. README.md gives both the full replay and certificate-only commands. Native binary C++ readers/writers require little-endian hosts, checked by the full replay script.

Historical unchunked provenance is retained for traceability. Its references to a directly stored gzip do not describe this variant. provenance/PACKAGING.diff records changes to existing text/script files; provenance/PACKAGING_CHANGES.json enumerates added, removed and changed files. The unchunked review is preserved at provenance/UNCHUNKED_RELEASE_REVIEW.md. No general completeness or exhaustive per-witness subset-DP claim has been introduced.

Current integrity commands:

    sha256sum -c SHA256SUMS
    sh replay.sh

The release contains no compiled binaries, uncompressed census data, or redundant original gzip. Nothing was edited in the staging repository or published by this packaging task.

The chunked release full replay completed successfully; see provenance/CHUNKED_REPLAY.out.
