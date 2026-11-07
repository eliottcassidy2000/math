# Session216 inverse advances: chunked transport release

This is the same independently audited finite two-apex certificate, packaged into31 lossless binary parts for smaller individual uploads. Mathematical programs and recorded results are unchanged. No general Chan–Pak completeness claim is made.

The original4,002,651-byte gzip is NOT stored in this variant. It is split in byte order into30 parts of131,072 bytes (128KiB) and one70,491-byte final part, named certificate/coverage2327673.bin.gz.part000 through part030. Each full part requires174,764 base64 characters if transported through a base64 blob API. Upload parts separately, not in one large inline tree.

The exact original gzip SHA256 is:

b0786da688454f1cad753b13ac1ee055c0210a06c67c700fa699e04679537394

certificate/MANIFEST.json records ordered filenames, exact byte lengths, individual SHA256 hashes and the reconstructed original hash. SHA256SUMS covers all release files including every part and the certificate manifest.

## Replay

Run from this directory:

    sha256sum -c SHA256SUMS
    sh replay.sh

The replay reconstructs the gzip into temporary space, verifies each part and the complete original hash, decompresses it, and checks the full census plus the independent certificate/mode/probe audit. Temporary large files and binaries are removed on exit. POSIX shell, g++, Python3, gzip, cmp and mktemp are required. The native C++ binary routines require a little-endian host, which replay.sh checks explicitly.

For certificate-only reconstruction and verification:

    tmp=$(mktemp -d)
    python3 reconstruct_certificate.py "$tmp/coverage2327673.bin.gz"
    gzip -cd "$tmp/coverage2327673.bin.gz" > "$tmp/coverage.bin"
    g++ -O3 -std=c++17 audit/check.cpp -o "$tmp/check"
    "$tmp/check" "$tmp/coverage.bin"
    rm -rf "$tmp"

The original packed certificate is9,310,696 decompressed bytes containing2,327,674 little-endian uint32 records. The complete set of1,163,835 stored witnesses was checked with the independent36-coordinate positive transfer. Independent subset DP checks297 sampled certificate witnesses, plus the boundary repair separately. This is not an exhaustive per-witness subset-DP claim.

## Provenance

provenance/CHUNKED_PROVENANCE.json pins the unchunked source release and the original gzip. provenance/PACKAGING.diff records the ADVANCES.md and replay.sh packaging-only changes. Historical source manifests and provenance/UNCHUNKED_RELEASE_REVIEW.md describe the unchunked predecessor; their direct gzip paths are superseded by this README. The current release manifest is SHA256SUMS. No raw9MB or16MB data files, compiled binaries, or redundant original gzip are included.
