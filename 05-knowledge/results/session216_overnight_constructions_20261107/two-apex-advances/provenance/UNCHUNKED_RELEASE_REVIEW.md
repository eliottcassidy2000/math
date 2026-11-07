# Portable release review

Scope: approved Session216 post-freeze inverse advances, independently audited. Mathematical claims and executable research sources are unchanged. No ambient coverage extension, all-odd theorem, or unrestricted minimum-order claim is made.

## Contents and provenance

Every payload file named in provenance/ADVANCE_FILES.txt is present. Except ADVANCES.md, each is byte-identical to the source verified by provenance/ADVANCE_SHA256SUMS. The original audited ADVANCES.md is preserved as provenance/ADVANCES.original.md; provenance/ADVANCES.diff records all wording changes. The only changes clarify map-counter units, binary endianness, and executable replay commands. Existing baseline helpers required by the programs are already in the approved list. audit/ contains the independent audit and essential reproducibility sources/receipts. There are no compiled binaries or uncompressed census/certificate files in this release.

Exact SHA256 pins:

- Original ADVANCES.md:54c5ea097e542b55e7738b7c3cbf3932f5fe4966f1d7bc809a24e01d7ce2aac5
- Revised ADVANCES.md:5afc94df68ec8e8e4199d29557ba4a3d42a35dc2bc24c002f40ade031b53810a
- coverage2327673.bin.gz:b0786da688454f1cad753b13ac1ee055c0210a06c67c700fa699e04679537394
- Original ADVANCE_SHA256SUMS:de1b23418d6b671b3972077eb35adb8f3d7b4ce59ff318890b76ae9032dec370

The certificate is4,002,651 compressed bytes and decompresses to9,310,696 bytes. It stores2,327,674 explicitly little-endian uint32 records. Valid records encode length in the top4 bits and increasing-order letters in successive low2-bit pairs;0xffffffff denotes absence. The census C++ writer and independent C++ reader use native-endian words. The replay script explicitly rejects a big-endian host rather than silently misreading the certificate. On such a host explicit byte conversion would be required.

## Verification and commands

From the release directory run:

    sha256sum -c SHA256SUMS
    sh replay.sh

Required tools: POSIX shell, g++ with C++17, Python3 standard library, gzip, cmp, mktemp. No network, external package, or user credentials are needed. The replay writes its large intermediate data and compiled programs only into an automatically cleaned temporary directory. It checks the exact full-depth14 census receipt, byte-exact equality of the full initial-range certificate, every saved witness with the independent36-coordinate transfer,297 sampled certificate witnesses by independent subset DP, the boundary repair by independent subset DP, the independent mode census and15 mode subset-DP checks, and all bounded map-probe/mode receipts. Timings are excluded only when comparing receipt fields containing elapsed seconds.

The complete certificate has1,163,835 witnesses. All were checked by the independent36-coordinate transfer; only297 sampled certificate witnesses were checked by subset DP. The boundary word and15 formula values were additionally checked by subset DP. No exhaustive per-witness subset-DP claim is intended.

The historical audit_large.out remains a receipt for1,200 sampled witnesses from the earlier cap1,000,000 run, not from this certificate. Its raw input is intentionally omitted. To replay that historical sample separately, compile census.cpp, run with arguments14 1000000 and an absolute output path, compile audit.cpp specifically to /tmp/session216-inverse-audit as required by audit_large.py, then pass that absolute census path to audit_large.py. The main replay does not depend on this historical receipt or its hard-coded helper path.

The included audit report establishes finite exact coverage through2,327,673 except7,21 and verifies the first extra gap2,327,675 at chain length<=14. The15-letter repair proves minimum chain length15 in this precise two-apex class only. All pure-mode formulas use the inherited homogeneous-spectrum proof, with m=0 checked separately. Map probes are finite restricted-action tests only.

SHA256SUMS is the authoritative release manifest. RELEASE_FILES.txt lists all manifested payload, audit, documentation, and provenance files; the SHA256SUMS file itself is excluded to avoid circular hashing. The preserved source manifest is historical and should not be applied directly to the revised release root.

Release-local end-to-end replay completed successfully. Its exact transcript is provenance/REPLAY.out. All independent audit source hashes also pass.
