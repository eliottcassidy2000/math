# Independent deletion-comparison audit

Read REVIEW.md for the precise verdict and scope. The full tournament census is reproduced from the pinned original source; the twelve named response polynomials and 108 cards are independently recalculated.

From this directory, with the frozen producer in the neighboring session189-deletion-incoming directory:

    python3 audit.py > normal.stdout
    cp receipt.json receipt.normal.json
    python3 -O audit.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.normal.json receipt.json
    sha256sum -c SHA256SUMS

A relocated producer directory can be supplied as the sole argument to audit.py. The script uses Python's standard library, gcc, and g++; no packages or downloads are needed. It compiles two independent literal C++ engines, retains the full polynomial/card data, reproduces the original census through order nine, and checks both producer verification modes in a temporary directory.

The PUBLICATION_FILES.json allowlist excludes generated executables and redundant raw fetch copies. SOURCE_PROVENANCE.json records the independently checked upstream source blobs and URLs. All twelve orientation comparisons, with complete linear coefficient multisets, are retained in exclusions.json.
