# Reproducing the pair-entry bridge checkpoint

Status: **PROVED / VERIFIED** at the precise scope in PROOF.md. The exact pair-entry transfer, exposed-clone root switch, nonadjacent two-birth construction, and finite inverse criterion have independent review. No all-target scalar coverage theorem is claimed.

All Python code and graph data needed by these checks are included. Only the Python standard library is required. Run from this directory:

    sha256sum -c RELEASE_FILES.sha256
    python3 verify.py > replay.json
    cmp replay.json VALIDATION.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-optimized.json VALIDATION_OPTIMIZED.json
    python3 audit/independent.py > independent-replay.json
    cmp independent-replay.json audit/INDEPENDENT.json
    python3 -O audit/independent.py > independent-optimized.json
    cmp independent-optimized.json audit/INDEPENDENT_OPTIMIZED.json

The independent checker uses literal permutations and cuts, not the source's subset dynamic programs. Its report and exact source hashes are in audit/REPORT.md. PUBLICATION.sha256 additionally pins every delivered file other than that manifest itself.

The proof and report preserve their audited bytes and original research-path references. The final section of PROOF.md pins the structural proof dependencies by SHA256. Those dependencies supply the complete C1/C2 classifications and physical two-port compiler; the numerical verifier and independent checker here have no runtime imports from them.

The later response-cone corollary is deliberately not part of this checkpoint or its audit scope.
