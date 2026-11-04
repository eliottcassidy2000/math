# Restore the exact audited response receipts

The normal and optimized receipts are byte-identical, each382,239 bytes. They are stored losslessly once as gzip/base64 chunks to keep the published transport files bounded. No coefficient or test result has been removed or rounded.

From this directory run:

    python3 reassemble.py
    sha256sum -c core.sha256
    python3 verify.py > replay.json
    cmp receipt.json replay.json
    python3 -O verify.py > replay.optimized.json
    cmp receipt.json replay.optimized.json

`REASSEMBLY.json` records every chunk's size/SHA-256 and both original output sizes/SHA-256. The script verifies the chunks before decoding and the restored bytes before writing. It refuses to replace an existing receipt with different bytes. The mathematical verifier is unchanged.

The independent response audit remains in `../response_audit/`. This packaging changes storage only. It does not extend the finite test universe or mathematical claims.
