# Independent order-eleven structural residual audit

`AUDIT.md` gives the independently derived exhaustive positive-error nonregular
classification. `audit.py` is standalone Python and imports no producer code.
The inherited regular-nine input is source-pinned; its atlas completeness is
inherited rather than inferred from reroot closure.

The exact reduction is

    E_s,11 = E_regular11 + 9166,
    E_11   = E_regular11 + 19480,
    U11 ≡10 (mod11).

The regular-eleven contribution is not computed by this packet. No unrestricted order-eleven
census, upper trace, or A000568(11) number is claimed.

Run normally and with optimization:

    python3 audit.py > normal.json
    python3 -O audit.py > optimized.json
    cmp normal.json optimized.json
    python3 hostiles.py > hostile-results.json
    python3 verify_manifest.py
    python3 replay.py

All semantic guards use explicit exceptions and survive `-O`. The12 hostile
runs challenge cell coalescence, missing matrix degrees, wrong isomorphism
quotients, erased full-orbit errors, and SCC coefficient indexing.
