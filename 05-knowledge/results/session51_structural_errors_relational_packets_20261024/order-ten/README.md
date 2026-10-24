# Order-ten structural error certificate

The exact strong-order-ten refinement error is3992, with disjoint contributions
3360 from5+5 stable cells,56 from7+1+1+1, and576 from3+3+3+1.
Together with the inherited Session50 strong-order-nine coefficient716, the SCC
transfer gives total error5508 and the reduction a_10=(U_10-5508)/10.
The upper trace U_10 is not computed by this packet.

Run from any directory, with Python3 and no external packages or inputs:

    python3 /path/to/check.py > normal.json
    python3 -O /path/to/check.py > optimized.json
    cmp normal.json optimized.json
    python3 /path/to/verify_manifest.py

`THEOREM.md` proves completeness, recognition of actual stable cells, exact
isomorphism quotients, all strongness conditions, the crossing counts, and a
prime-cell family generalization. `check.py` is independently written and uses
literal matrix and block presentation orbits, plus unrestricted degree-pruned
arc checks for automorphisms. All tests are explicit exceptions, never Python
assert statements. No full order-ten tournament atlas/count is an input.

`sources/` preserves the Session50 predecessor theorem and receipt supporting
E_s,9=716 and the SCC transfer. These files are not read by `check.py`.
Its only inherited arithmetic input is716; the strong order-ten coefficient
and regular-seven coefficient are rederived from small-sector enumerations.
`SOURCE_PINS.json` records exact upstream paths and hashes. `MANIFEST.json`
records all deliverable bytes except itself and verification receipts.

`distribution_check.py` additionally rebuilds the regular-nine error histogram
using only85 inherited almost-regular eight-vertex parents and eight-vertex
canonicalizations. Together with the literal order-ten family results this
certifies160 errorful order-ten classes, with distribution
12 of error14,46 of error18,4 of error36,84 of error40,14 of error72.
Its own input contract and the bivariate SCC law are explicit in THEOREM.md.
Run it normally and with `-O`; compare the two distribution receipts.
