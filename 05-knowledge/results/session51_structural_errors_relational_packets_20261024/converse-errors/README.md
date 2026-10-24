# Converse-merged refinement-error extraction

Read `THEOREM.md` for definitions, proof, inherited dependencies, and limits.
The exact positive-error order-ten merged spectrum is
`6 t^14 + 23 t^18 + 2 t^36 + 54 t^40 + 7 t^72`: 92 nodes, total node-weighted
error3234. Exactly24 of the160 ordinary positive-error classes are
self-converse, all strong rigid two-five-cell classes with error40.

The general identities are `F=(1+F_s)A(z^2,t^2)` and `M=(A+F)/2`, where the
series count ordinary, self-converse, and converse-merged classes by size and
refinement error. This is not a blue/black wiring claim or a total-count result.

Reproduce from this directory:

    python3 check.py --output normal.json
    python3 -O check.py --output optimized.json
    cmp normal.json optimized.json
    python3 matrix_check.py --output matrix-normal.json
    python3 -O matrix_check.py --output matrix-optimized.json
    cmp matrix-normal.json matrix-optimized.json
    python3 series_check.py --output series-normal.json
    python3 -O series_check.py --output series-optimized.json
    cmp series-normal.json series-optimized.json
    python3 prime_check.py --output prime-normal.json
    python3 -O prime_check.py --output prime-optimized.json
    cmp prime-normal.json prime-optimized.json
    python3 port_check.py --output port-normal.json
    python3 -O port_check.py --output port-optimized.json
    cmp port-normal.json port-optimized.json
    python3 verify_manifest.py

The complete actual representative list, converse partner permutation, and
arc-map witnesses are in `normal.json`. Normal and optimized receipts are
retained. Source pins include the exhaustive order-ten classification and
regular-nine parent data; main and matrix checkers import no producer code.
`series_check.py` imports only this packet's own low-level graph checker.

Inheritance: the anti-involution existence/odd-power mechanism is already in
THM-647; the two-C5 converse coefficient formulas are in Session40's regular
block packet. The new use is the error-weighted merged extraction and transfer.

Native incidence consequence:24 mixed blue-incident nodes and68 black-only nodes.
Literal path checks give4200 blue presentation ports across the24, but no
line pairing or induced-edge count is inferred. The all-prime two-carousel
subfamily extraction and p=7 coefficient controls are in Section4.1.
