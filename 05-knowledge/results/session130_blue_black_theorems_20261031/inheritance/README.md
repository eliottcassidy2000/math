# Session130 inheritance / equivariant OCF packet

Read **INHERITANCE_AND_THEORY.md** for the all-order length-refined
fixed-packing quotient, the n5 equivariant-conflict-graph hostile, and current
blue/black source corrections.

`normal.json` and `optimized.json` are byte-identical verified receipts from
`check_equivariant_ocf.py`; explicit exceptions remain active under `-O`.
They include the complete n3..6 reconstructed node/black-edge tables,
all reflection-fixed words through n7,128 deterministic n8 controls,
triangle-stack controls through order15, and literal n5 paths/cycle bijection.
No inherited producer or data import is used by this checker.

`SOURCE_PINS.json` records exact Git and SHA256 source pins. The fresh
THM-4854 comparison is treated as PREAUDIT, not silently promoted.

Reproduce from this directory:

    python3 check_equivariant_ocf.py > normal.json
    python3 -O check_equivariant_ocf.py > optimized.json
    cmp normal.json optimized.json

This is a portable research handoff, not a publication or theorem-ID reservation.
