# Session 127 volume lane

Read `THEORY.md` for the proved bridge and obstructions. Reproduce the finite exact algebra with:

```sh
python3 verify_volume_bridge.py
```

Outputs: `verification.json` and a one-line PASS summary. Requires Python 3 and SymPy (version recorded in the output). `sources.json` lists the primary references; `SHA256SUMS` pins this packet.

Main result: degree-k full norm lattices use `V_k`, while binary Thue slices have rank 2, Gram determinant `k S_k-1`, and a polynomial-dependent norm-body area. Norm collapse is not Euclidean smallness; the exact tribonacci input `(103,56)` already has norm `-1` and length `>180`.
