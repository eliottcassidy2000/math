"""Covering density D_k(t) of HYP-9470.

Usage: python3 density.py K mats1.txt ... matsK.txt
matsN.txt = sorted distinct transfer matrices of marked gadgets of order exactly N (9 integers per line, row-major,
SEAM convention: v(X (.) G, z) = M v(X, a), v = (H, S, R)), as printed by gadget_matrices.c.
Prints, for k <= K: the number of distinct S-row-zero matrices, the number of distinct nonsingular 2x2 restrictions
[[a, b], [c, d]] acting on (H, R) on the plane S = 0, and D_k(t) = sum over those restrictions of (a + b t)^(-2).
"""
import sys
K = int(sys.argv[1]); files = sys.argv[2:2 + K]
order = {}
for k, f in enumerate(files, start=1):
    for line in open(f):
        m = tuple(map(int, line.split()))
        order.setdefault(m, k)
for k in range(1, K + 1):
    srow0 = [m for m, o in order.items() if o <= k and m[3] == 0 and m[4] == 0 and m[5] == 0]
    restr = {}
    for m in srow0:
        A = (m[0], m[2], m[6], m[8])
        if A[0] * A[3] - A[1] * A[2] != 0:
            restr.setdefault(A, k)
    dets = sorted(set(A[0] * A[3] - A[1] * A[2] for A in restr))
    odd_det = [d for d in dets if d % 2]
    Ds = {t: sum((A[0] + A[1] * t) ** -2 for A in restr) for t in (2.0, 2.5, 3.0, 4.0)}
    print(f"k<={k}: S-row-zero matrices {len(srow0)}, nonsingular 2x2 restrictions {len(restr)}, "
          f"odd determinants {odd_det[:5]}, " + ", ".join(f"D({t})={v:.4f}" for t, v in Ds.items()))
