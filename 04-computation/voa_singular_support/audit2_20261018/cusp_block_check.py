#!/usr/bin/env python3
"""Locate the (weight, Li) block where HS J_infty(y^2=x^3) first exceeds the complete-intersection count, and
recompute that block's rank at a second prime and exactly over Z (fmpz_mat)."""
import sys
import flint
from arc_audit import Arc, cusp_ci
W = int(sys.argv[1]) if len(sys.argv) > 1 else 34
P1, P2 = 2147483647, 2147483629
A1 = Arc([{(0, 2): 1, (3, 0): -1}], W, P1)
A2 = Arc([{(0, 2): 1, (3, 0): -1}], W, P2)
# complete-intersection count per block: bigraded Euler characteristic.  Use the bigraded CI formula:
# HS_CI(t,s) = prod_j 1/((1-t^{2+j}s^j)(1-t^{3+j}s^j)) * prod_n (1 - t^{6+n} s^n)
def bigraded_ci(W):
    # dict (w,l) -> coeff
    ser = {(0, 0): 1}
    def mul_factor(ser, w, l, sign):  # multiply by (1 - t^w s^l)^{sign} ; sign=-1 -> 1/(1-..), sign=+1 -> (1-..)
        if sign == 1:
            out = dict(ser)
            for (a, b), c in ser.items():
                if a + w <= W:
                    out[(a + w, b + l)] = out.get((a + w, b + l), 0) - c
            return out
        out = dict(ser)
        # geometric series: process in increasing weight
        keys = sorted(out.keys())
        res = {}
        for (a, b) in sorted(set(k for k in out), key=lambda k: k):
            pass
        res = dict(ser)
        for a in range(0, W + 1):
            for (aa, bb) in [k for k in list(res.keys()) if k[0] == a]:
                c = res[(aa, bb)]
                if aa + w <= W and c:
                    res[(aa + w, bb + l)] = res.get((aa + w, bb + l), 0) + c
        return res
    for j in range(0, W + 1):
        if 2 + j <= W: ser = mul_factor(ser, 2 + j, j, -1)
        if 3 + j <= W: ser = mul_factor(ser, 3 + j, j, -1)
    for n in range(0, W + 1):
        if 6 + n <= W: ser = mul_factor(ser, 6 + n, n, 1)
    return ser
ci = bigraded_ci(W)
for L in range(0, W + 1):
    nc, r1 = A1.block_rank(W, L)
    if nc == 0:
        continue
    _, r2 = A2.block_rank(W, L)
    hs1 = nc - r1
    c = ci.get((W, L), 0)
    if hs1 != c or r1 != r2:
        print(f'W={W} L={L}: #mon={nc} rank(P1)={r1} rank(P2)={r2} HS={hs1} CI={c}', flush=True)
        # exact rank over Z
        P = P1
        cols = A1.monos(W, L); ci_ = {m: i for i, m in enumerate(cols)}
        # rebuild integer rows (F_n has small integer coefficients; recompute with a huge modulus-free approach)
        Aint = Arc([{(0, 2): 1, (3, 0): -1}], W, 2**61 - 1)  # Mersenne prime > coefficients: coefficients are exact integers mod it
        rows = []
        for (wf, F) in Aint.Fs:
            for n in range(0, L + 1):
                Wm, Lm = W - wf - n, L - n
                if Wm < 0: continue
                for m in Aint.monos(Wm, Lm):
                    row = [0] * nc
                    for mo, co in F[n].items():
                        k = ci_[tuple(sorted(mo + m))]
                        cc = co if co < 2**60 else co - (2**61 - 1)
                        row[k] += cc
                    rows.append(row)
        M = flint.fmpz_mat(rows)
        print('   exact integer rank:', M.rank(), ' -> HS over Q =', nc - M.rank(), flush=True)
print('total HS(P1) at W', A1.hs(W), ' CI total', sum(v for (a, b), v in ci.items() if a == W))
