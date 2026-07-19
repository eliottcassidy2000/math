#!/usr/bin/env python3
"""runs_parity_D3P_kps_S128c83c.py -- kind-pasteur-2026-07-18-S128 (cont.83)

A note on death-star-S58c's residue (HYP-7737).  They reduce the maximiser bound to

    #runs <= 2D/3        (D = max component),

verified for primitive d with components <= 30.  This script checks the restatement:

    mirror symmetry (death-star-S58b) pairs the runs, so #runs = 2P with P = #mirror
    pairs -- PROVIDED no run is self-mirror.  If #runs is always even, their inequality
    is exactly

    D >= 3P .

That is the same shape as the ceiling itself (2/21 = 2/(7*3)), and it is the form in
which my independent mirror-charged route (THM-1211) landed: rho >= 3P with
rho = max(d_exit(r), d_exit(r*)) <= D.  So rho >= 3P => D >= 3P: my version is
strictly stronger, theirs is the one that needs proving.

Also records the box/simplex ratio: death-star's run bound uses B's BOUNDING BOX
[1/7,2/7]x[3/7,4/7]x[5/7,6/7] (side 1/7).  B is the inscribed simplex, exactly 1/6 of it.
"""
import sys
from fractions import Fraction as F
from math import gcd

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30


def bad_intervals(DD):
    """death-star-S58b's exact affine-cell decomposition (verbatim logic)."""
    bps = {F(0), F(1)}
    for d in DD:
        for k in range(d + 1):
            bps.add(F(k, d))
    for i in range(3):
        for j in range(i + 1, 3):
            dd = abs(DD[i] - DD[j])
            if dd > 0:
                for k in range(dd + 1):
                    bps.add(F(k, dd))
    bps = sorted(bps)
    ivs = []
    for idx in range(len(bps) - 1):
        lo, hi = bps[idx], bps[idx + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        gmid = [(-DD[i] * mid) % 1 for i in range(3)]
        a = [gmid[i] + DD[i] * mid for i in range(3)]
        order = sorted(range(3), key=lambda i: gmid[i])
        s0, s1, s2 = order
        cons = [("ge", a[s0] - F(2, 7), DD[s0]), ("le", a[s2] - F(5, 7), DD[s2])]
        feas = True
        for (x, y) in [(s1, s0), (s2, s1)]:
            c = a[x] - a[y]
            dc = DD[x] - DD[y]
            if dc > 0:
                cons.append(("ge", c - F(2, 7), dc))
            elif dc < 0:
                cons.append(("le", c - F(2, 7), dc))
            elif not (c <= F(2, 7)):
                feas = False
        if not feas:
            continue
        ulo, uhi = lo, hi
        for typ, const, dd in cons:
            b = const / dd
            if typ == "ge":
                ulo = max(ulo, b)
            else:
                uhi = min(uhi, b)
        if uhi > ulo:
            ivs.append((ulo, uhi))
    ivs.sort()
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


odd = 0
odd_ex = None
selfmirror = 0
viol = 0
worstP = {}
tested = 0
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            ivs = bad_intervals(DD)
            if not ivs:
                continue
            tested += 1
            D = max(DD)
            k = len(ivs)
            if k % 2:
                odd += 1
                if odd_ex is None:
                    odd_ex = (tuple(DD), k)
            S = set(ivs)
            for lo, hi in ivs:
                if (1 - hi, 1 - lo) == (lo, hi):
                    selfmirror += 1
            P = k // 2
            if P and D < 3 * P:
                viol += 1
            rec = worstP.setdefault(P, [10 ** 9, None])
            if D < rec[0]:
                rec[0] = D
                rec[1] = tuple(DD)

print("primitive directions with a nonempty bad set, components 1..%d : %d" % (N, tested))
print()
print("  directions with an ODD number of runs : %d %s"
      % (odd, "" if odd_ex is None else "(first: d=%s, #runs=%d)" % odd_ex))
print("  self-mirror runs found                : %d" % selfmirror)
print("  -> #runs = 2P is %s, so death-star's #runs <= 2D/3 is exactly  D >= 3P"
      % ("VALID" if odd == 0 else "NOT valid in general"))
print()
print("  violations of D >= 3P : %d" % viol)
print()
print("  P  |  min D observed  |  witness        |  3P")
for P in sorted(worstP):
    if P == 0:
        continue
    lo, wit = worstP[P]
    print("  %-3d|  %-16d|  %-16s|  %d" % (P, lo, wit, 3 * P))
print()
print("  BOX vs SIMPLEX: death-star's run bound uses B's bounding box")
print("    [1/7,2/7] x [3/7,4/7] x [5/7,6/7], side 1/7, volume (1/7)^3 per ordering.")
print("    B itself is the inscribed simplex (THM-1211): volume (1/7)^3/6 per ordering,")
print("    so |B| = 6*(1/7)^3/6 = 1/343 and the box overstates B by exactly 6x.")
print("    The run bound 1/(7D) is TIGHT at (1,2,3), so the simplex cannot improve it --")
print("    the place the simplex may bite is the residual COUNTING inequality.")
