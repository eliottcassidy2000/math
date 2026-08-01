"""Diagnose HOW the final row fails at R=128 (and how it SUCCEEDS at R=64).

The final row admits NO freedom: B_{d,k} = x^{d-k}(1-x)^k has lowest term
x^{d-k}, so reading res[0], res[1], ... determines delta_d, delta_{d-1}, ...
uniquely.  Either the forced reduction runs to completion or it does not.
Failure modes: (F0) residual DEGREE exceeds d, state discarded before we look;
(F1) |want| > cap; (F2) parity wrong; (F3) leftover nonzero remainder.
"""
from math import comb
from itertools import product
import amm12592_gamma35_beam_deathstar as beam

def final_diagnose(R, D0=0, bw=250, ctrl=2, span=2):
    d = [(3 * (R + i)) // 5 + D0 for i in range(R)]
    states = [([], beam.qpow(R - 1))]
    opts = [None] + list(range(-span, span + 1))
    for i in range(R - 1):
        nxt = []
        for acc, sig in states:
            for tg in product(opts, repeat=ctrl):
                if tg[0] not in (1, -1): continue
                r = beam.step(sig, d[i], tg)
                if r is None: continue
                de, ns = r
                if not ns or abs(ns[0]) != 1: continue
                nxt.append((acc + [de], ns))
        if not nxt: return f"died at row {i}", None
        nxt.sort(key=lambda s: (len(s[1]), sum(abs(c) for c in s[1][:6])))
        seen, uniq = set(), []
        for a, s in nxt:
            k = tuple(s[:10])
            if k in seen: continue
            seen.add(k); uniq.append((a, s)); 
            if len(uniq) >= bw: break
        states = uniq
    di = d[R - 1]
    modes = {"F0_degree": 0, "F1_capacity": 0, "F2_parity": 0, "F3_remainder": 0, "OK": 0}
    worst_over = 0.0
    degs = []
    for acc, sig in states:
        degs.append(len(sig) - 1)
        if len(sig) - 1 > di:
            modes["F0_degree"] += 1; continue
        res = list(sig) + [0] * (di + 2)
        failed = None
        for k in range(di, -1, -1):
            cap = comb(di, k); want = res[di - k]
            if abs(want) > cap:
                failed = "F1_capacity"; worst_over = max(worst_over, abs(want) / cap); break
            if (want - cap) % 2:
                failed = "F2_parity"; break
            for t, c in enumerate(beam.basis_poly(di, k)): res[t] -= want * c
        if failed: modes[failed] += 1; continue
        modes["OK" if not beam.trim(res) else "F3_remainder"] += 1
    return modes, (di, min(degs), max(degs), worst_over)

for (R, D0) in [(64, 0), (128, 0), (128, 2)]:
    modes, info = final_diagnose(R, D0=D0)
    if info is None:
        print(f"R={R} D0={D0}: {modes}"); continue
    di, dmin, dmax, wo = info
    print(f"R={R:3d} D0={D0}:  final-row degree budget d={di},  residual degrees in "
          f"[{dmin},{dmax}]")
    print(f"          modes over the {sum(modes.values())} surviving states: {modes}"
          f"{'  worst |want|/cap = %.2f' % wo if wo else ''}")
