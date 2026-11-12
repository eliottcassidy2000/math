#!/usr/bin/env python3
"""S93 lane "seeds": independent recount of the certificates certs_plus.txt / certs_minus.txt of b_census.py.

Each line: c v order colcode Omask M  (a pair (B, O) with offset +1 or -1 and slope M = c 2^v).
The closed tournaments T_t = (B => TT_t => closer, closer -> O), t = 1, 2, are built explicitly and their
Hamiltonian paths counted by a subset recursion that shares no code with the census (seedlib.hpaths: paths by end
vertex only).  Expected: H(T_t) = sign + 2^t M for t = 1, 2 (two values of t determine offset and slope).
Also: coverage of the odd N by the certified unit rays.
"""
import sys

import numpy as np

from seedlib import hpaths, closed_out, out_from_colcode


def main():
    tot = 0
    bad = 0
    summary = {}
    for sign, fn in ((1, "certs_plus.txt"), (-1, "certs_minus.txt")):
        rows = [tuple(int(x) for x in line.split()) for line in open(fn) if not line.startswith("#")]
        best = {}
        for c, v, n, code, O, M in rows:
            out = np.zeros(n, np.int64)
            out_from_colcode(n, code, out)
            if M != c << v:
                bad += 1
            for t in (1, 2):
                res = np.zeros(n + t + 1, np.int64)
                N = closed_out(n, out, O, t, res)
                if hpaths(N, res) != sign + (1 << t) * M:
                    bad += 1
                tot += 1
            best[c] = v
        summary[sign] = best
        nc = sum(1 for c in best if c <= 4097)
        print(f"sign {sign:+d}: {len(rows)} certificates (odd c <= 65535), {nc} of the 2049 odd c <= 4097; first odd c without a certificate: {min(c for c in range(1, 65537, 2) if c not in best)}; orders {sorted(set(r[2] for r in rows))}; "
              f"largest v {max(r[1] for r in rows)}; odd slopes (v = 0): {sum(1 for r in rows if r[1] == 0)}")
    print(f"closed tournaments recounted: {tot}; mismatches: {bad}")
    # coverage of odd N by the certified rays: N = 1 + c 2^s, s >= v_plus(c) + 1;  N = c 2^s - 1, s >= v_minus(c) + 1
    for k in range(4, 13):
        lo, hi = 1 << k, 1 << (k + 1)
        nodd = (hi - lo) // 2
        p = m = b = 0
        for N in range(lo + 1, hi, 2):
            x = N - 1
            s = (x & -x).bit_length() - 1
            c = x >> s
            onp = c in summary[1] and s >= summary[1][c] + 1
            y = N + 1
            s2 = (y & -y).bit_length() - 1
            c2 = y >> s2
            onm = c2 in summary[-1] and s2 >= summary[-1][c2] + 1
            p += onp
            m += onm
            b += (onp or onm)
        print(f"  odd N in [2^{k}, 2^{k + 1}): on a ray of offset +1: {p} ({p / nodd:.3f}); offset -1: {m} ({m / nodd:.3f}); either: {b} ({b / nodd:.3f}) of {nodd}")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
