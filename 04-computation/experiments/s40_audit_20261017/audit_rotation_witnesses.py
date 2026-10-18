"""S40 audit: explicit rotation-order witnesses against HYP-9287.

For a prime p = 3 mod 4 and a real alpha, the rotation order O_alpha lists x in {0, ..., p-1} by increasing frac(alpha x)
(ties impossible for alpha strictly between consecutive Farey fractions of order p - 1). This script scans alpha over
all midpoints of consecutive Farey fractions of order p - 1 (every distinct rotation order), reports the best raw
rotation order (no local search), the Farey interval it comes from, its excess recomputed by direct pair counting with
Legendre symbols from Euler's criterion, and the natural excess p h(-p). Witness orders are stored in the JSON.
Usage: python audit_rotation_witnesses.py P1 [P2 ...]
"""
import sys, os, math, json
from fractions import Fraction
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import legendre, h_forms


@njit(cache=True)
def leg_table(p):
    t = np.zeros(p, dtype=np.int64)
    for a in range(1, p):
        t[a] = legendre(a, p)
    return t


@njit(cache=True)
def excess(order, chi, p):
    s = 0
    for i in range(p):
        u = order[i]
        for j in range(i + 1, p):
            s += chi[(order[j] - u) % p]
    return s


def farey(n):
    a, b, c, d = 0, 1, 1, n
    out = [(0, 1)]
    while c <= n:
        out.append((c, d))
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
    return out


def rot_order(p, num, den):
    """Order of x by frac(x * num/den) computed exactly with integers: key = (x*num) mod den, scaled."""
    keys = [((x * num) % den, x) for x in range(p)]
    keys.sort()
    return [x for _, x in keys]


def main():
    res = []
    for p in [int(a) for a in sys.argv[1:]]:
        chi = leg_table(p)
        h = h_forms(-p)
        nat = int(excess(np.arange(p), chi, p))
        F = farey(p - 1)
        best = (-10 ** 9, None, None)
        for (a, b), (c, d) in zip(F[:-1], F[1:]):
            # exact midpoint alpha = (a/b + c/d)/2 = (a d + b c)/(2 b d)
            num, den = a * d + b * c, 2 * b * d
            o = np.array(rot_order(p, num, den), dtype=np.int64)
            e = int(excess(o, chi, p))
            if e > best[0]:
                best = (e, (a, b, c, d), o)
        e, (a, b, c, d), o = best
        assert sorted(o.tolist()) == list(range(p))
        row = dict(p=p, h=h, L=round(2 * math.pi * h / ((6 if p == 3 else 2) * math.sqrt(p)), 4), natural=nat, best_rotation_excess=e,
                   beaten=bool(e > nat), farey_interval='(%d/%d, %d/%d)' % (a, b, c, d),
                   alpha_midpoint='%d/%d' % ((a * d + b * c), 2 * b * d), MAS_lower_bound=(p * (p - 1) // 2 + e) // 2,
                   natural_forward=(p * (p - 1) // 2 + nat) // 2, order=o.tolist())
        res.append(row)
        print({k: v for k, v in row.items() if k != 'order'}, flush=True)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_rotation_witnesses.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
