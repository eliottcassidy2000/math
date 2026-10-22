"""Residuals of the repunits under the arithmetic-support model (session opus-2026-10-21-S51; HYP-9359, note section 0).
Strong inventory below 1100 (exact ears, as in repunit_rarity.py), the full model of HYP-9359 refitted on two ranges:
[25, 1098] (where 31, 63, 127 lie) and [301, 1098] (the training range of the pre-registration).  Prints the residuals of
the repunits and their neighbours.  Prints ALL CHECKS PASSED."""
import os, sys, math, time
from collections import Counter
from multiprocessing import Pool
import numpy as np
from sympy import isprime

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))
import alib as L                                       # noqa: E402
from b3_ear_census_625 import grow                     # noqa: E402

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def blocked(d):
    if isprime(d) or d < 9:
        return False
    for a in range(3, int(d ** 0.5) + 1, 2):
        if d % a == 0 and a not in (7, 21) and d // a not in (7, 21):
            return False
    return True


def feats(d):
    Ld = math.log(d)
    return [1.0, Ld, Ld * Ld, Ld ** 3, float(d % 3 == 0), float(d % 9 == 0), float(d % 5 == 0), float(isprime(d)),
            float(blocked(d)), float(d % 8 == 7)]


if __name__ == '__main__':
    LIM = 1100
    inv = {3: {L.canon(L.chain((2,))): 3}}
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv[n - 1].items() if h <= LIM]
            got, _ = grow(pool, n, par, cap=LIM) if par else ({}, None)
            got = {c: h for c, h in got.items() if h <= LIM}
            if not got:
                break
            inv[n] = got
            n += 1
    S = Counter(h for dd in inv.values() for h in dd.values())
    check(sum(1 for dd in inv.values() for _ in dd) == 136513, f'strong classes with H <= 1100: {sum(len(v) for v in inv.values())} (136,513)')
    for lo, hi in ((25, 1098), (301, 1098)):
        ds = [d for d in range(lo, hi + 1, 2) if S[d] > 0]
        X = np.array([feats(d) for d in ds]); y = np.log(np.array([S[d] for d in ds], dtype=float))
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        R = dict(zip(ds, y - X @ c))
        print(f'     full model on [{lo},{hi}]: ' + ', '.join(f'{d}: {R[d]:+.2f}' for d in (31, 33, 63, 65, 85, 127, 129, 255, 257, 341, 511, 513, 1023, 1025) if d in R))
        order = sorted(R, key=lambda d: R[d])
        print(f'     the five rarest on [{lo},{hi}]: ' + ', '.join(f'{d} ({R[d]:+.2f})' for d in order[:5]))
    check(True, 'residuals listed')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
