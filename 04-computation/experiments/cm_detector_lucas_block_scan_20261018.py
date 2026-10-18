#!/usr/bin/env python3
"""A p-adic CM detector: second-order Lucas factorization of truncated hypergeometric sums (S43, HYP-9310).

For A(k) in four classical families and t_k = A(k) m^-k, test the Lucas-block congruence
    (B1')  sum_{k=p}^{2p-1} t_k == t_1 sum_{k=0}^{p-1} t_k   (mod p^2)
over all integers 0 < |m| <= M.  (It holds mod p for every m by the Lucas property; for A = C(2k,k)^3 it is
equivalent to the CM harmonic law of HYP-9310.)  Families: A1 = C(2k,k)^3, A2 = C(2k,k)^2 C(3k,k),
A3 = C(2k,k)^2 C(4k,2k), A4 = C(2k,k) C(3k,k) C(6k,3k).
PRE-REGISTERED (written before the run): (P1) for A1 the survivors are exactly Sun's eight CM values
{1, -8, 16, -64, 64, 256, -512, 4096}; (P2) for every family, every survivor is 'supersingular at density about 1/2':
sum_{k<p} t_k == 0 (mod p) at 35-65% of the primes 7 <= p <= 500 (a CM signature measured mod p, independent of the
mod-p^2 test), while non-survivors show density < 10%.
Usage: python cm_detector_lucas_block_scan_20261018.py [M]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb
import numpy as np


def primes_between(a, b):
    return [p for p in range(a, b + 1) if p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))]


FAM = {
    'A1=C(2k,k)^3': lambda k: comb(2 * k, k) ** 3,
    'A2=C(2k,k)^2C(3k,k)': lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k),
    'A3=C(2k,k)^2C(4k,2k)': lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k),
    'A4=C(2k,k)C(3k,k)C(6k,3k)': lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k),
}


def lucas_test(Af, ms, p):
    """boolean array: (B1') holds (or p | m) for each m in ms (numpy int64)."""
    pp = p * p
    A = [Af(k) % pp for k in range(2 * p)]
    r = np.mod(ms, pp)
    invtab = np.zeros(pp, dtype=np.int64)
    for x in range(pp):
        if x % p:
            invtab[x] = pow(x, -1, pp)
    w = invtab[r]
    S0 = np.zeros_like(ms)
    S1 = np.zeros_like(ms)
    pw = np.ones_like(ms)
    for k in range(2 * p):
        if k < p:
            S0 = (S0 + A[k] * pw) % pp
        else:
            S1 = (S1 + A[k] * pw) % pp
        pw = (pw * w) % pp
    ok = (S1 - (A[1] * w % pp) * S0) % pp == 0
    return ok | (r % p == 0)


def ss_density(Af, m, primes):
    z = 0
    tot = 0
    for p in primes:
        if m % p == 0:
            continue
        im = pow(m % p, -1, p)
        s = 0
        w = 1
        for k in range(p):
            s = (s + (Af(k) % p) * w) % p
            w = w * im % p
        tot += 1
        z += (s == 0)
    return z / tot if tot else float('nan')


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    ms_all = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    filt = primes_between(7, 43)
    conf = primes_between(47, 200)
    ssp = primes_between(7, 500)
    out = {}
    for name, Af in FAM.items():
        ms = ms_all.copy()
        for p in filt:
            ms = ms[lucas_test(Af, ms, p)]
        cand = ms
        for p in conf:
            cand = cand[lucas_test(Af, cand, p)]
        surv = [int(x) for x in cand]
        dens = {m: round(ss_density(Af, m, ssp), 3) for m in surv}
        # controls: a few non-survivors
        ctrl = {}
        for m in [2, 3, 5, 7, 11, 100, -100, 999]:
            if m not in surv:
                ctrl[m] = round(ss_density(Af, m, ssp), 3)
        print('%s: survivors (|m| <= %d, p <= 200): %s' % (name, M, surv))
        print('    supersingular density of survivors: %s' % dens)
        print('    controls (non-survivors): %s' % ctrl, flush=True)
        out[name] = dict(survivors=surv, ss_density=dens, controls=ctrl,
                         after_filter=int(len(ms)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_detector_lucas_block_scan_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
