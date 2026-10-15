#!/usr/bin/env python3
"""Audit C: test HYP-9276 (M_q = -b(q^2) - 2q - 1) with this audit's own M_q evaluators and own b(p).

M_q evaluators (audit_common.py):
  FFT   : multiplicative-group FFT (discrete logs), rounding gap recorded         -> every prime 5 <= q <= 100,000
                                                                                     + random sample of 1,000 primes in (1e5, 3e5]
  EXACT : direct int64 correlation/convolution (no FFT), O(q^2)                    -> every prime 5 <= q <= 3,000
                                                                                     + random sample of 300 primes in (3e3, 3e4]
                                                                                     + random sample of 30 primes in (3e4, 6e4]
b(p): independent eta/theta construction (audit_newform.py -> audit_newform_bp.json), cross-checked against PARI.
Also: comparison with the audited session's k6p6_fft.json (M_q, W_q at the 3,243 primes to 30,000); zeros of M_q;
min v_2(M_q) per class mod 8; moments of x_q = (|b(q)|^2 - q^2)/q^2 by class mod 4 (q <= 30,000 and q <= 100,000).

Usage: python3 -u audit_identity.py > audit_identity.out
"""
import sys
import json
import time
import numpy as np
import sympy
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014')
from audit_common import M_mult_fft, M_exact_direct, chi_prime, W_from_M

D = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/'
bp = {int(k): tuple(v) for k, v in json.load(open(D + 'audit_C_20261014/audit_newform_bp.json')).items()}


def predicted_M(q):
    a, b = bp[q]                     # b(q) = a + b i
    chi4 = 1 if q % 4 == 1 else -1
    re_bq2 = a * a - b * b - chi4 * q * q
    im_bq2 = 2 * a * b
    assert im_bq2 == 0
    return -re_bq2 - 2 * q - 1


def main():
    t0 = time.time()
    rng = np.random.default_rng(14102026)
    allp = list(sympy.primerange(5, 300001))
    fft_set = [p for p in allp if p <= 100000]
    big = [p for p in allp if p > 100000]
    fft_set += sorted(rng.choice(big, size=1000, replace=False).tolist())
    ex_set = [p for p in allp if p <= 3000]
    mid = [p for p in allp if 3000 < p <= 30000]
    ex_set += sorted(rng.choice(mid, size=300, replace=False).tolist())
    hi = [p for p in allp if 30000 < p <= 60000]
    ex_set += sorted(rng.choice(hi, size=30, replace=False).tolist())

    Mfft = {}
    maxdev = 0.0
    for q in fft_set:
        M, dev, _, _ = M_mult_fft(q)
        Mfft[q] = M
        maxdev = max(maxdev, dev)
    print(f"FFT: {len(Mfft)} primes (all 5..100000 and 1000 random in (1e5, 3e5]); worst rounding gap {maxdev:.2e} "
          f"[{time.time() - t0:.1f}s]", flush=True)
    Mex = {}
    for q in ex_set:
        Mex[q] = M_exact_direct(q)[0]
    print(f"EXACT: {len(Mex)} primes (all 5..3000, 300 random in (3e3, 3e4], 30 random in (3e4, 6e4]) "
          f"[{time.time() - t0:.1f}s]", flush=True)
    both = [q for q in Mex if q in Mfft]
    dis = [q for q in both if Mex[q] != Mfft[q]]
    print(f"EXACT vs FFT agree at {len(both) - len(dis)} of {len(both)} common primes; disagreements {dis}")
    # HYP-9276 test
    fails_fft = [q for q in Mfft if Mfft[q] != predicted_M(q)]
    fails_ex = [q for q in Mex if Mex[q] != predicted_M(q)]
    print(f"HYP-9276 M_q = -b(q^2) - 2q - 1: FFT set {len(Mfft)} primes, failures {fails_fft[:10]}; "
          f"EXACT set {len(Mex)} primes, failures {fails_ex[:10]}")
    print(f"   largest prime tested: {max(Mfft)} (FFT), {max(Mex)} (EXACT)")
    # comparison with the audited session's stored values
    sess = {int(k): v for k, v in json.load(open(D + 'paley_k6p6_20261014/k6p6_fft.json')).items()}
    sm = [q for q in sess if Mfft.get(q) != sess[q][0]]
    sw = [q for q in sess if W_from_M(q, 1 if q % 4 == 1 else -1, Mfft[q]) != sess[q][1]]
    print(f"session k6p6_fft.json: {len(sess)} primes; M_q mismatches {sm[:10]}; W_q mismatches {sw[:10]}")
    # zeros, 2-adic, moments
    for lim in (30000, 100000):
        Q = [q for q in Mfft if q <= lim]
        zeros = [q for q in Q if Mfft[q] == 0]
        v2 = {}
        v2_1mod4_not1 = []
        for q in Q:
            m = Mfft[q]
            if m:
                v = (abs(m) & -abs(m)).bit_length() - 1
                v2[q % 8] = min(v2.get(q % 8, 99), v)
                if q % 4 == 1 and v != 1:
                    v2_1mod4_not1.append(q)
        absb_eq = [q for q in Q if q % 4 == 3 and bp[q][1] ** 2 == (q + 1) ** 2]
        print(f"q <= {lim}: {len(Q)} primes ({sum(q % 4 == 3 for q in Q)} = 3 mod 4, {sum(q % 4 == 1 for q in Q)} = 1 mod 4)")
        print(f"   M_q = 0 at {zeros}  (all = 7 mod 8: {all(q % 8 == 7 for q in zeros)}); "
              f"Paley primes with |b(q)| = q+1: {absb_eq == zeros and 'same list' or absb_eq}")
        print(f"   min v_2(M_q) by class mod 8: {dict(sorted(v2.items()))}; q = 1 mod 4 with v_2(M_q) != 1: {v2_1mod4_not1[:10]}")
        ri = [1, 0, 1, 1, 3, 6, 15, 36, 91, 232, 603, 1585, 4213, 11298, 30537, 83097, 227475]
        for cls in (3, 1):
            xs = np.array([(bp[q][0] ** 2 + bp[q][1] ** 2 - q * q) / (q * q) for q in Q if q % 4 == cls])
            n = len(xs)
            mom = [float(np.mean(xs ** k)) for k in range(0, 9)]
            se = [np.sqrt((ri[2 * k] - ri[k] ** 2) / n) for k in range(0, 9)]
            z = [(mom[k] - ri[k]) / se[k] if se[k] > 0 else 0.0 for k in range(9)]
            print(f"   class {cls} mod 4, N = {n}: moments k=0..8 " + ", ".join(f"{m:.3f}" for m in mom)
                  + "  | z-scores k=1..8 " + ", ".join(f"{v:+.2f}" for v in z[1:]))
    json.dump({'fft': {str(k): v for k, v in Mfft.items()}, 'exact': {str(k): v for k, v in Mex.items()}},
              open(D + 'audit_C_20261014/audit_identity_M.json', 'w'))
    print(f"done [{time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
