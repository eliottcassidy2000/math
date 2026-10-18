#!/usr/bin/env python3
"""S41 audit A7: which logarithmic budgets does HYP-9283(B) actually certify, in the normalisation of HYP-9288 (B)?

In collatz_pole_certificate_margins_20261015.py the energy hypothesis enters only through
    S(a, b) = sum_(j < 400) sqrt(a + b j) q^j + sqrt(a + b) q^400 (400 - 399 q)/(1 - q)^2,   q = 2^-(sigma_0 - 1/2),
for budgets E_X <= X^2 (a + b j), X = N 2^j (j counted from the census size N); every other term of the certificate is
independent of (a, b), and the certificate passes iff S(a, b) < S*_N. The published margins (HYP-9283, (B)) therefore
lie on one level set, which is checked below; then the largest b is computed for budgets written as in HYP-9288 (B),
a + b log2 X = (a + b log2 N) + b j.
Usage: python a7_calibration.py
"""
import os, sys, json
sys.dont_write_bytecode = True
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    mp.mp.dps = 30
    Dfun = lambda s: 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3
    s0 = mp.findroot(Dfun, mp.mpc(0.99821, 240.248))
    sigma = s0.real
    q = mp.power(2, -(sigma - mp.mpf(1) / 2))
    J0 = 400

    def S(a, b):
        a = mp.mpf(a)
        b = mp.mpf(b)
        tot = mp.mpf(0)
        qj = mp.mpf(1)
        for j in range(J0):
            tot += mp.sqrt(a + b * j) * qj
            qj *= q
        tot += mp.sqrt(a + b) * q ** J0 * (J0 - (J0 - 1) * q) / (1 - q) ** 2
        return tot

    print('Mercator zero s_0 = %s, sigma_0 = %.8f, q = %.8f' % (mp.nstr(s0, 12), float(sigma), float(q)))
    out = dict(s0=str(s0), q=float(q))
    margins = {37: (1806.6, 1194.9, 292.86), 40: (14348.0, 9616.0, 2325.9)}
    for logN, (amax, bmax, cmax) in margins.items():
        lv = [S(amax, 0), S(2, bmax), S(4 * cmax, cmax)]
        Sstar = min(lv)
        print('N = 2^%d: S(a_max, 0) = %.4f, S(2, b_max) = %.4f, S(4c, c) = %.4f  (one level set: spread %.2e)'
              % (logN, lv[0], lv[1], lv[2], float((max(lv) - min(lv)) / min(lv))))

        def bmax_abs(a0):
            lo, hi = mp.mpf(0), mp.mpf(1e6)
            for _ in range(80):
                mid = (lo + hi) / 2
                if S(a0 + logN * mid, mid) < Sstar:
                    lo = mid
                else:
                    hi = mid
            return float(lo)
        b51 = bmax_abs(mp.mpf(51) / 7)
        b0 = bmax_abs(0)
        r1000 = float(S(51 / 7 + logN * 1000, 1000) / Sstar)
        out[logN] = dict(level_set=[float(x) for x in lv], b_max_abs_log_a51_7=b51, b_max_abs_log_a0=b0,
                         ratio_S_at_b1000_to_threshold=r1000)
        print('   budgets a + b log2 X: largest certified b with a = 51/7: %.2f; with a = 0: %.2f; at b = 1000 S/S* = %.2f'
              % (b51, b0, r1000))
    json.dump(out, open(os.path.join(HERE, 'a7_calibration.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
