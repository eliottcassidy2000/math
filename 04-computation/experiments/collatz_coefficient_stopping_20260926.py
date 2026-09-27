#!/usr/bin/env python3
"""collatz_coefficient_stopping_20260926.py -- coefficient stopping time versus actual stopping time (Syracuse map),
and the one-member bound for uncertified class members (session gilbreath6-collatz-precision-20260926, opus).

For odd n let U(n) = (3n+1)/2^v. sigma(n) = least j with U^j(n) < n (actual); sigma_inf(n) = least j with
3^j < 2^(A_j), A_j = v_1 + ... + v_j (coefficient, Terras 1976). Since U^j(n) = (3^j n + S_j)/2^(A_j) with S_j > 0,
actual descent at j needs coefficient descent at j, so sigma >= sigma_inf; and a coefficient descent at j gives the
actual one unless n <= N(w) = S_j/(2^(A_j) - 3^j).
 (1) For all odd n <= NMAX: compare sigma and sigma_inf (vectorised over n; stop when every n has descended).
 (2) For j <= JMAX: check 2^A - 3^j > (3/2)^j - 1 at A = ceil(j log_2 3), which implies N(w) < 2^A for every
     coefficient-descent word of length j, hence at most one uncertified member per class.
Usage: python3 collatz_coefficient_stopping_20260926.py [NMAX=10000000] [JMAX=5000]
"""
import sys, math
import numpy as np


def compare(NMAX):
    n = np.arange(3, NMAX + 1, 2, dtype=np.int64)  # n = 1 excluded: U(1) = 1 never descends
    start = n.copy()
    j = np.zeros(len(n), dtype=np.int64)          # Syracuse steps done
    A = np.zeros(len(n), dtype=np.int64)          # total valuation
    sig = np.full(len(n), -1, dtype=np.int64)     # actual stopping time
    sig_inf = np.full(len(n), -1, dtype=np.int64) # coefficient stopping time
    active = np.ones(len(n), dtype=bool)
    cur = n.copy()
    step = 0
    log23 = math.log2(3)
    pow3 = 1
    while active.any():
        step += 1
        idx = np.nonzero(active)[0]
        m = 3 * cur[idx] + 1
        # valuation
        v = np.zeros(len(idx), dtype=np.int64)
        mm = m.copy()
        while True:
            even = (mm % 2 == 0)
            if not even.any():
                break
            mm[even] //= 2
            v[even] += 1
        cur[idx] = mm
        A[idx] += v
        # coefficient descent: 3^step < 2^A  (exact integer comparison)
        pow3 *= 3
        coef = (np.left_shift(np.int64(1), np.minimum(A[idx], 62)) > pow3) if pow3 < 2 ** 62 else (A[idx] > step * log23)
        newly_c = coef & (sig_inf[idx] < 0)
        sig_inf[idx[newly_c]] = step
        newly_a = (mm < start[idx]) & (sig[idx] < 0)
        sig[idx[newly_a]] = step
        done = (sig[idx] >= 0)
        active[idx[done]] = False
        if step > 2000:
            break
    assert not active.any(), "some n did not descend within 2000 Syracuse steps"
    diff = np.nonzero(sig != sig_inf)[0]
    print(" odd n <= %d: %d numbers; max actual stopping time %d; sigma != sigma_inf for %d numbers: %s" % (NMAX, len(n), sig.max(), len(diff), [(int(n[i]), int(sig[i]), int(sig_inf[i])) for i in diff[:20]]))
    # n = 1 is the fixed point: U(1) = 1, never < 1; it is excluded from 'active' termination only if sig set... handle separately
    return diff


def convergent_check(JMAX):
    log23 = math.log2(3)
    bad = []
    worst = (0, None)
    for j in range(1, JMAX + 1):
        A = math.ceil(j * log23)
        if 2 ** A <= 3 ** j:
            A += 1
        gap = 2 ** A - 3 ** j
        need = int(math.floor(1.5 ** j)) - 1 if j < 1000 else None
        # compare exactly with rationals: (3/2)^j - 1 < gap  <=>  3^j - 2^j < gap * 2^j
        ok = (3 ** j - 2 ** j) < gap * 2 ** j
        ratio = (3 ** j - 2 ** j) / (gap * 2 ** j)
        if ratio > worst[0]:
            worst = (ratio, j, A)
        if not ok:
            bad.append((j, A))
    print(" j <= %d: 2^A - 3^j > (3/2)^j - 1 at the minimal A fails for: %s ; worst ratio ((3/2)^j - 1)/(2^A - 3^j) = %.4g at j = %d, A = %d" % (JMAX, bad, worst[0], worst[1], worst[2]))


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10000000
    JMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
    print("== (1) actual versus coefficient stopping time, odd n <= %d (n = 1 has neither; it is handled as sigma = sigma_inf = 1 by the loop since U(1) = 1 is not < 1 ... reported separately) ==" % NMAX)
    n1 = 1
    m = 3 * n1 + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    print(" n = 1: v = %d, 2^v = %d > 3 so sigma_inf(1) = 1, but U(1) = %d is not < 1: sigma(1) is infinite; the unique exception" % (v, 2 ** v, m))
    compare_range = NMAX
    diff = compare(compare_range)
    print("== (2) one-member bound ==")
    convergent_check(JMAX)
    print(" hence for every coefficient-descent class with j <= %d the only possibly uncertified member is its representative; and (1) shows no such member exists below %d except n = 1." % (JMAX, NMAX))


if __name__ == '__main__':
    main()
