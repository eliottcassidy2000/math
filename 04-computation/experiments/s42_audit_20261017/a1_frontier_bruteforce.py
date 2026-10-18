#!/usr/bin/env python3
"""S42 audit a1: brute-force frontier coefficients c_Y(u), E_Y(N), R(N) from the DEFINITION, compared with the session's
incremental scripts.

Definition used here (independent of the session's K_N formula and three-entry update):
  R(N) = sum_{N < n <= 4N} w(n) f(n),  w(n) = 2 if n = 2 mod 3 else -1,  f(n) = (-1)^(first-hit shortcut clock of n).
  Each n is folded to its odd part m = n / 2^k (f(n) = (-1)^k f(m)), then to seed(m) (strip m -> (m-1)/4 while m = 5 mod 8,
  which preserves f), then transported s -> P(s) = seed((3s+1)/2^a) with sign (-1)^a until the first seed <= Y.
  c_Y(u) = sum over n landing at u of w(n) (-1)^k sigma(path).  E_Y = sum c_Y(u)^2,  R_Y = sum f(u) c_Y(u).

Compared with root_transport_frontier_energy_20261017.run (per-N arrays) for every N < 2^12 and for random N < 2^20,
and with root_transport_frontier_energy_big_20261017.run (octave statistics) for octaves 2^8..2^11.
Numba on-disk caching of the session modules is disabled (_nocache.py), so no session file is touched.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
import numpy as np
from numba import njit

EXP = os.path.dirname(HERE)


@njit(cache=False)
def fpar(n):
    c = 0
    while n != 1:
        if n & 1:
            n = (3 * n + 1) >> 1
        else:
            n >>= 1
        c += 1
    return 1 - 2 * (c & 1)


@njit(cache=False)
def wt(n):
    return 2 if n % 3 == 2 else -1


@njit(cache=False)
def to_seed(m):
    while m % 8 == 5:
        m = (m - 1) // 4
    return m


@njit(cache=False)
def brute(N, Y):
    c = np.zeros(Y + 1, dtype=np.int64)
    rdir = 0
    bad_a = 0
    for n in range(N + 1, 4 * N + 1):
        wn = wt(n)
        m = n
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        sg = 1 - 2 * (k & 1)
        x = to_seed(m)
        while x > Y:
            v = 3 * x + 1
            a = 0
            while v % 2 == 0:
                v //= 2
                a += 1
            if a != 1 and a != 2:
                bad_a += 1
            if a % 2 == 1:
                sg = -sg
            x = to_seed(v)
        c[x] += wn * sg
        rdir += wn * fpar(n)
    E = 0
    R = 0
    for u in range(1, Y + 1):
        if c[u] != 0:
            E += c[u] * c[u]
            R += c[u] * fpar(u)
    return c, E, R, rdir, bad_a


def main():
    sys.path.insert(0, EXP)
    import root_transport_frontier_energy_20261017 as S1
    import root_transport_frontier_energy_big_20261017 as S2
    sys.path.insert(0, HERE)
    from _nocache import nocache
    nocache(S1); nocache(S2)
    rng = np.random.default_rng(42)
    report = {}
    t0 = time.time()
    for Y in (30, 100, 1000):
        # session per-N arrays up to 2^20
        Nmax = 2 ** 20 + 1
        M = 4 * Nmax
        fu, fs = S1.frontier_map(M, Y)
        fY = np.zeros(Y + 1, dtype=np.int64)
        for u in range(1, Y + 1):
            fY[u] = fpar(u)
        EY, RY = S1.run(Nmax, Y, fu, fs, fY)
        mism_E = 0; mism_R = 0; mism_Rdir = 0; bad_total = 0
        allN = list(range(1, 2 ** 12))
        randN = sorted(set(int(x) for x in rng.integers(2 ** 12, 2 ** 20, 60))) + [4155, 4351, 7910, 66463, 1048576]
        oct_stats = {}
        for N in allN + randN:
            c, E, R, rdir, bad = brute(N, Y)
            bad_total += bad
            if abs(EY[N] - E / N) > 1e-9 * max(1.0, E / N):
                mism_E += 1
            if RY[N] != R:
                mism_R += 1
            if R != rdir:
                mism_Rdir += 1
            if N < 2 ** 12:
                k = N.bit_length() - 1
                st = oct_stats.setdefault(k, [0.0, 0, 0.0, 1e300, 0])
                e = E / N
                if e > st[0]:
                    st[0] = e; st[1] = N
                st[2] += e; st[4] += 1
                st[3] = min(st[3], e)
        # session big-script octave stats on N < 2^12
        par = S2.clock_parity8(4 * 2 ** 12)
        fu2, fs2 = S2.frontier_map(4 * 2 ** 12, Y)
        ok, emax, eargmax, esum, emin, rmax, rargmax, cnt = S2.run(2 ** 12, Y, fu2, fs2, par)
        oct_cmp = []
        for k in range(8, 12):
            st = oct_stats[k]
            oct_cmp.append(dict(k=k, brute_max=st[0], brute_argmax=st[1], brute_mean=st[2] / st[4], brute_min=st[3],
                                big_max=float(emax[k]), big_argmax=int(eargmax[k]), big_mean=float(esum[k] / cnt[k]),
                                big_min=float(emin[k]),
                                agree=bool(abs(st[0] - emax[k]) < 1e-12 and st[1] == eargmax[k] and
                                           abs(st[2] / st[4] - esum[k] / cnt[k]) < 1e-9 and abs(st[3] - emin[k]) < 1e-12)))
        report[Y] = dict(N_checked=len(allN) + len(randN), mismatches_E=mism_E, mismatches_R_session_vs_brute=mism_R,
                         mismatches_R_brute_vs_direct_sum=mism_Rdir, transport_steps_with_a_not_1_or_2=bad_total,
                         big_script_invariant_ok=bool(ok), octave_comparison=oct_cmp,
                         sample_values={str(N): float(EY[N]) for N in (4155, 4351, 7910, 66463, 1048576)})
        print('Y = %d: %d N checked (all N < 4096 + %d random/special N < 2^20 + 1): E mismatches %d, R(session) vs brute %d, '
              'brute R vs direct H(4N)-H(N) %d, a not in {1,2}: %d [%.0f s]'
              % (Y, len(allN) + len(randN), len(randN), mism_E, mism_R, mism_Rdir, bad_total, time.time() - t0), flush=True)
        for r in oct_cmp:
            print('   octave 2^%d: brute max %.4f (N=%d) mean %.4f min %.4f | big-script max %.4f (N=%d) mean %.4f min %.4f | agree %s'
                  % (r['k'], r['brute_max'], r['brute_argmax'], r['brute_mean'], r['brute_min'], r['big_max'], r['big_argmax'],
                     r['big_mean'], r['big_min'], r['agree']), flush=True)
        print('   E_Y/N at N = 4155, 4351, 7910, 66463, 2^20: ' + ', '.join('%.4f' % v for v in report[Y]['sample_values'].values()))
    # N = 1 check of the E_Y(1) = 10 claim for Y >= 3 and the Y = 1, 2 values
    small = {}
    for Y in (1, 2, 3, 10, 100):
        c, E, R, rdir, bad = brute(1, Y)
        small[Y] = dict(E=int(E), R=int(R), R_direct=int(rdir))
    print('N = 1: ' + ', '.join('Y=%d: E=%d R=%d' % (Y, v['E'], v['R']) for Y, v in small.items()))
    report['N1'] = small
    json.dump(report, open(os.path.join(HERE, 'a1_frontier_bruteforce.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
