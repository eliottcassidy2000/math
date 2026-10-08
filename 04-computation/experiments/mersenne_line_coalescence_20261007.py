"""Coalescence profile of the Mersenne line and the deep barriers (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Uses the engines of the two companion scripts.
  (L) local hardness below the t = 23 source: odd exponents in [TOP - 19999, TOP] with no exit (D <= 64) within W = 4096.
  (Q) deep exits for the 30 progression members that need W > 16384 (fast parity, D <= 256, W <= 2^22), and the
      negative for t = 3525143 within 2^26.
  (S) coalescence profiles below six random odd exponents in [10^10, 10^11] (seed 77): contiguous absorbed extent S(T)
      at T = 2^j, j = 4..21, for deletions D <= 2048 (NUMERICAL).
  (V) level 2 below the fan: the chains for D = 1911, 1912, 2500, 4000 of the fan bottom are not absorbed within 2^27.
Prints ALL CHECKS PASSED.  Runtime about 12 minutes; memory about 3 GB."""
import sys, time, random, importlib.util, os
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))


sys.path.insert(0, HERE)
import mersenne_line_barriers_20261007 as B          # noqa: E402
import mersenne_fan_escape_20261007 as F             # noqa: E402
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


TOP = 99708993705


def job_local(j):
    E = TOP - j
    if E % 2 == 0:
        return j, True                       # reset rule
    return j, bool(B.exits(E, 4096, 64))


def deep_exit(E, logWs=(20, 22), DMAX=256):
    for LOGW in logWs:
        cb = F.coins_of(E, (1 << LOGW) - 8)
        alive = {D: (-D, 1 - F.P3[D]) for D in range(1, DMAX + 1)}
        groups = {D: [D] for D in alive}
        n = 0; N = len(cb)
        while alive and n < N:
            for b in cb[n:n + 65536].tolist():
                new = {}
                for D, (k, av) in alive.items():
                    if k == 0 and av == 0:
                        return LOGW, n, sorted(groups[D])
                    new[D] = F.stepf(k, av, b)
                seen = {}; alive = {}
                for D, st in new.items():
                    if st in seen:
                        groups[seen[st]].extend(groups.pop(D))
                    else:
                        seen[st] = D; alive[D] = st
                n += 1
                if n >= N:
                    break
    return None, None, None


def job_deep(t):
    E = 924745897 + (1 << 32) * t
    return t, deep_exit(E)


def job_profile(E, LOGW=21, DMAX=2048):
    cb = F.coins_of(E, (1 << LOGW) - 8).tolist()
    alive = {D: (-D, 1 - F.P3[D]) for D in range(1, DMAX + 1)}
    members = {D: [D] for D in alive}
    abs_time = {}
    out = []; nextck = 16
    for n, b in enumerate(cb):
        new = {}
        for D, (k, av) in alive.items():
            if k == 0 and av == 0:
                for m in members[D]:
                    abs_time[m] = n
                continue
            new[D] = F.stepf(k, av, b)
        seen = {}; alive = {}
        for D, st in new.items():
            if st in seen:
                members[seen[st]].extend(members.pop(D))
            else:
                seen[st] = D; alive[D] = st
        if n + 1 == nextck:
            S = 0
            while S + 1 <= DMAX and (S + 1) in abs_time:
                S += 1
            out.append((n + 1, S, len(alive)))
            nextck *= 2
    return E, out


COINS2 = None


def init_level2():
    global COINS2
    COINS2 = F.coins_of(99708993677, (1 << 27) - 8)


def job_level2(D):
    k, av = -D, 1 - F.P3[D]; kmin = k
    n = 0; N = len(COINS2)
    while n < N:
        for b in COINS2[n:n + 65536].tolist():
            if k == 0 and av == 0:
                return D, n, kmin
            k, av = F.stepf(k, av, b)
            if k < kmin:
                kmin = k
            n += 1
    return D, None, kmin


if __name__ == '__main__':
    t0 = time.time()
    with Pool(12) as pool:
        print('(L) local hardness below the t = 23 source')
        loc = pool.map(job_local, range(20000), chunksize=50)
        hard = [j for j, ok in loc if not ok]
        check(all((TOP - j) % 2 == 1 for j in hard) and 0.02 < len(hard) / 20000 < 0.04,
              f'{len(hard)} of 20000 exponents ({len(hard) / 20000:.4f}) have no exit with D <= 64 within 4096, all odd; first offsets {hard[:8]}')
        print('(Q) deep exits for the progression members needing W > 16384')
        mem = pool.map(B.job_member, range(2000), chunksize=10)
        need = [m[1] for m in mem if m[2] is None]
        deep = dict(pool.map(job_deep, need))
        got = {t: r for t, r in deep.items() if r[0] is not None}
        for t in need:
            r = deep[t]
            print(f'     t = {t}: ' + (f'absorbed at {r[1]} (W = 2^{r[0]}), absorbed group max D {r[2][-1]}, size {len(r[2])}' if r[0] else 'no exit with D <= 256 within 2^22'))
        miss = [t for t in need if deep[t][0] is None]
        check(len(need) == 30 and len(got) == 29 and miss == [3525143],
              f'{len(got)} of the {len(need)} members needing W > 16384 are certified within 2^22; the remaining one is t = {miss}')
        print('(S) coalescence profiles below six random odd exponents')
        rng = random.Random(77)
        Es = [rng.randrange(10 ** 10, 10 ** 11) | 1 for _ in range(6)]
        profs = pool.map(job_profile, Es)
        meds = {}
        for E, out in profs:
            print(f'     E = {E}: (T, S(T), live groups) = {out}')
            for T, S, g in out:
                meds.setdefault(T, []).append(S)
        med = {T: sorted(v)[len(v) // 2] for T, v in meds.items()}
        print('     median S(T):', {T: med[T] for T in sorted(med) if T >= 4096})
        check(med[1 << 20] >= 2 * med[1 << 12] and med[1 << 20] <= 64 * med[1 << 12] + 64,
              'median absorbed extent grows between T = 2^12 and 2^20 at a sub-linear rate (NUMERICAL)')
    print('(Q2) the remaining member at 2^26')
    t = 3525143
    r = deep_exit(924745897 + (1 << 32) * t, logWs=(26,))
    check(r[0] is None, f't = {t}: no exit with D <= 256 within 2^26')
    print('(V) level 2 below the fan, W = 2^27')
    with Pool(4, initializer=init_level2) as pool:
        lv = pool.map(job_level2, [1911, 1912, 2500, 4000])
    print('     ', lv)
    check(all(n is None for D, n, kmin in lv), 'the chains for D = 1911, 1912, 2500, 4000 of the fan bottom are not absorbed within 2^27 steps')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
