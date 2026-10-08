"""Coalescence profile of the Mersenne line and the deep barriers (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Uses the engines of the two companion scripts.
  (L) local hardness below the t = 23 source: odd exponents in [TOP - 19999, TOP] with no exit (D <= 64) within W = 4096.
  (Q) deep exits for the 30 progression members that need W > 16384 (fast parity, D <= 256, W <= 2^22), and the
      negative for t = 3525143 within 2^26.
  (S) coalescence profiles below six random odd exponents in [10^10, 10^11] (seed 77), deletions D <= 2048, at
      T = 2^j, j = 4..20: contiguous absorbed extent S(T), absorbed cluster size N(T) (number of D <= 2048 whose chains are
      absorbed by T, i.e. the source's cluster below it), and the number of unabsorbed groups (NUMERICAL).
  (H2) exact horizon clusters for K <= 12800 from mac-mini's runcompress_20261007/mersenne_sigma_12800.txt: by the
      stopping-time invariant (note section 4) the clusters are the level sets of (sigma_T(M_K) - K, odd count of M_K);
      per window: median horizon T = sigma_T(x_K), median contiguous extent S and cluster size below K over sqrt(T),
      orphan fraction per odd K against 2 / (mean cluster size), and the even cluster bottoms.
  Level 2 below the fan is in mersenne_fan_block_20261007.py (merge at 4,474,989) and mersenne_fan_level2_20261007.py (2^29).
A negative "within W" means: no absorption at any Terras time <= W - 9 (<= W - 8 in deep_exit, which also checks the final state).
Prints ALL CHECKS PASSED.  Runtime about 6 minutes; memory about 3 GB."""
import sys, time, random, os, math, statistics
from collections import defaultdict
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
        for D, (k, av) in alive.items():          # the state after the last coin
            if k == 0 and av == 0:
                return LOGW, n, sorted(groups[D])
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
            out.append((n + 1, S, len(abs_time), len(alive)))
            nextck *= 2
    return E, out


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
        meds = {}; medN = {}
        for E, out in profs:
            print(f'     E = {E}: (T, S(T), N(T), unabsorbed groups) = {out}')
            for T, S, Nc, g in out:
                meds.setdefault(T, []).append(S); medN.setdefault(T, []).append(Nc)
        med = {T: sorted(v)[len(v) // 2] for T, v in meds.items()}
        mdN = {T: sorted(v)[len(v) // 2] for T, v in medN.items()}
        print('     upper median S(T):', {T: med[T] for T in sorted(med) if T >= 4096})
        print('     upper median N(T):', {T: mdN[T] for T in sorted(mdN) if T >= 4096})
        print('     upper median N(T)/sqrt(T):', {T: round(mdN[T] / math.sqrt(T), 3) for T in sorted(mdN) if T >= 4096})
        print('     unabsorbed groups at 2^20:', [out[-1][3] for E, out in profs])
        check(med[1 << 20] >= 2 * med[1 << 12] and med[1 << 20] <= 64 * med[1 << 12] + 64 and mdN[1 << 20] > mdN[1 << 12],
              'median absorbed extent and cluster size grow between T = 2^12 and 2^20 at a sub-linear rate (NUMERICAL)')
    print('(Q2) the remaining member at 2^26')
    t = 3525143
    r = deep_exit(924745897 + (1 << 32) * t, logWs=(26,))
    check(r[0] is None, f't = {t}: no exit with D <= 256 within 2^26')
    print('(H2) exact horizon clusters for K <= 12800 (mac-mini table)')
    rows = {}
    with open(os.path.join(HERE, 'runcompress_20261007', 'mersenne_sigma_12800.txt')) as f:
        for line in f:
            q = line.split()
            if len(q) >= 3:
                rows[int(q[0])] = (int(q[1]), int(q[2]))
    key = {K: (t - K, o) for K, (o, t) in rows.items()}
    cl = defaultdict(list)
    for K in sorted(rows):
        cl[key[K]].append(K)
    bottom = {min(v) for v in cl.values()}
    even_bottoms = sorted(K for K in bottom if K % 2 == 0 and K >= 4)
    print(f'     {len(cl)} clusters among K = 2..12800; even cluster bottoms with K >= 4: {even_bottoms}')
    edges = [400, 800, 1600, 3200, 6400, 12801]
    print('     window | odd K | median T | median S/sqrtT | median N_below/sqrtT | orphan fraction x sqrtK | 2/(mean size) x sqrtK')
    for lo, hi in zip(edges, edges[1:]):
        Ss = []; Ns = []; Ts = []
        oddK = [K for K in range(lo | 1, hi, 2) if K in rows]
        for K in oddK:
            Tk = rows[K][1] - (K - 1)
            S = 0
            while (K - S - 1) in key and key[K - S - 1] == key[K]:
                S += 1
            Nb = sum(1 for K2 in cl[key[K]] if K2 < K)
            Ss.append(S / math.sqrt(Tk)); Ns.append(Nb / math.sqrt(Tk)); Ts.append(Tk)
        orph = sum(1 for K in oddK if K in bottom) / len(oddK)
        sizes = [len(v) for v in cl.values() if lo <= min(v) < hi]
        kmid = math.sqrt(lo * hi)
        print(f'     [{lo},{hi}) {len(oddK):5d} {statistics.median(Ts):9.0f} {statistics.median(Ss):8.3f} {statistics.median(Ns):8.3f} '
              f'{orph * math.sqrt(kmid):8.3f} {2 / statistics.mean(sizes) * math.sqrt(kmid):8.3f}')
    check(not even_bottoms, 'every cluster bottom with K >= 4 is odd (so the orphan fraction per odd K is about 2 / (mean cluster size))')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
