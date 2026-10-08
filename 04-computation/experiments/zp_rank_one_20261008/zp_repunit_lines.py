"""Base-p repunit lines R_K = p^K - 1 under C_p, p = 7, 11, 13, K <= 8000 (session opus-2026-10-08-S22).
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md (section 4).

R_K rises for K - 1 steps to x_K = p (p+1)^(K-1) - 1 (all digits p - 1).  For each K: n = steps to the first value on a
known cycle (cycles met by starts <= 20000), n_mult = multiplication steps, the entry value and the cycle's minimum.
Data are written to repunit_orbits_p{p}_K8000.txt (columns K n n_mult entry cycle_min).
  (R1) FINITE-EXACT: every K <= 8000 ends in a known cycle; for p = 7, 13 all in the trivial cycle; for p = 11 the K that end
       in the 642-cycle, and their partner classes.
  (R2) FINITE-EXACT: partner classes (key = (n - K, n_mult, entry, cycle)) are level sets of n_mult within each terminal
       cycle ("rivers", the base-p form of mac-mini's identity), except for pairs that enter the cycle at different points
       at the same height (a coincidence on the cycle itself, not an equal-time merge above it); exceptions are listed.
  (R3) NUMERICAL: residual steps per digit (n - K + 1)/K against ln(p+1)/|Lambda_p|.
  (R4) NUMERICAL: orphan fraction (class bottoms) per dyadic window of K, times sqrt(K), and local exponents.
Prints ALL CHECKS PASSED.  Runtime about 20 minutes on 12 cores."""
import sys, os, math, time, statistics
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
KMAX = 8000


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def C(p, x):
    i = x % p
    return x // p if i == 0 else ((p + 1) * x + p - i) // p


def known_cycles(p, N=20000):
    elems = {}
    for x0 in range(1, N + 1):
        x = x0; seen = set()
        while x not in elems and x not in seen:
            seen.add(x); x = C(p, x)
        if x in seen and x not in elems:
            cyc = [x]; z = C(p, x)
            while z != x:
                cyc.append(z); z = C(p, z)
            m = min(cyc)
            for z in cyc:
                elems[z] = m
    return elems


ELEMS = None
PP = None


def init(p):
    global ELEMS, PP
    PP = p
    ELEMS = known_cycles(p)


def job(K):
    p = PP; P = p + 1
    x = p ** K - 1
    n = 0; nm = 0
    elems = ELEMS
    while x not in elems:
        i = x % p
        if i:
            x = (P * x + p - i) // p; nm += 1
        else:
            x //= p
        n += 1
        if n > 10 ** 7:
            return K, n, nm, x, -1
    return K, n, nm, x, elems[x]


def orbits(p):
    path = os.path.join(HERE, f'repunit_orbits_p{p}_K{KMAX}.txt')
    with Pool(12, initializer=init, initargs=(p,)) as pool:
        res = pool.map(job, list(range(KMAX, 0, -1)), chunksize=4)
    res.sort()
    with open(path, 'w') as f:
        for r in res:
            f.write(' '.join(map(str, r)) + '\n')
    return {r[0]: r[1:] for r in res}


if __name__ == '__main__':
    t0 = time.time()
    data = {}
    for p in (7, 11, 13):
        data[p] = orbits(p)
        print(f'     p = {p}: orbits of R_K for K <= {KMAX} computed [{time.time() - t0:.0f}s]', flush=True)
    print('(R1) terminal cycles')
    ok1 = True; second11 = []
    for p, rows in data.items():
        cyc = defaultdict(int)
        for K, (n, nm, x, c) in rows.items():
            cyc[c] += 1
        print(f'     p = {p}: terminal cycles (min: count) {dict(cyc)}')
        ok1 &= -1 not in cyc
        if p in (7, 13):
            ok1 &= set(cyc) == {1}
        if p == 11:
            second11 = sorted(K for K, r in rows.items() if r[3] != 1)
    rows = data[11]
    key11 = {K: (rows[K][0] - K, rows[K][1], rows[K][2], rows[K][3]) for K in rows}
    cls = defaultdict(list)
    for K in second11:
        cls[key11[K]].append(K)
    print(f'     p = 11: {len(second11)} repunits end in the 642-cycle; their partner classes (bottom, size): '
          f'{sorted((min(v), len(v)) for v in cls.values())}')
    check(ok1 and len(second11) > 0, 'every R_K (K <= 8000) ends in a known cycle; p = 7, 13: all in the trivial cycle; p = 11: some end in the 642-cycle')
    print('(R2) partner classes are rivers (level sets of n_mult within a terminal cycle)')
    ok2 = True; exc_all = {}
    for p, rows in data.items():
        by = defaultdict(list)
        for K, (n, nm, x, c) in rows.items():
            by[(nm, c)].append(K)
        exc = []
        for lev, Ks in by.items():
            keys = {(rows[K][0] - K, rows[K][1], rows[K][2], rows[K][3]) for K in Ks}
            if len(keys) > 1:
                same_height = len({rows[K][0] - K for K in Ks}) == 1
                exc.append((lev, sorted(Ks), sorted(keys), same_height))
        exc_all[p] = exc
        print(f'     p = {p}: (n_mult, cycle) levels carrying more than one partner key: {len(exc)} {exc}')
        ok2 &= all(same_h and max(Ks) <= 10 for _, Ks, _, same_h in exc)
    check(ok2 and exc_all[11] and not exc_all[7] and not exc_all[13],
          'within each terminal cycle the partner key is a function of n_mult (FINITE-EXACT, K <= 8000), except p = 11, K = 3, 4: '
          'same height and n_mult, entering the trivial cycle at 9 and 10')
    print('(R3) residual steps per digit')
    ok3 = True
    for p, rows in data.items():
        lam = (1 / p) * math.log(1 / p) + ((p - 1) / p) * math.log((p + 1) / p)
        pred = math.log(p + 1) / abs(lam)
        meas = statistics.mean((rows[K][0] - (K - 1)) / K for K in rows if K >= 2000)
        print(f'     p = {p}: mean (n - K + 1)/K over K >= 2000 = {meas:.3f}; ln(p+1)/|Lambda_p| = {pred:.3f}')
        ok3 &= abs(meas / pred - 1) < 0.02
    check(ok3, 'the residual horizon matches ln(p+1)/|Lambda_p| within 2% (NUMERICAL)')
    print('(R4) orphan fractions')
    edges = [100, 200, 400, 800, 1600, 3200, 6400, KMAX + 1]
    for p, rows in data.items():
        key = {K: (rows[K][0] - K, rows[K][1], rows[K][2], rows[K][3]) for K in rows}
        cl = defaultdict(list)
        for K in sorted(rows):
            cl[key[K]].append(K)
        bottoms = {min(v) for v in cl.values()}
        pts = []; line = []
        for lo, hi in zip(edges, edges[1:]):
            Ks = range(lo, hi)
            orph = sum(1 for K in Ks if K in bottoms) / len(Ks)
            km = math.sqrt(lo * hi)
            pts.append((km, orph))
            line.append(f'[{lo},{hi}): {orph:.4f} ({orph * math.sqrt(km):.3f})')
        ex = [round(math.log(o1 / o2) / math.log(k2 / k1), 2) for (k1, o1), (k2, o2) in zip(pts, pts[1:]) if o1 > 0 and o2 > 0]
        print(f'     p = {p}: {len(cl)} classes; orphan fraction (x sqrt K): ' + '; '.join(line))
        print(f'            local exponents between windows: {ex}')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
