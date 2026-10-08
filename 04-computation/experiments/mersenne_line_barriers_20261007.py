"""Mersenne-line exit certificates at scale (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Deep parts: mersenne_fan_escape_20261007.py,
mersenne_fan_block_20261007.py, mersenne_line_coalescence_20261007.py, mersenne_fan_level2_20261007.py.

M_E = 2^E - 1.  After E - 1 Terras steps it is x = 2*3^(E-1) - 1; the child M_(E-D) is y = 2*3^(E-D-1) - 1 and
y + 1 = 3^-D (x + 1).  Source-reference pair chain (THM-4581's table, reference orbit x): state (k, a), a = 3^max(0,-k) e,
start (-D, 1 - 3^D); absorption (0, 0) at Terras time n proves T^n(x) = T^n(y), an exact common future, so
M_E ~> M_(E-D) (M_E reaches 1 iff M_(E-D) does).  Only x mod 2^W is needed (W > n).
Checks:
  (A) the chain against direct big-integer orbits (large E: no cycle coincidences in the window); small-E disagreements
      are exactly trivial-cycle value coincidences, which the relation chain correctly does not claim.
  (A2) no-meeting lemma (PROVED in the note; checked here): for E in [1500, 1520), D <= 40, N = 600 (so E >= 2.27 N + D + 2),
      if the chain is not absorbed at any time <= N the first N + 1 iterates of x and y are disjoint, and if it is first
      absorbed at n every common value occurs at equal times >= n.
  (B) reset rule (THM-4601 (i) / THM-4556): for even E, T^3(x) = T^3(y_1) (M_E ~> M_(E-1) in 3 steps); never for odd E.
  (C) the supplied residual t = 23 (E = 99708993705): exits are exactly D in {3, 4, 7, ..., 28}, all at Terras time 6553
      (codex-grounding's 24 children); the original parent M_99708993713 exits to M_99708993677 (D = 36) at 6553 (all its
      exits with D <= 64 within W = 8000 are printed with their times).
  (D) the fan: within W = 16384, D <= 128 the exits of each supplied child are exactly the fan members below it; the bottom
      M_99708993677 has none (FINITE-EXACT negative within the bound); the internal merge times are printed.
  (E) the unpaid progression t = 13847 mod 27648 (E = 924745897 + 2^32 t): certified exits for its first 2000 members,
      W <= 16384 by direct coins (split by the first successful W); members needing more are listed (deep script).
  (F) random giant odd exponents: tail of the cheapest exit depth (NUMERICAL).
  (M) stopping-time invariant (PROVED in the note; checked here): for 4 <= E <= 160 and 1 <= D <= E - 2, a chain absorption
      with merge value >= 3 occurs iff sigma_T(x_E) = sigma_T(x_(E-D)) and M_E, M_(E-D) have equal odd-step counts to 1;
      absorptions at the trivial cycle are counted separately.  Exact sigma_T(x_K) for K <= 2000 (cross-checked against
      mac-mini's runcompress_20261007/mersenne_sigma_12800.txt) and the exponent E* above which no route of such
      certificates can reach any K <= 2000 (resp. K <= 12800, from mac-mini's table); no pair K' < K <= 12800 qualifies for an
      absorption on the trivial cycle (equal 2o - s with different s).
Prints ALL CHECKS PASSED.  Runtime about 30 seconds on 12 cores."""
import sys, os, time, random, math, statistics
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


class _P3(list):
    def __getitem__(self, i):
        while i >= len(self):
            self.append(super().__getitem__(len(self) - 1) * 3)
        return super().__getitem__(i)


P3 = _P3([3 ** i for i in range(2000)])


def T(z):
    return (3 * z + 1) // 2 if z & 1 else z // 2


def src_coins(E, W):
    z = (2 * pow(3, E - 1, 1 << W) - 1) % (1 << W)
    out = bytearray(W - 8)
    for n in range(W - 8):
        b = z & 1
        out[n] = b
        z = ((3 * z + 1) >> 1) if b else (z >> 1)
    return out


def chain(D, cb, nmax):
    k = -D; a = 1 - P3[D]
    for n in range(nmax):
        if k == 0 and a == 0:
            return n
        if abs(k) > nmax - n:
            return None
        b = cb[n]
        if a & 1 == 0:
            if b == 0:
                a >>= 1
            elif k >= 0:
                a = (3 * a + 1 - P3[k]) >> 1
            else:
                a = (3 * a + P3[-k] - 1) >> 1
        else:
            if b == 0:
                a = ((3 * a + 1) >> 1) if k >= 0 else ((a + P3[-k - 1]) >> 1)
                k += 1
            else:
                a = ((a - P3[k - 1]) >> 1) if k >= 1 else ((3 * a - 1) >> 1)
                k -= 1
    return None


def exits(E, W, Dmax):
    cb = src_coins(E, W)
    return {D: n for D in range(1, min(Dmax, E - 2) + 1) for n in [chain(D, cb, W - 8)] if n is not None}


# ---------------------------------------------------------------- workers
def job_direct(args):
    E, D, W = args
    x = 2 * 3 ** (E - 1) - 1; y = 2 * 3 ** (E - D - 1) - 1
    n = chain(D, src_coins(E, W), W - 8)
    a, b = x, y; dn = None; val = None
    for i in range(W - 8):
        if a == b:
            dn = i; val = a; break
        a, b = T(a), T(b)
    return E, D, n, dn, val


def job_nomeet(args):
    E, D, N = args
    x = 2 * 3 ** (E - 1) - 1; y = 2 * 3 ** (E - D - 1) - 1
    n = chain(D, src_coins(E, N + 9), N + 1)          # absorption checked at every time 0..N
    px = {}; py = {}
    a, b = x, y
    for i in range(N + 1):
        px.setdefault(a, i); py.setdefault(b, i)
        a, b = T(a), T(b)
    common = set(px) & set(py)
    if n is None:
        return E, D, n, not common
    return E, D, n, bool(common) and all(px[v] == py[v] >= n for v in common) and min(px[v] for v in common) == n


def sig_odd(z):
    n = 0; o = 0
    while z != 1:
        if z & 1:
            z = (3 * z + 1) >> 1; o += 1
        else:
            z >>= 1
        n += 1
    return n, o


def job_sigma(K):
    return K, sig_odd(2 * 3 ** (K - 1) - 1)


def job_invariant(E):
    x = 2 * 3 ** (E - 1) - 1
    sx, _ = sig_odd(x)
    orbit = [x]; z = x
    for _ in range(sx + 12):
        z = T(z); orbit.append(z)
    cb = bytes(v & 1 for v in orbit)
    out = []
    for D in range(1, E - 1):
        n = chain(D, cb, sx + 10)
        out.append((D, n, orbit[n] if n is not None else None))
    return E, out


def job_child(D):
    E = 99708993705
    return D, exits(E - D, 16384, 128)


def job_member(i):
    t = 13847 + 27648 * i
    E = 924745897 + (1 << 32) * t
    for W in (256, 1024, 4096, 16384):
        r = exits(E, W, 64)
        if r:
            n = min(r.values())
            return i, t, W, n, max(D for D, m in r.items() if m == n)
    return i, t, None, None, None


def job_random(seed):
    rng = random.Random(1000 + seed)
    E = rng.randrange(10 ** 10, 10 ** 11) | 1
    r = exits(E, 4096, 64)
    return min(r.values()) if r else None


if __name__ == '__main__':
    t0 = time.time()
    with Pool(12) as pool:
        print('(A) pair chain vs direct orbits')
        tasks = [(E, D, 400) for E in range(200, 260) for D in range(1, 41)]
        res = pool.map(job_direct, tasks, chunksize=40)
        check(all(n == dn for E, D, n, dn, v in res), f'large E (200..259, D <= 40, 392 steps): chain absorption time = first equal-time coincidence of the direct orbits in all {len(res)} cases ({sum(1 for r in res if r[2] is not None)} merges)')
        small = pool.map(job_direct, [(E, D, 400) for E in range(8, 40) for D in range(1, E - 2)], chunksize=40)
        mism = [r for r in small if r[2] != r[3]]
        check(all(r[4] is not None and r[4] <= 2 and (r[2] is None or r[2] > r[3]) for r in mism),
              f'small E (8..39): all {len(mism)} disagreements are trivial-cycle value coincidences (merge value <= 2) that the relation chain does not claim')

        print('(A2) no-meeting lemma')
        N = 600
        nm = pool.map(job_nomeet, [(E, D, N) for E in range(1500, 1520) for D in range(1, 41)], chunksize=20)
        cond = all(E >= N * (1 + math.log(4, 3)) + D + 2 for E, D, n, ok in nm)
        check(cond and all(ok for E, D, n, ok in nm),
              f'E in [1500, 1520), D <= 40, N = {N} (E >= 2.27 N + D + 2): the {sum(1 for r in nm if r[2] is None)} chains not absorbed by N have '
              f'disjoint orbit segments (first N + 1 iterates of x and y); the {sum(1 for r in nm if r[2] is not None)} absorbed chains meet '
              f'only at equal times, first at the absorption time')

        print('(B) reset rule on the Mersenne line')
        T3 = lambda z: T(T(T(z)))
        check(all(T3(2 * 3 ** (E - 1) - 1) == T3(2 * 3 ** (E - 2) - 1) for E in range(4, 3000, 2)) and
              not any(T3(2 * 3 ** (E - 1) - 1) == T3(2 * 3 ** (E - 2) - 1) for E in range(5, 3000, 2)),
              'for even E < 3000, M_E ~> M_(E-1) at Terras time 3 (w = 3^(E-2) = 1 mod 8: both reach (9w-1)/4); for odd E never')

        print('(C) the supplied residual t = 23')
        r23 = exits(99708993705, 7200, 40)
        check(r23 == {D: 6553 for D in [3, 4] + list(range(7, 29))}, 'E = 99708993705: exits exactly D in {3,4,7,...,28}, all absorbed at Terras time 6553 (reproduces the 24 children)')
        rp = exits(99708993713, 8000, 64)
        check(rp.get(36) == 6553 and rp.get(8) is not None,
              f'original parent E = 99708993713: D = 36 (to the fan bottom 99708993677) absorbed at 6553; D = 8 (to the source) at {rp.get(8)}')
        print('     all exits of the parent with D <= 64 within W = 8000 (D: Terras time):', dict(sorted(rp.items())))

        print('(D) the fan of the 24 children')
        CH = [3, 4] + list(range(7, 29))
        fan = set(99708993705 - D for D in CH) | {99708993705}
        cres = dict(pool.map(job_child, CH))
        inside = all(all((99708993705 - D) - D2 in fan for D2 in r) for D, r in cres.items())
        below = all(set(r) == {D2 - D for D2 in CH if D2 > D} for D, r in cres.items())
        check(inside and below and cres[28] == {},
              'within W = 16384 and D <= 128 the exits of each child are exactly the fan members below it, and the bottom M_99708993677 has no exit')
        print('     exits per child (D: count):', {D: len(r) for D, r in cres.items()})
        times = Counter(n for r in cres.values() for n in r.values())
        print('     merge times of the exits among the children (Terras time: number of exits):', dict(sorted(times.items())))

        print('(E) the unpaid progression t = 13847 mod 27648: first 2000 members')
        mem = pool.map(job_member, range(2000), chunksize=10)
        ok = [m for m in mem if m[2] is not None]
        hard = [(m[0], m[1]) for m in mem if m[2] is None]
        ns = sorted(m[3] for m in ok)
        check(len(ok) >= 1950, f'{len(ok)} of 2000 members certified with W <= 16384 (median cheapest merge time {ns[len(ns) // 2]}, max {ns[-1]}); {len(hard)} need the deep script')
        print('     members needing W > 16384 (index, t):', hard)
        print('     first members:', [(m[1], m[4], m[3]) for m in mem[:4]], '(t, D, merge time)')
        print('     certified members by first successful W:', dict(sorted(Counter(m[2] for m in ok).items())))

        print('(F) random giant odd exponents: cheapest exit depth (W = 4096, D <= 64)')
        rr = pool.map(job_random, range(600), chunksize=10)
        for Wc in (128, 256, 512, 1024, 2048, 4096):
            print(f'     P(no exit by {Wc}) = {sum(1 for v in rr if v is None or v > Wc) / len(rr):.3f}')

        print('(M) the stopping-time invariant')
        sg = dict(pool.map(job_sigma, range(2, 2001), chunksize=20))
        inv = pool.map(job_invariant, range(4, 161))
        nabs = nhi = ncyc = bad = maxcyc = 0
        for E, out in inv:
            for D, n, v in out:
                K = E - D
                same = sg[E][0] == sg[K][0] and (E - 1) + sg[E][1] == (K - 1) + sg[K][1]
                hi = n is not None and v >= 3
                if n is not None:
                    nabs += 1
                    if v >= 3:
                        nhi += 1
                    else:
                        ncyc += 1; maxcyc = max(maxcyc, abs(sg[E][0] - sg[K][0]))
                bad += hi != same
        check(bad == 0, f'4 <= E <= 160, D <= E - 2: {nabs} chain absorptions; the {nhi} with merge value >= 3 are exactly the pairs with '
                        f'sigma_T(x_E) = sigma_T(x_(E-D)) and equal odd-step counts of M_E, M_(E-D) to 1; {ncyc} absorptions at the trivial cycle '
                        f'(sigma_T differs there by up to {maxcyc})')
        mm = {}
        with open(os.path.join(HERE, 'runcompress_20261007', 'mersenne_sigma_12800.txt')) as f:
            for line in f:
                q = line.split()
                if len(q) >= 3:
                    mm[int(q[0])] = (int(q[1]), int(q[2]))
        okm = all(mm[K] == ((K - 1) + sg[K][1], (K - 1) + sg[K][0]) for K in range(2, 2001))
        smax = max(sg[K][0] for K in range(2, 2001))
        smax2 = max(t - (K - 1) for K, (o, t) in mm.items())
        Es1 = 1 + smax / math.log2(3); Es2 = 1 + smax2 / math.log2(3)
        check(okm, f'sigma_T(x_K) for K <= 2000 agrees with the mac-mini table; max sigma_T(x_K) = {smax} for K <= 2000 and {smax2} for K <= 12800; '
                   f'since sigma_T(x_E) > (E - 1) log2 3, no route of certificates with merge values >= 3 from any E > {Es1:.0f} '
                   f'(resp. E > {Es2:.0f}) reaches any K <= 2000 (resp. K <= 12800)')
        big = [K for K in mm if K >= 1000]
        c = statistics.mean((mm[K][1] - (K - 1)) / K for K in big)
        sd = statistics.pstdev([((mm[K][1] - (K - 1)) - c * K) / math.sqrt(K) for K in big])
        cls = {}
        for K, (o, t) in mm.items():
            cls.setdefault((t - K, o), []).append(K)
        spans = [(max(v) - min(v) + 1) / math.sqrt(min(v)) for v in cls.values() if min(v) >= 1000]
        c33 = sorted(set(min(v) for v in cls.values() if any(33 <= K <= 42 for K in v)))
        print(f'     residual horizon: mean sigma_T(x_K)/K = {c:.4f} over 1000 <= K <= 12800 (6.95212 ln 3 = {6.95212 * math.log(3):.4f}); '
              f'std of (sigma_T(x_K) - cK)/sqrt(K) = {sd:.3f}; cluster span/sqrt(bottom) for bottoms >= 1000: mean {statistics.mean(spans):.2f}, max {max(spans):.2f} '
              f'({len(spans)} clusters); clusters meeting K = 33..42 have bottoms {c33}')
        grp = {}
        for K, (o, t) in mm.items():
            grp.setdefault(2 * o - (t - (K - 1)), []).append(t - (K - 1))
        ncyc_pairs = 0
        for v in grp.values():
            n_all = len(v) * (len(v) - 1) // 2
            ncyc_pairs += n_all - sum(m * (m - 1) // 2 for m in Counter(v).values())
        check(ncyc_pairs == 0, f'trivial-cycle absorptions (equal 2o - s, different s) among K <= 12800: {ncyc_pairs} pairs')
        z = (1 << 100000) - 1; st = od = 0
        while z != 1:
            if z & 1:
                z = (3 * z + 1) >> 1; od += 1
            else:
                z >>= 1
            st += 1
        check(st == 863323 and od == 481603, f'M_100000 = 2^100000 - 1 (Ren 2018): {st} Terras steps (= halvings), {od} odd steps; sigma_T/K = {st / 100000:.4f}')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
