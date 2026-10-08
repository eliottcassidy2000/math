"""Mersenne-line exit certificates at scale (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Deep part: mersenne_fan_escape_20261007.py.

M_E = 2^E - 1.  After E - 1 Terras steps it is x = 2*3^(E-1) - 1; the child M_(E-D) is y = 2*3^(E-D-1) - 1 and
y + 1 = 3^-D (x + 1).  Source-reference pair chain (THM-4581's table, reference orbit x): state (k, a), a = 3^max(0,-k) e,
start (-D, 1 - 3^D); absorption (0, 0) at Terras time n proves T^n(x) = T^n(y), an exact common future, so
M_E ~> M_(E-D) (M_E reaches 1 iff M_(E-D) does).  Only x mod 2^W is needed (W > n).
Checks:
  (A) the chain against direct big-integer orbits (large E: no cycle coincidences in the window); small-E disagreements
      are exactly trivial-cycle value coincidences, which the relation chain correctly does not claim.
  (B) reset rule (THM-4601 (i) / THM-4556): for even E, T^3(x) = T^3(y_1) (M_E ~> M_(E-1) in 3 steps); never for odd E.
  (C) the supplied residual t = 23 (E = 99708993705): exits are exactly D in {3, 4, 7, ..., 28}, all at Terras time 6553
      (codex-grounding's 24 children); the original parent M_99708993713 exits to M_99708993677 (D = 36) at 6553.
  (D) the fan: every exit of every supplied child (W = 16384, D <= 128) lands on another fan member; the bottom
      M_99708993677 has none (FINITE-EXACT negative within the bound).
  (E) the unpaid progression t = 13847 mod 27648 (E = 924745897 + 2^32 t): certified exits for its first 2000 members,
      W <= 16384 by direct coins; members needing more are listed (handled by the deep script).
  (F) random giant odd exponents: tail of the cheapest exit depth (NUMERICAL).
Prints ALL CHECKS PASSED.  Runtime about 3 minutes on 12 cores."""
import sys, time, random
from multiprocessing import Pool

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

        print('(B) reset rule on the Mersenne line')
        T3 = lambda z: T(T(T(z)))
        check(all(T3(2 * 3 ** (E - 1) - 1) == T3(2 * 3 ** (E - 2) - 1) for E in range(4, 3000, 2)) and
              not any(T3(2 * 3 ** (E - 1) - 1) == T3(2 * 3 ** (E - 2) - 1) for E in range(5, 3000, 2)),
              'for even E < 3000, M_E ~> M_(E-1) at Terras time 3 (w = 3^(E-2) = 1 mod 8: both reach (9w-1)/4); for odd E never')

        print('(C) the supplied residual t = 23')
        r23 = exits(99708993705, 7200, 40)
        check(r23 == {D: 6553 for D in [3, 4] + list(range(7, 29))}, 'E = 99708993705: exits exactly D in {3,4,7,...,28}, all absorbed at Terras time 6553 (reproduces the 24 children)')
        rp = exits(99708993713, 8000, 64)
        check(rp.get(36) == 6553, f'original parent E = 99708993713: D = 36 (to the fan bottom 99708993677) absorbed at 6553; its exits {sorted(rp)[:8]}...')

        print('(D) the fan of the 24 children')
        CH = [3, 4] + list(range(7, 29))
        fan = set(99708993705 - D for D in CH) | {99708993705}
        cres = dict(pool.map(job_child, CH))
        inside = all(all((99708993705 - D) - D2 in fan for D2 in r) for D, r in cres.items())
        check(inside and cres[28] == {}, 'within W = 16384 and D <= 128 every exit of every child lands on another fan member, and the bottom M_99708993677 has no exit')
        print('     exits per child (D: count):', {D: len(r) for D, r in cres.items()})

        print('(E) the unpaid progression t = 13847 mod 27648: first 2000 members')
        mem = pool.map(job_member, range(2000), chunksize=10)
        ok = [m for m in mem if m[2] is not None]
        hard = [(m[0], m[1]) for m in mem if m[2] is None]
        ns = sorted(m[3] for m in ok)
        check(len(ok) >= 1950, f'{len(ok)} of 2000 members certified with W <= 16384 (median cheapest merge time {ns[len(ns) // 2]}, max {ns[-1]}); {len(hard)} need the deep script')
        print('     members needing W > 16384 (index, t):', hard)
        print('     first members:', [(m[1], m[4], m[3]) for m in mem[:4]], '(t, D, merge time)')

        print('(F) random giant odd exponents: cheapest exit depth (W = 4096, D <= 64)')
        rr = pool.map(job_random, range(600), chunksize=10)
        for Wc in (128, 256, 512, 1024, 2048, 4096):
            print(f'     P(no exit by {Wc}) = {sum(1 for v in rr if v is None or v > Wc) / len(rr):.3f}')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
