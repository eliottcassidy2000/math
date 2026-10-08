"""The t = 23 fan escapes its barrier (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Companion: mersenne_line_barriers_20261007.py.

Fast Terras parity vectors by divide and conquer: for x known mod 2^n, par(x, n) = (V, a, c) with V the parity vector,
a its weight, and T^n(x') = (3^a x' + c)/2^n for all x' = x mod 2^n; halves combine by x_mid = (3^a1 x + c1) >> n1,
a = a1 + a2, c = 3^a2 c1 + 2^n1 c2.  O(M(n) log n) with python-flint (FLINT fmpz).
Checks:
  (P) par agrees with naive iteration (random x, n <= 5000; the fan bottom's source, n = 2^17).
  (G) the fan bottom M_99708993677: source-reference chains for D = 1..256 with collapse of identical states, W = 2^25:
      all 256 chains collapse into one by step 2^22 and that chain is absorbed at Terras time 19,000,765
      (so M_99708993677 ~> M_(99708993677 - D) for every D <= 256).
  (H) independent residue check from affine maps (not the chain): T^n(x) = T^n(y_D) mod 2^(2^20) at n = 19,000,765 for
      D = 1 and D = 256, with odd-step counts differing by exactly D; the same residues differ at n - 1.
  (I) the contiguous edge of the absorbed block: the chains for D = 1500, 1910 collapse into the D = 1 chain before
      19,000,765; the chain for D = 1911 does not (it stays a separate chain through 19,000,765).  The absorbed group is not
      an interval: the full block D <= 4000 is run in mersenne_fan_block_20261007.py (members 1..1910, 1929..1932, 1935, 1936).
Prints ALL CHECKS PASSED.  Runtime about 2 minutes (single process plus 3 workers); memory about 2 GB."""
import sys, time, random
import numpy as np
from multiprocessing import Pool
from flint import fmpz

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
BASE = 64
E0 = 99708993677
NABS = 19000765


def par_small(x, n):
    V = 0; a = 0; c = 0; p2 = 1; z = x
    for i in range(n):
        if z & 1:
            V |= 1 << i; z = (3 * z + 1) >> 1; c = 3 * c + p2; a += 1
        else:
            z >>= 1
        p2 <<= 1
    return V, a, c


def par(x, n):
    if n <= BASE:
        return par_small(int(x) % (1 << n), n)
    n1 = n // 2; n2 = n - n1
    xf = fmpz(x) & ((fmpz(1) << n) - 1)
    V1, a1, c1 = par(xf & ((fmpz(1) << n1) - 1), n1)
    xm = ((fmpz(3) ** a1 * xf + fmpz(c1)) >> n1) & ((fmpz(1) << n2) - 1)
    V2, a2, c2 = par(xm, n2)
    return int(fmpz(V1) | (fmpz(V2) << n1)), a1 + a2, int(fmpz(3) ** a2 * fmpz(c1) + (fmpz(c2) << n1))


def pow3_mod2(e, W):
    mask = (fmpz(1) << W) - 1
    r = fmpz(1); b = fmpz(3)
    while e:
        if e & 1:
            r = (r * b) & mask
        e >>= 1
        if e:
            b = (b * b) & mask
    return r


def xres(E, W):
    return (2 * pow3_mod2(E - 1, W) - 1) & ((fmpz(1) << W) - 1)


def coins_of(E, nsteps):
    V, a, c = par(xres(E, nsteps + 8), nsteps)
    nb = (nsteps + 7) // 8
    return np.unpackbits(np.frombuffer(int(V).to_bytes(nb, 'little'), dtype=np.uint8), bitorder='little')[:nsteps]


def stepf(k, av, b):
    if av & 1 == 0:
        if b == 0:
            av >>= 1
        elif k >= 0:
            av = (3 * av + 1 - P3[k]) >> 1
        else:
            av = (3 * av + P3[-k] - 1) >> 1
    else:
        if b == 0:
            av = ((3 * av + 1) >> 1) if k >= 0 else ((av + P3[-k - 1]) >> 1); k += 1
        else:
            av = ((av - P3[k - 1]) >> 1) if k >= 1 else ((3 * av - 1) >> 1); k -= 1
    return k, av


COINS = None


def init_worker():
    global COINS
    COINS = coins_of(E0, NABS + 8).tolist()


def lockstep(D):
    k1, a1 = -1, 1 - P3[1]
    k2, a2 = -D, 1 - P3[D]
    for n in range(NABS + 1):
        if (k1, a1) == (k2, a2):
            return D, 'collapsed', n
        if k2 == 0 and a2 == 0:
            return D, 'absorbed separately', n
        b = COINS[n]
        k1, a1 = stepf(k1, a1, b)
        k2, a2 = stepf(k2, a2, b)
    return D, 'separate through the absorption time', NABS


if __name__ == '__main__':
    t0 = time.time()
    print('(P) fast parity vectors')
    rng = random.Random(3)
    okp = True
    for n in (100, 1000, 5000):
        x = rng.getrandbits(n) | 1
        okp &= par(x, n) == par_small(x, n)
    W = (1 << 17) + 8
    z = (2 * pow(3, E0 - 1, 1 << W) - 1) % (1 << W)
    Vn, _, _ = par(z, W - 8)
    bits = bytearray(W - 8); zz = z
    for i in range(W - 8):
        b = zz & 1
        bits[i] = b
        zz = ((3 * zz + 1) >> 1) if b else (zz >> 1)
    Vnaive = int.from_bytes(np.packbits(np.frombuffer(bytes(bits), dtype=np.uint8), bitorder='little').tobytes(), 'little')
    okc = Vnaive == int(Vn)
    check(okp and okc, 'divide-and-conquer parity vectors equal naive iteration (random x, n <= 5000; the fan bottom, n = 2^17 = 131072)')

    print('(G) the fan bottom, chains D = 1..256 to W = 2^25 with state collapse')
    LOGW = 25
    cb = coins_of(E0, (1 << LOGW) - 8)
    print(f'     coins computed ({len(cb)} steps, odd density {cb.mean():.5f}) [{time.time() - t0:.0f}s]', flush=True)
    alive = {D: (-D, 1 - P3[D]) for D in range(1, 257)}
    groups = {D: [D] for D in alive}
    absorbed = None; collapse_one = None; kmin = 0
    n = 0; N = len(cb)
    while alive and n < N:
        blk = cb[n:n + 65536].tolist()
        for b in blk:
            new = {}
            for D, (k, av) in alive.items():
                if k == 0 and av == 0:
                    absorbed = (n, sorted(groups[D]))
                    break
                new[D] = stepf(k, av, b)
            if absorbed:
                break
            seen = {}; alive = {}
            for D, st in new.items():
                if st in seen:
                    groups[seen[st]].extend(groups.pop(D))
                else:
                    seen[st] = D; alive[D] = st
            n += 1
            if collapse_one is None and len(alive) == 1:
                collapse_one = n
            if collapse_one is not None and len(alive) == 1:
                kk = next(iter(alive.values()))[0]
                if kk < kmin:
                    kmin = kk
        if absorbed:
            break
    check(absorbed is not None and absorbed[0] == NABS and absorbed[1] == list(range(1, 257)) and collapse_one is not None and collapse_one <= (1 << 22),
          f'all 256 chains collapse into one by step {collapse_one} (<= 2^22) and it is absorbed at Terras time {absorbed[0] if absorbed else None}: '
          f'M_99708993677 ~> M_(99708993677 - D) for every D <= 256; minimum level of the collapsed chain {kmin}')

    print('(H) independent residue check from affine maps')
    m = 1 << 20

    def Tn_mod(E, n):
        Wt = n + m
        zt = xres(E, Wt)
        Vt, at, ct = par(zt & ((fmpz(1) << n) - 1), n)
        return ((fmpz(3) ** at * zt + fmpz(ct)) >> n) & ((fmpz(1) << m) - 1), at
    rx, ax = Tn_mod(E0, NABS)
    okh = True; cnts = []
    for D in (1, 256):
        ry, ay = Tn_mod(E0 - D, NABS)
        okh &= (rx == ry) and (ay - ax == D)
        cnts.append(ay - ax)
    rx1, _ = Tn_mod(E0, NABS - 1)
    ry1, _ = Tn_mod(E0 - 1, NABS - 1)
    check(okh and rx1 != ry1, f'T^n(x) = T^n(y_D) mod 2^(2^20) at n = {NABS} for D = 1, 256 (odd-step count differences {cnts}); the residues differ at n - 1')

    print('(I) the edge of the absorbed block')
    with Pool(3, initializer=init_worker) as pool:
        out = dict((D, (st, nn)) for D, st, nn in pool.map(lockstep, [1500, 1910, 1911]))
    print('     ', out)
    check(out[1500][0] == 'collapsed' and out[1910][0] == 'collapsed' and out[1911][0] != 'collapsed',
          'D = 1500 and 1910 collapse into the absorbed chain before 19,000,765; D = 1911 stays separate: '
          'M_99708993677 ~> M_99708991767 (D = 1910) is certified, M_99708991766 is not reached by this event')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
