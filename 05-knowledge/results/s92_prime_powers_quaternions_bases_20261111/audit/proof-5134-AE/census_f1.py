# census_f1.py  -- auditor proof-5134-AE, independent exact census of the F1 family.
# Written BEFORE reading the session scripts (odd_modulus_census.py etc.).
#
# F1: N = M - [s = 0],  M = x1 x2 + 2^s y1 y2,  x(h) = 2 psi(h) + 1, y(h) = 2h + 1,
# ordered triples (h1, h2, s), h1, h2 >= 0, s >= 0.  Sheet 0: s = 0 (c_0 = -1); sheet 1: s >= 1 (c_1 = 0).
# T_sigma(X) = triples of the sheet with N <= X.
#
# For every modulus m in the list the script records, for each sheet, the exact number of
#   unordered pairs h1 < h2 (per s summed) by (type d = gcd(g(h1), m) of the SMALLER word, class of M mod m)
#   diagonal triples h1 = h2 by (type of the state, class of M mod m).
# All counts are exact integers (int64, no floating point in the census).
#
# usage: python3 census_f1.py LOG2X m1,m2,...  OUT.json
import sys, math, time, json
import numpy as np
import numba as nb


@nb.njit(cache=True)
def build_psi(H):
    psi = np.zeros(H + 1, np.int32)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit(parallel=True, cache=True)
def census(X, psi, mods, tidx, nchunk, maxm, maxt):
    nm = mods.shape[0]
    off = np.zeros((nchunk, nm, 2, maxt, maxm), np.int64)
    dia = np.zeros((nchunk, nm, 2, maxt, maxm), np.int64)
    smax = 0
    while (1 << (smax + 1)) <= X:
        smax += 1
    for ch in nb.prange(nchunk):
        for s in range(0, smax + 1):
            Xp = X
            sh = 1
            if s == 0:
                Xp = X + 1
                sh = 0
            p2 = np.int64(1) << np.int64(s)
            h1 = np.int64(ch)
            while True:
                y1 = np.int64(2 * h1 + 1)
                x1 = np.int64(2 * psi[h1] + 1)
                if p2 * y1 * y1 > Xp:
                    break
                Md = x1 * x1 + p2 * y1 * y1
                if Md <= Xp:
                    for j in range(nm):
                        dia[ch, j, sh, tidx[j, h1], Md % mods[j]] += 1
                den = x1 + 2 * p2 * y1
                y2max = (2 * Xp - 1) // den
                h2max = (y2max - 1) // 2
                A = p2 * y1
                for h2 in range(h1 + 1, h2max + 1):
                    M = x1 * np.int64(2 * psi[h2] + 1) + A * np.int64(2 * h2 + 1)
                    if M <= Xp:
                        for j in range(nm):
                            off[ch, j, sh, tidx[j, h1], M % mods[j]] += 1
                h1 += np.int64(nchunk)
    return off, dia


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def main():
    L = int(sys.argv[1])
    mods = [int(t) for t in sys.argv[2].split(",")]
    outp = sys.argv[3]
    X = 1 << L
    t0 = time.time()
    H = ((2 * (X + 1) - 1) // 3) // 2 + 2
    psi = build_psi(H)
    # sanity of psi against the definition (first values) : psi(1)=1, psi(2)=2, psi(3)=2, psi(4)=4, psi(5)=4
    assert [int(psi[h]) for h in range(8)] == [0, 1, 2, 2, 4, 4, 5, 4], [int(psi[h]) for h in range(8)]
    h1max = int(math.isqrt(X + 1)) // 2 + 2
    nm = len(mods)
    divs = [divisors(m) for m in mods]
    maxt = max(len(d) for d in divs)
    maxm = max(mods)
    tidx = np.zeros((nm, h1max + 1), np.int64)
    for h in range(h1max + 1):
        x = 2 * int(psi[h]) + 1
        y = 2 * h + 1
        g = math.gcd(x, y)
        for j, m in enumerate(mods):
            tidx[j, h] = divs[j].index(math.gcd(g, m))
    nchunk = 240
    off, dia = census(np.int64(X), psi, np.array(mods, np.int64), tidx, np.int64(nchunk), np.int64(maxm), np.int64(maxt))
    off = off.sum(axis=0)
    dia = dia.sum(axis=0)
    res = {"log2X": L, "mods": mods, "data": {}}
    for j, m in enumerate(mods):
        res["data"][str(m)] = {
            "divisors": divs[j],
            "off": [[[int(off[j, sh, t, r]) for r in range(m)] for t in range(len(divs[j]))] for sh in range(2)],
            "diag": [[[int(dia[j, sh, t, r]) for r in range(m)] for t in range(len(divs[j]))] for sh in range(2)],
        }
    with open(outp, "w") as f:
        json.dump(res, f)
    print("X = 2^%d  moduli %s  time %.1f s  -> %s" % (L, mods, time.time() - t0, outp))
    for sh in range(2):
        j = 0
        m = mods[0]
        T = 2 * int(off[j, sh, :, :m].sum()) + int(dia[j, sh, :, :m].sum())
        print("sheet %d: #T = %d  (diagonal %d)" % (sh, T, int(dia[j, sh, :, :m].sum())))


if __name__ == "__main__":
    main()
