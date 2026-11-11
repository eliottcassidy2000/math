# Auditor numerics-prereg: common helpers (written before reading the session scripts).
# Definitions from THM-5128 / THM-5134 (Setting):
#   psi(0)=0, psi(2h)=psi(h)+h, psi(2h+1)=psi(h)+floor(h/2)+1
#   x(h)=2 psi(h)+1, y(h)=2h+1, g(h)=gcd(x,y)
#   F1: N = x1 x2 + 2^s y1 y2 - [s=0], ordered triples (h1,h2,s); sheet 0: s=0, sheet 1: s>=1
import math
import json
import numpy as np
import numba as nb

MODULI = [9, 27, 81, 25, 125, 49, 45, 75, 63, 225, 15, 21]
LCM = 81 * 125 * 49  # 496125, a multiple of every modulus above


@nb.njit
def psi_table(H):
    psi = np.zeros(H + 1, np.int64)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit
def gcd_nb(a, b):
    while b != 0:
        a, b = b, a % b
    return a


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def c_d(m, d):
    # share of residue vectors (alpha,beta) mod m with gcd(alpha,beta,m)=d
    v = 1.0 / (d * d)
    md = m // d
    p = 2
    t = md
    while t > 1:
        if t % p == 0:
            v *= (1.0 - 1.0 / (p * p))
            while t % p == 0:
                t //= p
        p += 1
    return v


def L_m(m, e):
    g = math.gcd(m, e % m) if e % m != 0 else m
    return sum(d * c_d(m, d) for d in divisors(g))


def fold(histL, typL, diagL, m):
    """histL[sheet, r] counts of N mod LCM; typL[sheet, dL] off-diagonal ordered triples whose
    smaller state has gcd(g, LCM) = dL; diagL likewise for h1 = h2.  Returns for modulus m:
    cls[sheet][c], typ[sheet][d], diag[sheet][d] as python ints."""
    out = {}
    L = histL.shape[1]
    idx = np.arange(L) % m
    cls = []
    for sh in range(2):
        acc = np.zeros(m, dtype=np.int64)
        np.add.at(acc, idx, histL[sh])
        cls.append([int(v) for v in acc])
    typ = []
    diag = []
    for sh in range(2):
        t = {d: 0 for d in divisors(m)}
        dg = {d: 0 for d in divisors(m)}
        nz = np.nonzero(typL[sh])[0]
        for dL in nz:
            t[math.gcd(int(dL), m)] += int(typL[sh, dL])
        nz = np.nonzero(diagL[sh])[0]
        for dL in nz:
            dg[math.gcd(int(dL), m)] += int(diagL[sh, dL])
        typ.append(t)
        diag.append(dg)
    return cls, typ, diag


def remainder(m, cls, typ, diag, sheet, include_diag):
    """m * max_c | #{N=c} - sum_{d | gcd(m, c - c_sigma)} (d/m) D^(d) | / #T for one sheet."""
    csig = -1 if sheet == 0 else 0
    T = sum(cls[sheet])
    worst = 0.0
    for c in range(m):
        e = (c - csig) % m
        g0 = math.gcd(m, e) if e != 0 else m
        pred = 0.0
        for d in divisors(g0):
            D = typ[sheet][d] + (diag[sheet][d] if include_diag else 0)
            pred += d * D / m
        worst = max(worst, abs(cls[sheet][c] - pred))
    return m * worst / T


def groups(m, cls, typ, diag, sheet, include_diag):
    """per g0: (mean, min, max of m*share over classes with gcd(m,c-c_sigma)=g0, first-order Phi(g0), limit)"""
    csig = -1 if sheet == 0 else 0
    T = sum(cls[sheet])
    res = {}
    for g0 in divisors(m):
        fac = []
        for c in range(m):
            e = (c - csig) % m
            gg = math.gcd(m, e) if e != 0 else m
            if gg == g0:
                fac.append(m * cls[sheet][c] / T)
        phi = 0.0
        for d in divisors(g0):
            D = typ[sheet][d] + (diag[sheet][d] if include_diag else 0)
            phi += d * D / T
        res[g0] = (sum(fac) / len(fac), min(fac), max(fac), phi, L_m(m, g0), len(fac))
    return res
