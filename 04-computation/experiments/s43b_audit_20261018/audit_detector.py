#!/usr/bin/env python3
"""S43b referee: reproduction of the Lucas-block CM detector scan (HYP-9310(c)) by a different method, plus an
independent CM signature over ALL parameters, plus a j-invariant CM test of every survivor.

(B1')  sum_{k=p}^{2p-1} t_k == t_1 sum_{k<p} t_k (mod p^2), t_k = A(k) m^-k.
Method: for each prime p the truth table of (B1') on all units r mod p^2 is computed (vectorized), then every
0 < |m| <= M is looked up (m mod p^2); primes dividing m are skipped.  Survivors are then re-tested on all primes
p <= PCONF by direct computation.
Independent CM signature: s_p(m) = [sum_{k<p} A(k) m^-k == 0 mod p] ("supersingular"), density over 7 <= p <= PSS,
computed for EVERY 0 < |m| <= M.
CM test of survivors: j-invariant from the Hauptmodul relation, tau = j^-1(j) by mpmath, tau quadratic?
Usage: python audit_detector.py [M] [PSCAN] [PCONF] [PSS]
"""
import sys
from math import comb, isqrt
import numpy as np

FAM = {
    'A1=C(2k,k)^3': lambda k: comb(2 * k, k) ** 3,
    'A2=C(2k,k)^2C(3k,k)': lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k),
    'A3=C(2k,k)^2C(4k,2k)': lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k),
    'A4=C(2k,k)C(3k,k)C(6k,3k)': lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k),
}


def primes_in(a, b):
    return [p for p in range(max(a, 2), b + 1) if all(p % d for d in range(2, isqrt(p) + 1))]


def vec_pow(base, e, mod):
    res = np.ones_like(base)
    b = base % mod
    while e:
        if e & 1:
            res = (res * b) % mod
        b = (b * b) % mod
        e >>= 1
    return res


def truth_table(Af, p):
    pp = p * p
    r = np.arange(pp, dtype=np.int64)
    unit = (r % p) != 0
    w = vec_pow(r, p * (p - 1) - 1, pp)          # r^{-1} mod p^2 on units
    A = [Af(k) % pp for k in range(2 * p)]
    S0 = np.zeros(pp, dtype=np.int64)
    S1 = np.zeros(pp, dtype=np.int64)
    pw = np.ones(pp, dtype=np.int64)
    for k in range(2 * p):
        if k < p:
            S0 = (S0 + A[k] * pw) % pp
        else:
            S1 = (S1 + A[k] * pw) % pp
        pw = (pw * w) % pp
    ok = ((S1 - (A[1] * w % pp) * S0) % pp) == 0
    return ok | ~unit


def direct_ok(Af, m, p):
    pp = p * p
    if m % p == 0:
        return True
    w = pow(m % pp, -1, pp)
    S0 = 0
    S1 = 0
    pw = 1
    for k in range(2 * p):
        t = Af(k) % pp * pw % pp
        if k < p:
            S0 += t
        else:
            S1 += t
        pw = pw * w % pp
    return (S1 - Af(1) * w * S0) % pp == 0


def ss_tables(Af, p):
    """boolean over residues w mod p: sum_{k<p} A(k) w^k == 0 (mod p)."""
    w = np.arange(p, dtype=np.int64)
    A = [Af(k) % p for k in range(p)]
    S = np.zeros(p, dtype=np.int64)
    pw = np.ones(p, dtype=np.int64)
    for k in range(p):
        S = (S + A[k] * pw) % p
        pw = (pw * w) % p
    return S == 0


def jinv(fam, m, mp):
    """the two j-invariants attached to integer m (complex mpf)."""
    if fam.startswith('A1'):
        return [mp.mpf((m - 16) ** 3) / m]
    if fam.startswith('A4'):
        return [mp.mpf(m)]
    if fam.startswith('A2'):
        u = mp.mpf(27) / m
        xs = [(1 + s * mp.sqrt(mp.mpc(1 - 4 * u))) / 2 for s in (1, -1)]
        return [27 * (1 + 8 * x) ** 3 / (x * (1 - x) ** 3) for x in xs]
    if fam.startswith('A3'):
        u = mp.mpf(64) / m
        xs = [(1 + s * mp.sqrt(mp.mpc(1 - 4 * u))) / 2 for s in (1, -1)]
        return [64 * (1 + 3 * x) ** 3 / (x * (1 - x) ** 2) for x in xs]


def tau_of_j(J, mp):
    """tau in the upper half plane with j(tau) = J (Newton from a coarse grid of the fundamental domain)."""
    target = mp.mpc(J) / 1728
    best = None
    for xr in [i / 20.0 - 0.5 for i in range(21)]:
        for yi in [0.87, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
            t0 = mp.mpc(xr, yi)
            try:
                d = abs(mp.kleinj(t0) - target)
            except Exception:
                continue
            if best is None or d < best[0]:
                best = (d, t0)
    t = best[1]
    for _ in range(200):
        f = mp.kleinj(t) - target
        h = mp.mpf(10) ** (-mp.mp.dps // 2)
        df = (mp.kleinj(t + h) - mp.kleinj(t - h)) / (2 * h)
        step = f / df
        t = t - step
        if t.imag <= 0:
            t = mp.mpc(t.real, abs(t.imag) + 0.1)
        if abs(step) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
    return t


def quad_disc(t, mp):
    """if tau = u + iv is imaginary quadratic with small height, return its discriminant b^2 - 4ac."""
    for _ in range(100):                      # reduce to the standard fundamental domain
        t = t - mp.nint(t.real)
        if abs(t) < 1 - mp.mpf(10) ** (-40):
            t = -1 / t
        else:
            break
    u = t.real
    n2 = t.real ** 2 + t.imag ** 2
    for a in range(1, 3000):
        b = -2 * a * u
        c = a * n2
        bi = mp.nint(b)
        ci = mp.nint(c)
        if abs(b - bi) < mp.mpf(10) ** (-25) and abs(c - ci) < mp.mpf(10) ** (-25):
            from math import gcd
            A_, B_, C_ = a, int(bi), int(ci)
            g = gcd(gcd(A_, abs(B_)), abs(C_))
            A_, B_, C_ = A_ // g, B_ // g, C_ // g
            return B_ * B_ - 4 * A_ * C_
    return None


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    PSCAN = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    PCONF = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    PSS = int(sys.argv[4]) if len(sys.argv) > 4 else 1200
    ms = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    scan_primes = primes_in(7, PSCAN)
    ss_primes = primes_in(101, PSS)
    import mpmath as mp
    mp.mp.dps = 60
    for name, Af in FAM.items():
        mask = np.ones(len(ms), dtype=bool)
        for p in scan_primes:
            ok = truth_table(Af, p)
            mask &= ok[ms % (p * p)]
        surv = [int(x) for x in ms[mask]]
        conf = {m: all(direct_ok(Af, m, p) for p in primes_in(7, PCONF)) for m in surv}
        # independent supersingular density for ALL m
        cnt = np.zeros(len(ms), dtype=np.int64)
        tot = np.zeros(len(ms), dtype=np.int64)
        for p in ss_primes:
            z = ss_tables(Af, p)
            r = ms % p
            unit = r != 0
            inv = np.zeros(p, dtype=np.int64)
            for x in range(1, p):
                inv[x] = pow(x, -1, p)
            cnt += (unit & z[inv[r]]).astype(np.int64)
            tot += unit.astype(np.int64)
        dens = cnt / np.maximum(tot, 1)
        high = [int(x) for x in ms[dens >= 0.30]]
        nons = dens[~mask]
        print('%s:' % name)
        print('   (B1\') survivors over 0<|m|<=%d, primes 7..%d: %s' % (M, PSCAN, surv))
        print('   re-confirmed on all primes 7..%d: %s' % (PCONF, all(conf.values())))
        print('   supersingular density >= 0.30 over primes 101..%d (independent, all m): %s' % (PSS, high))
        print('   same set as survivors: %s; max density among non-survivors: %.3f; survivor densities: %s'
              % (sorted(high) == sorted(surv), nons.max(), [round(float(dens[np.where(ms == m)[0][0]]), 3) for m in surv]))
        cmres = []
        for m in surv:
            Ds = []
            for J in jinv(name, m, mp):
                t = tau_of_j(J, mp)
                Ds.append(quad_disc(t, mp))
            cmres.append((m, Ds))
        print('   CM test (discriminant of tau with j(tau) = j attached to m): %s' % cmres, flush=True)


if __name__ == '__main__':
    main()
