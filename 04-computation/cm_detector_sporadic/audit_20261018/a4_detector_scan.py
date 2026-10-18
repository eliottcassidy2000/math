"""Audit checks 4 and 8: the (B1') and (B1'_ord) scans, by residue tables.

(B1') at p (p not dividing m):  sum_{k=p}^{2p-1} u_k m^{-k} == u_1 m^{-1} sum_{k<p} u_k m^{-k}  (mod p^2).
The truth value depends only on m mod p^2, so for each prime we tabulate it over all residues r (p not | r) once,
with x = r^{-1} mod p^2 and numpy int64 Horner evaluation, then index by m mod p^2.  A prime dividing m is skipped.
(B1'_ord): a prime only counts if F_p(1/m) = sum_{k<p} u_k m^{-k} != 0 (mod p).
Survivors: m with no failing prime among 7 <= p <= PMAX.
Density: fraction of primes 7 <= p <= 500, p not | m, with F_p(1/m) == 0 (mod p).

usage: python a4_detector_scan.py M [families...]      (default: all 15 families)
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in, FAMILIES

PMAX = 200
PRIMES = primes_in(7, PMAX)
DPRIMES = primes_in(7, 500)
HERE = os.path.dirname(os.path.abspath(__file__))


def tables(u, p):
    p2 = p * p
    r = np.arange(p2, dtype=np.int64)
    good = (r % p) != 0
    x = np.zeros(p2, dtype=np.int64)
    x[good] = [pow(int(v), -1, p2) for v in r[good]]
    uk = [int(v) % p2 for v in u[:2 * p]]
    # F = sum_{k<p} u_k x^k ; G = sum_{k=p}^{2p-1} u_k x^k   (Horner, mod p^2)
    F = np.zeros(p2, dtype=np.int64)
    for k in range(p - 1, -1, -1):
        F = (F * x + uk[k]) % p2
    G = np.zeros(p2, dtype=np.int64)
    for k in range(2 * p - 1, p - 1, -1):
        G = (G * x + uk[k]) % p2
    xp = np.ones(p2, dtype=np.int64)
    for _ in range(p):
        xp = (xp * x) % p2
    lhs = (G * xp) % p2
    rhs = (int(uk[1]) * x % p2 * F) % p2
    passB1 = (lhs == rhs) | ~good
    Fzero = (F % p == 0) & good
    passOrd = passB1 | Fzero
    return passB1, passOrd


def ftable(u, p):  # F_p(1/m) == 0 mod p, as a table over m mod p
    r = np.arange(p, dtype=np.int64)
    good = r != 0
    x = np.zeros(p, dtype=np.int64)
    x[good] = [pow(int(v), -1, p) for v in r[good]]
    uk = [int(v) % p for v in u[:p]]
    F = np.zeros(p, dtype=np.int64)
    for k in range(p - 1, -1, -1):
        F = (F * x + uk[k]) % p
    return (F == 0) & good, good


def scan(name, M):
    u = seq(name, 500)
    ms = np.concatenate([np.arange(-M, 0), np.arange(1, M + 1)]).astype(np.int64)
    aliveB = np.ones(ms.size, dtype=bool)
    aliveO = np.ones(ms.size, dtype=bool)
    for p in PRIMES:
        pB, pO = tables(u, p)
        idx = ms % (p * p)
        aliveB &= pB[idx]
        aliveO &= pO[idx]
    survB = ms[aliveB].tolist()
    survO = ms[aliveO].tolist()
    # density for every m (supersingular density over 7 <= p <= 500)
    zero_cnt = np.zeros(ms.size, dtype=np.int32)
    tot_cnt = np.zeros(ms.size, dtype=np.int32)
    for p in DPRIMES:
        z, g = ftable(u, p)
        idx = ms % p
        zero_cnt += z[idx]
        tot_cnt += g[idx]
    dens = zero_cnt / np.maximum(tot_cnt, 1)
    hi = ms[dens >= 0.3].tolist()
    keep = set(survB) | set(survO) | set(hi)
    dens_of = {int(m): round(float(dens[i]), 3) for i, m in enumerate(ms.tolist()) if m in keep}
    nonsurv_max = float(np.max(np.where(aliveB | aliveO, 0, dens)))
    return dict(B1=survB, B1ord=survO, dens_ge_0p3=hi, density=dens_of, max_density_nonsurvivors=round(nonsurv_max, 3))


if __name__ == '__main__':
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    fams = sys.argv[2:] or list(FAMILIES)
    out = {}
    for name in fams:
        t0 = time.time()
        res = scan(name, M)
        out[name] = res
        print(f'{name:8s} w={FAMILIES[name][1]}  M={M}: (B1\') survivors {res["B1"]}; (B1\'_ord) {res["B1ord"]}; '
              f'density>=0.3 {res["dens_ge_0p3"]}; max density among non-survivors {res["max_density_nonsurvivors"]}  '
              f'[{time.time() - t0:.1f}s]', flush=True)
        print('          densities:', res['density'], flush=True)
    with open(os.path.join(HERE, f'a4_detector_scan_M{M}.json'), 'w') as f:
        json.dump(out, f, indent=1)
