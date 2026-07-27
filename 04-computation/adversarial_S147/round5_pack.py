#!/usr/bin/env python3
"""Round 5: pack-at-top (k>=7-far duty-stacked, twin-peak) scan.

V = smalls (6 elements <= 14) u pack {N+d_0,...,N+d_6}, N in [65, 6500].
Per modulus q in [14,60], the admissible translations T(q) = {tau : smalls u
{tau+d_k} covers-or-pins q} are precomputed; N must satisfy N mod q in T(q)
for all q, plus leftover covering duties (some pack slot ≡ 0 mod the missing
q <= 13).  Survivors get the fold-anatomy check at the three bands and
exact M.  Kill reasons recorded by band:
   PACKSUM  q = 2N + d_i + d_j   (pack pair sums, the mirror-fold band)
   MIXED    q = N + d_k + s      (pack-small sums; also diffs N + d_k - s)
   DOUBLE   q = 2(N + d_k)
   INTERNAL q <= 120 small moduli (via exact M)
"""
import sys, time
import numpy as np
from math import gcd
from fractions import Fraction
from itertools import combinations
sys.path.insert(0, ".")
from lrc_engine import exact_M, covering_check, d_gap, W_of_q, THR, FLOOR

NMAX = 6500
QLO, QHI = 14, 60

MENUS = [
    (8, 9, 10, 11, 12, 13),   # leftover duty: 7
    (7, 8, 9, 10, 11, 12),    # leftover: 13
    (1, 8, 9, 10, 11, 12),    # leftover: 7, 13
    (1, 9, 10, 11, 12, 13),   # leftover: 7, 8
    (2, 3, 10, 11, 12, 13),   # leftover: 7, 8, 9
    (1, 2, 9, 10, 11, 13),    # leftover: 7, 8, 12(->covered by 12? no 12 absent: leftover 12? 4|8? 4: absent-> 4,7,8,12... recompute in code
    (1, 2, 3, 4, 5, 6),       # leftover: 7,8,9,10,11,12,13 (all to pack)
    (5, 8, 9, 11, 12, 13),    # leftover: 7, 10
    (1, 7, 9, 11, 12, 13),    # leftover: 8, 10? 10=2*5: 5 absent: leftover 5,8,10
    (1, 7, 8, 11, 12, 13),    # leftover: 5,9,10
]

PATTERNS = [
    (0, 1, 2, 3, 4, 5, 6),
    (0, 1, 2, 3, 4, 5, 7),
    (0, 1, 2, 3, 4, 6, 8),
    (0, 1, 2, 3, 5, 7, 9),
    (0, 1, 2, 4, 6, 8, 10),
    (0, 1, 3, 5, 7, 9, 11),
    (0, 2, 4, 6, 8, 10, 12),
    (0, 1, 2, 3, 4, 6, 10),
    (0, 1, 2, 5, 8, 9, 12),
    (0, 1, 3, 4, 9, 10, 12),
    (0, 1, 2, 3, 7, 8, 9),
    (0, 1, 4, 9, 15, 22, 31),
    (0, 2, 3, 4, 5, 6, 8),
    (0, 3, 4, 5, 6, 7, 10),
]

def leftover(smalls):
    return [q for q in range(2, 14) if all(v % q for v in smalls)]

def tau_sets(smalls, pat):
    """T(q) for q in [QLO,QHI]: translations tau making the full set pass q."""
    T = {}
    for q in range(QLO, QHI + 1):
        d = d_gap(q)
        units = [a for a in range(1, q) if gcd(a, q) == 1]
        good = np.zeros(q, dtype=bool)
        for tau in range(q):
            res = [s % q for s in smalls] + [(tau + dd) % q for dd in pat]
            if 0 in res:
                good[tau] = True
                continue
            ok = True
            for a in units:
                m = min(min((a * r) % q, q - (a * r) % q) for r in res)
                if m > d:
                    ok = False
                    break
            good[tau] = ok
        if not good.any():
            return None  # template dead at q
        T[q] = good
    return T

def main():
    t0 = time.time()
    from collections import Counter
    killband = Counter()
    results = []
    n_crt = 0
    for smalls in MENUS:
        need = leftover(smalls)
        for pat in PATTERNS:
            T = tau_sets(smalls, pat)
            if T is None:
                continue
            # CRT admissibility over N
            Ns = np.arange(NMAX + 1)
            adm = np.ones(NMAX + 1, dtype=bool)
            adm[:65] = False
            for q, good in T.items():
                adm &= good[Ns % q]
            # leftover covering: each missing q needs some pack slot ≡ 0
            for qq in need:
                okq = np.zeros(NMAX + 1, dtype=bool)
                for dd in pat:
                    okq |= ((Ns + dd) % qq == 0)
                adm &= okq
            cand = np.nonzero(adm)[0]
            for N in cand:
                N = int(N)
                V = sorted(set(smalls) | {N + dd for dd in pat})
                if len(V) != 13:
                    continue
                g = 0
                for x in V:
                    g = gcd(g, x)
                if g != 1 or covering_check(V):
                    continue
                n_crt += 1
                # fold anatomy, cheap first: packsums, mixed, doubles
                kill = None
                for i, j in combinations(range(7), 2):
                    q = 2 * N + pat[i] + pat[j]
                    W, k = W_of_q(V, q)
                    if W > d_gap(q):
                        kill = ("PACKSUM", q, k, W)
                        break
                if not kill:
                    for dd in pat:
                        for s in smalls:
                            for q in (N + dd + s, N + dd - s):
                                if q >= 55:
                                    W, k = W_of_q(V, q)
                                    if W > d_gap(q):
                                        kill = ("MIXED", q, k, W)
                                        break
                            if kill: break
                        if kill: break
                if not kill:
                    for dd in pat:
                        q = 2 * (N + dd)
                        W, k = W_of_q(V, q)
                        if W > d_gap(q):
                            kill = ("DOUBLE", q, k, W)
                            break
                if kill:
                    killband[kill[0]] += 1
                    if n_crt <= 40 or n_crt % 200 == 0:
                        print(f"  kill {kill[0]} q={kill[1]} W={kill[3]}"
                              f">{d_gap(kill[1])} N={N} pat={pat} "
                              f"smalls={smalls}", flush=True)
                    continue
                M, (q, k) = exact_M(V)
                tag = ("*** IN GAP ***" if FLOOR < M < THR else
                       ("TIGHT" if M == FLOOR else ""))
                results.append((M, tuple(V), q, k))
                print(f"FULL M={M}={float(M):.6f} wit {k}/{q} V={V} {tag}",
                      flush=True)
    print(f"\n[{time.time()-t0:.0f}s] CRT-passing (N,template): {n_crt}")
    print("fold-kill bands:", dict(killband))
    results.sort()
    for M, V, q, k in results[:20]:
        tag = "*** IN GAP ***" if FLOOR < M < THR else ""
        print(f"M={str(M):>9}={float(M):.6f} wit {k}/{q} V={list(V)} {tag}")

if __name__ == "__main__":
    main()
