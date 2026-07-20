#!/usr/bin/env python3
"""gmc2_pattern_batch_kps_S128c127.py -- kind-pasteur-2026-07-20-S128c127
Batch the finite Groebner emptiness test over ALL genuinely-two-sided charge patterns of
bounded span -> a systematic 'GMC(2) on span <= s' certificate (extends opus THM-1685's
per-pattern TNC to the full GMC moments)."""
import sys
from math import factorial
from itertools import combinations
import sympy as sp

w = sp.Symbol('w')

def mul_sym(p, q):
    out = {}
    for (a, b), u in p.items():
        for (c, d), v in q.items():
            k = (a + c, b + d)
            out[k] = out.get(k, 0) + u * v
    return {k: sp.expand(v) for k, v in out.items() if v != 0}

def Emom_sym(p):
    return sp.expand(sum(v * factorial(a) for (a, b), v in p.items() if a == b))

def gmc_empty(charges, K):
    coeffs = sp.symbols('a0:%d' % len(charges))
    P = {}
    for idx, c in enumerate(charges):
        a, b = (c, 0) if c >= 0 else (0, -c)
        P[(a, b)] = coeffs[idx]
    cmax, cmin = max(charges), min(charges)
    ext = [coeffs[i] for i, c in enumerate(charges) if c in (cmax, cmin)]
    cur = {(0, 0): sp.Integer(1)}
    moms = []
    for m in range(1, K + 1):
        cur = mul_sym(cur, P)
        e = Emom_sym(cur)
        if e != 0:
            moms.append(e)
    gens = moms + [1 - w * sp.prod(ext)]
    G = sp.groebner(gens, *coeffs, w, order='grevlex')
    return list(G) == [sp.Integer(1)]

# enumerate two-sided patterns: min charge = -M, max charge = N, span = M+N <= SMAX
SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
total = ok = 0
fails = []
print("Batch GMC(2) finite-Groebner over all two-sided charge patterns, span <= %d" % SMAX)
print("(both extremes present; interior charges optional)")
print()
for span in range(2, SMAX + 1):
    for M in range(1, span):
        N = span - M
        interior = [c for c in range(-M, N + 1) if c not in (-M, N)]
        pats = 0; pok = 0
        for r in range(len(interior) + 1):
            for sub in combinations(interior, r):
                charges = sorted([-M] + list(sub) + [N])
                K = 2 * span + 4
                total += 1; pats += 1
                try:
                    if gmc_empty(charges, K):
                        ok += 1; pok += 1
                    else:
                        fails.append(charges)
                except Exception as e:
                    fails.append((charges, str(e)[:40]))
        print("  span=%d (M=%d,N=%d): %d/%d patterns EMPTY (GMC(2) holds)" % (span, M, N, pok, pats))
    sys.stdout.flush()
print()
print("TOTAL: %d/%d two-sided patterns closed by finite Groebner (span<=%d)" % (ok, total, SMAX))
print("failures: %s" % (fails if fails else "NONE"))
