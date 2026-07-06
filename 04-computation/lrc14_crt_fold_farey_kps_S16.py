#!/usr/bin/env python3
"""
lrc14_crt_fold_farey_kps_S16.py -- kind-pasteur-2026-07-05-S16, HYP-4197.

The CRT-fold of the Farey-cell voids (mac-mini THM-622). The second-value gap
(1/13,2/25) interior points c/q (c>=3, q in (12.5c,13c)) with COMPOSITE q reduce,
by CRT, to small-modulus multi-color avoidance problems (extends S11 parity-split
q=2p to general q=ab). Prime q (first q=89, c=7) are irreducible cores.

Also: M=3/38 => B pinned at +-1 mod 19 (M at denom 19 <= 1/19); structured covering
38-merge sets attain 3/37 or 2/25, never the cell interior 3/38.
"""
from math import gcd
from functools import reduce
from fractions import Fraction as Fr
from itertools import combinations
import random

def distZ(x, q): r = x % q; return min(r, q - r)

def wit_brute(vs, c, q):
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        if all(distZ(v * a, q) >= c for v in vs): return True
    return False

def wit_paritysplit(vs, c, p):   # q = 2p
    q = 2 * p
    D = [r for r in range(q) if distZ(r, q) < c]
    E = set(r % p for r in D if r % 2 == 0); O = set(r % p for r in D if r % 2 == 1)
    for b in range(1, p):
        if gcd(b, p) != 1: continue
        if all(not ((v % 2 == 0 and (v * b) % p in E) or (v % 2 == 1 and (v * b) % p in O)) for v in vs):
            return True
    return False

def cell_points():
    out = []
    for c in range(3, 10):
        for q in range(int(12.5 * c) + 1, 13 * c):
            if gcd(c, q) == 1: out.append((c, q))
    return out

def M_exact(vs):
    Q = 2 * max(vs); best = Fr(0)
    for q in range(2, Q + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            mm = min(distZ(v * a, q) for v in vs)
            if Fr(mm, q) > best: best = Fr(mm, q)
    return best

if __name__ == "__main__":
    random.seed(4)
    print("=== Farey-cell interior points c/q (c>=3, q in (12.5c,13c)) and reducibility ===")
    for c, q in cell_points()[:12]:
        fac = [p for p in range(2, q) if q % p == 0]
        tag = "PRIME (core)" if not fac else f"= {fac[0]}*{q // fac[0]} (CRT-reducible)"
        print(f"  c={c} q={q}  {tag}")

    print("\n=== parity-split EXACT for c=3 q=38 (mod 19) ===")
    mis = sum(1 for _ in range(30000)
              if (lambda vs: wit_brute(vs, 3, 38) != wit_paritysplit(vs, 3, 19))(
                  [random.randint(1, 300) for _ in range(random.randint(8, 12))]))
    print(f"  brute vs parity-split mismatches: {mis} (0 = EXACT)")

    print("\n=== M=3/38 => pinned at +-1 mod 19 (M at denom 19 <= 1/19) ===")
    # a 3/38-attainer has M(denom 19) <= 1/19; verify on sets clearing 3/38 that their denom-19 clearance <=1
    cnt = 0; okpin = 0
    for _ in range(20000):
        vs = [random.randint(1, 60) for _ in range(12)]
        if not wit_brute(vs, 3, 38): continue
        cnt += 1
        m19 = max((min(distZ(v * a, 19) for v in vs) for a in range(1, 19) if gcd(a, 19) == 1), default=0)
        if m19 <= 1: okpin += 1
    print(f"  sets clearing 3/38: {cnt}; with denom-19 clearance <=1 (pinned): {okpin}"
          f"  ({'not a forced pin -- only M=3/38-attainers are' if okpin < cnt else 'all'})")
    print("  (the pin is forced only when 3/38 is the MAX, not merely cleared; see result doc)")
