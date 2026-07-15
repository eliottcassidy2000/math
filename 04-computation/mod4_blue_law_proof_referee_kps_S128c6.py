#!/usr/bin/env python3
"""
mod4_blue_law_proof_referee_kps_S128c6.py
=========================================
kind-pasteur-2026-07-14-S128 (cont.6).  REFEREE for the PROOF of THM-791(a) upgraded:

THE EVEN-n GRID-SYM STRATUM LAW: at even n, every iso class containing a grid-symmetric tiling
has H = 1 (mod 4).  Hence the blue mod-4 law (both ends of a blue line are such classes).

Proof chain being refereed, exhaustively at n = 4, 5, 6:
 (1) OCF (canon, Grinberg-Stanley): H(T) = I(Omega(T), 2) = Sum_k alpha_k 2^k, alpha_k = # size-k
     sets of pairwise vertex-disjoint directed odd cycles.  =>  H = 1 + 2*c_odd (mod 4),
     c_odd = alpha_1 = # directed odd cycles.                       [check on EVERY class]
 (2) grid-sym tiling t  <=>  rho(i) = n+1-i is an ANTI-automorphism of T(t).   [check per tiling]
 (3) sigma(C) = reverse(rho(C)) is an involution on directed odd cycles; at EVEN n rho is
     fixed-point-free => no rho-invariant odd vertex set => sigma fixed-point-free => c_odd EVEN.
     [check: c_odd even for every grid-sym tiling at n=4,6; find odd-c_odd grid-sym at n=5]
 (4) => H = 1 (mod 4) for every grid-sym-containing class at even n; blue lines preserve H mod 4.
     [check: tabulate (type, H mod 4) per class; PB/MX all 1 mod 4 at n=4,6]
 BONUS certificate: at even n, H = 3 (mod 4) forces PURE BLACK.     [tabulate]
"""
import sys
from math import comb
from itertools import permutations, combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)

def run(n):
    print("=" * 88)
    print("n = %d" % n)
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    m = len(tiles)
    tidx = {t: i for i, t in enumerate(tiles)}
    gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
    pairs = [(i, j) for i in range(n) for j in range(n) if i < j]
    pidx = {p: k for k, p in enumerate(pairs)}

    def beats_of(tv):
        B = [[False] * n for _ in range(n)]
        for k in range(2, n + 1):        # base: k beats k-1  (vertices 1..n -> 0-indexed)
            B[k - 1][k - 2] = True
        for (x, y), i in tidx.items():
            if tv[i] == 1:
                B[x - 1][y - 1] = True
            else:
                B[y - 1][x - 1] = True
        return B

    def ham(B):
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1 << n):
            for v in range(n):
                c = dp[S][v]
                if not c or not (S >> v) & 1:
                    continue
                for u in range(n):
                    if not (S >> u) & 1 and B[v][u]:
                        dp[S | 1 << u][u] += c
        return sum(dp[(1 << n) - 1][v] for v in range(n))

    def c_odd(B):
        tot = 0
        for L in range(3, n + 1, 2):
            for S in combinations(range(n), L):
                u = S[0]
                rest = S[1:]
                # count directed Ham cycles of induced tournament on S: paths u -> ... -> v, close v->u
                for perm in permutations(rest):
                    ok = True
                    prev = u
                    for w in perm:
                        if not B[prev][w]:
                            ok = False
                            break
                        prev = w
                    if ok and B[prev][u]:
                        tot += 1
        return tot

    def canon_key(B):
        best = None
        for pm in permutations(range(n)):
            x = 0
            for k, (i, j) in enumerate(pairs):
                pi, pj = pm[i], pm[j]
                b = 1 if (B[i][j] if pi < pj else not B[i][j]) else 0
                # encode arc (min,max) direction: bit=1 iff lower-image beats higher-image
                lo, hi = (pi, pj) if pi < pj else (pj, pi)
                bit = B[i][j] if pm[i] == lo else (not B[i][j])
                x |= (1 if bit else 0) << pidx[(lo, hi)]
            if best is None or x < best:
                best = x
        return best

    rho = [n - 1 - i for i in range(n)]   # 0-indexed rho
    ok1 = ok2 = ok3 = True
    cls_data = {}
    gs_codd_parities = Counter()
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        B = beats_of(tv)
        gs = all(tv[i] == tv[gmap[i]] for i in range(m))
        # (2) gs <=> rho anti-automorphism
        anti = all((B[u][v] == B[rho[v]][rho[u]]) for u in range(n) for v in range(n) if u != v)
        if (gs != anti):
            ok2 = False
            print("  (2) FAIL at tiling %d: gs=%s anti=%s" % (t, gs, anti))
        H = ham(B)
        co = c_odd(B)
        # (1) OCF mod-4 corollary
        if (H - 1 - 2 * co) % 4 != 0:
            ok1 = False
            print("  (1) FAIL: H=%d c_odd=%d at tiling %d" % (H, co, t))
        if gs:
            gs_codd_parities[co % 2] += 1
            if n % 2 == 0 and co % 2 != 0:
                ok3 = False
                print("  (3) FAIL: even n grid-sym tiling %d has ODD c_odd=%d" % (t, co))
        key = canon_key(B)
        op_key = canon_key([[B[v][u] for u in range(n)] for v in range(n)])
        mkey = min(key, op_key)
        d = cls_data.setdefault(mkey, dict(H=H, gs=0, tot=0))
        d["gs"] += gs
        d["tot"] += 1
        assert d["H"] == H
    print("(1) OCF mod-4 corollary H=1+2*c_odd (mod 4): %s (all %d tilings)" % (ok1, 1 << m))
    print("(2) grid-sym <=> rho anti-automorphism: %s" % ok2)
    print("(3) grid-sym c_odd parities {even,odd}: %s%s"
          % (dict(gs_codd_parities), "  [even-n law: %s]" % ok3 if n % 2 == 0 else "  [odd n: both allowed]"))
    tab = Counter()
    for mk, d in cls_data.items():
        typ = "PB" if d["gs"] == d["tot"] else ("PBk" if d["gs"] == 0 else "MX")
        tab[(typ, d["H"] % 4)] += 1
    print("(4) classes by (type, H mod 4): %s" % dict(sorted(tab.items())))
    if n % 2 == 0:
        bad = [(ty, r) for (ty, r), c in tab.items() if ty in ("PB", "MX") and r != 1]
        print("    even-n stratum law (PB/MX => H=1 mod 4): %s" % ("HOLDS" if not bad else "FAILS %s" % bad))
        print("    bonus certificate: H=3 mod 4 classes are all PBk: %s"
              % ("HOLDS" if all(ty == "PBk" for (ty, r) in tab if r == 3) else "FAILS"))

for n in [4, 5, 6]:
    run(n)
print("\nREFEREE COMPLETE")
