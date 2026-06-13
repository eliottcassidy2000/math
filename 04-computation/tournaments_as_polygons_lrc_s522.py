#!/usr/bin/env python3
"""
tournaments_as_polygons_lrc_s522.py    oracle-2026-06-01-S522

Seeding/thinking: model tournaments as POLYGONS inscribed in the unit circle and
ask how the REGULAR polygon (roots of unity) sits inside the LRC source menu of
S520/HYP-1987 -- and how that meets the twin-Goldbach mod-6 wheel of S521.

Dictionary
----------
m runner positions on the circle -> m points x_i = e^{2 pi i theta_i}.
HALF-TURN tournament: i -> j iff frac(theta_j - theta_i) in (0, 1/2).
Equally spaced points theta_i = i/m  = REGULAR m-GON = m-th roots of unity.

Claims tested
-------------
(A) The regular m-gon half-turn tournament is the ROTATIONAL (circulant)
    tournament R_m with connection set {1,..,(m-1)/2}  for m ODD; for m EVEN the
    antipodal pairs land at frac=1/2 exactly -> TIE (no arc) -> the regular even-
    gon is DEGENERATE, sitting on a wall of the arrangement (not a tournament).
(B) Compute H (Hamiltonian-path count) + score of R_m and of the Paley
    tournament P_m (m prime, m=3 mod 4), and cross-check against the S520 LRC
    reachable source-menu H-values:
        n=7 (m=6): H in {1,17,23,41,45}
        n=8 (m=7): H in {1,33,47,51,70,105,123,137,151,175}
    i.e. is the regular-polygon tournament a REACHABLE LRC source class?
"""
from itertools import permutations
from functools import lru_cache
from fractions import Fraction

ONE = Fraction(1)
def frac(x): return x - (x.numerator // x.denominator) if isinstance(x, Fraction) else x - int(x)

def half_turn_regular(m):
    """Half-turn tournament of the regular m-gon (theta_i = i/m). Returns
    (adj, ties): adj[i][j]=1 if i->j; ties=list of antipodal tie pairs."""
    adj = [[0]*m for _ in range(m)]
    ties = []
    for i in range(m):
        for j in range(m):
            if i == j: continue
            f = frac(Fraction(j - i, m))
            if 0 < f < Fraction(1,2):
                adj[i][j] = 1
            elif f == Fraction(1,2):
                ties.append((i,j))
    return adj, ties

def rotational(m):
    """R_m: i->j iff (j-i) mod m in {1,..,floor((m-1)/2)} (well-defined, m odd)."""
    half = (m-1)//2
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for k in range(1, half+1):
            adj[i][(i+k)%m] = 1
    return adj

def paley(p):
    """Paley tournament on p vertices (p prime, p=3 mod4): i->j iff j-i is QR."""
    qr = set((x*x) % p for x in range(1,p))
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in qr: adj[i][j]=1
    return adj

def Hc(adj):
    m = len(adj); full=(1<<m)-1
    @lru_cache(None)
    def dp(mask,last):
        if mask==full: return 1
        return sum(dp(mask|(1<<x),x) for x in range(m) if not (mask>>x)&1 and adj[last][x])
    return sum(dp(1<<s,s) for s in range(m))

def scores(adj): return tuple(sorted(sum(r) for r in adj))
def is_tournament(adj):
    m=len(adj)
    return all(adj[i][j]+adj[j][i]==1 for i in range(m) for j in range(m) if i!=j)
def canon(adj):
    m=len(adj); best=None
    for p in permutations(range(m)):
        flat=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or flat<best: best=flat
    return best

MENU = {6:{1,17,23,41,45}, 7:{1,33,47,51,70,105,123,137,151,175}}

def main():
    print("tournaments as polygons on the unit circle — regular polygon vs LRC menu (oracle-S522)\n")
    for m in range(3,8):
        adj,ties = half_turn_regular(m)
        n = m+1
        tour = is_tournament(adj)
        line = f"m={m} (n={n}): regular {m}-gon half-turn -> "
        if not tour:
            line += f"DEGENERATE: {len(ties)//2} antipodal TIE pair(s) at frac=1/2 (a WALL, not a tournament)."
            print(line)
            continue
        H = Hc(adj); sc = scores(adj)
        # compare to rotational R_m
        same_as_R = canon(adj)==canon(rotational(m))
        line += f"tournament H={H}, score={sc}, == R_{m}(rotational): {same_as_R}"
        print(line)
        inmenu = (m in MENU and H in MENU[m])
        print(f"        regular-polygon H={H} in S520 LRC source menu for m={m}? {inmenu}"
              + (f"   <-- REACHABLE LRC SOURCE" if inmenu else ""))
        # Paley, if applicable
        if m in (3,7):
            P = paley(m); HP=Hc(P)
            print(f"        Paley P_{m}: H={HP}, score={scores(P)}, ==R_{m}? {canon(P)==canon(rotational(m))}, in menu? {HP in MENU.get(m,set())}")
    print("\nNOTE: m even => regular polygon is a wall (antipodal ties). m=6 (the hexagon,"
          "\nthe twin-Goldbach mod-6 wheel) is exactly n=7, where the S520 menu COLLAPSES.")

if __name__=="__main__":
    main()
