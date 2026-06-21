#!/usr/bin/env python3
"""
lrc14_routeA_tournament_functional_probe_opus_0620.py   (opus-2026-06-20)

ROUTE A PROBE for HYP-2693/2694 (consec maximizes U4): is U4(E) a TOURNAMENT-CYCLE
FUNCTIONAL of the difference-winding tournament T(x), so THM-554/555 (score partition
function) close the compression-extremality?

SETUP (canon): E subset Z, 0 in E, |E|=k. Phases {frac(e x)} on R/Z. 7 sectors [j/7,(j+1)/7).
  N(x) = #{ j in 1..6 : sector j EMPTY of phases }.   p_t = meas{x: N(x)=t}.
  U4(E) = p_0 + p_5 + 5 p_6 = 1 - S_1 + S_2 - S_3 + S_4   (THM-556, re-VERIFIED here).
  T(x): i->j iff frac((e_i-e_j) x) in (0,1/2)  (round tournament, HYP-2605).
  R4: c3(T(x)) = C(k,3) - sum_i C(s_i,2).

RESULT (this script, all EXACT Fractions):
  * THM-556 identity re-VERIFIED: 1-S1+S2-S3+S4 == p0+p5+5p6 (k=8: 353/735).
    [NOTE: differs from a prompt-stated 2633/7350; that is a DIFFERENT g-dual L_y, not U4.]
  * TEST 1  N(x) is NOT determined by the score multiset of T(x)  (nor by (scores,c3)).
  * TEST 3  N(x) is NOT determined by the sorted GAP SPECTRUM (= the cut-space data).
  * TEST 4  N(x) is NOT determined even by the finer SCALE-1/7 SEPTILE DIGRAPH
            (out-set = clockwise 1/7-arc, the R3 object), labeled by phase order.
  * MECHANISM  N requires WALL-ALIGNED empty 1/7-arcs [j/7,(j+1)/7). Translating all
    phases by < 1/7 leaves every pairwise <1/2 (and <1/7) relation -- hence T(x) and the
    septile digraph -- UNCHANGED, yet slides the {j/7} grid across phases and CHANGES N.
    This translation is exactly the symmetry T(x) quotients out. So U4 lives at a STRICTLY
    FINER (continuous, wall-alignment) scale than any discrete tournament invariant.
  * SALVAGE / PARTIAL  consec IS U4-maximal over the full spread<=11 bank (k=8, 330 shapes).
    BUT consec does NOT dominate the tournament-lonely-measure G_c=meas{maxgap>c}: SIX
    two-block shapes (e.g. [0,1,2,3,7,8,9,10]) have LARGER G_{1/2} (more Condorcet-winner /
    "more hierarchical" tournaments) yet U4 ~0.20 << consec 0.48. So the U4-extremum is a
    SCALE-1/7 ALIGNMENT effect, NOT a cut-space (score/maxgap) effect.

CONCLUSION: Route A as 'U4 = closed-form in tournament score/cycle moments' is a DEAD END.
The extremality does not live in the tournament's cut-space (scores/maxgap) nor in its
cycle-space (c3); it lives in the wall-alignment (continuous, scale-1/7) degree of freedom
that T(x) is blind to. Any proof of consec-extremality must keep the 7-grid alignment data.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

# ---------------- core geometry ----------------
def phases(E, x):
    return sorted(set((e * x) % 1 for e in E))

def N_at(E, x):
    P = phases(E, x); cnt = 0
    for j in range(1, 7):
        if not any(F(j, 7) <= p < F(j + 1, 7) for p in P): cnt += 1
    return cnt

def maxgap(E, x):
    P = phases(E, x); g = []
    for i in range(len(P)):
        nxt = P[(i + 1) % len(P)] + (F(1) if i + 1 == len(P) else F(0))
        g.append(nxt - P[i])
    return max(g)

def round_tournament(E, x):
    n = len(E); A = [[0]*n for _ in range(n)]; tf = True
    for i in range(n):
        for j in range(n):
            if i == j: continue
            rel = ((E[i]-E[j]) * x) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2):   A[i][j] = 0
            else:                 A[i][j] = 1 if E[i] < E[j] else 0; tf = False
    return A, tf

def scores_ms(A): return tuple(sorted(sum(r) for r in A))
def c3_count(A):
    n = len(A); c = 0
    for a, b, cc in itertools.combinations(range(n), 3):
        if (A[a][b] and A[b][cc] and A[cc][a]) or (A[a][cc] and A[cc][b] and A[b][a]): c += 1
    return c

def septile_digraph(E, x):
    P = phases(E, x); n = len(P); A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if F(0) < (P[j]-P[i]) % 1 < F(1, 7): A[i][j] = 1
    return tuple(tuple(r) for r in A)

# ---------------- exact measures ----------------
def avoid_measure(E, Aset):
    E = sorted(set(E)); forb = [(F(j, 7), F(j+1, 7)) for j in Aset]
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for (lo, hi) in forb:
            for t in (lo, hi):
                for m in range(e): bps.add((t + m) / e)
    bps = sorted(z for z in bps if 0 <= z <= 1); tot = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; ok = True
        for e in E:
            p = (e * xm) % 1
            for (lo, hi) in forb:
                if lo <= p < hi: ok = False; break
            if not ok: break
        if ok: tot += x1 - x0
    return tot

def S_r(E, r):
    if r == 0: return F(1)
    return sum(avoid_measure(E, A) for A in itertools.combinations(range(1, 7), r))

def dist_p(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1): bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1); p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        p[N_at(E, (lo + hi) / 2)] += (hi - lo)
    return p

def U4(E):
    p = dist_p(E); return p[0] + p[5] + 5*p[6], p

def bps_full(E):
    """breakpoints capturing BOTH N (a/7e) and maxgap (a/14d) transitions."""
    Es = sorted(set(E)); bps = set([F(0), F(1)])
    for e in Es:
        if e == 0: continue
        for a in range(0, 7*e + 1): bps.add(F(a, 7*e))
    for a in range(len(Es)):
        for b in range(a+1, len(Es)):
            d = Es[b] - Es[a]
            if d:
                for m in range(0, 14*d + 1): bps.add(F(m, 14*d))
    return sorted(z for z in bps if 0 <= z <= 1)

def U4_and_maxgapCDF(E, thr=(F(1,2), F(4,7), F(5,7), F(6,7))):
    bps = bps_full(E); p = [F(0)]*7; G = {c: F(0) for c in thr}
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; w = x1 - x0
        p[N_at(E, xm)] += w; mg = maxgap(E, xm)
        for c in thr:
            if mg > c: G[c] += w
    return p[0] + p[5] + 5*p[6], G, p

def breakpoints7(E):
    bps = set([F(0), F(1)])
    for e in sorted(set(E)):
        if e == 0: continue
        for a in range(0, 7*e + 1): bps.add(F(a, 7*e))
    return sorted(x for x in bps if 0 <= x <= 1)

# ============================ DRIVER ============================
if __name__ == "__main__":
    print("="*84)
    print("ROUTE A: is U4 a tournament-score/cycle functional of T(x)?  (EXACT)")
    print("="*84)

    print("\n[THM-556 re-VERIFY]  1-S1+S2-S3+S4  ==  p0+p5+5p6")
    for k in (8, 9, 10):
        E = list(range(k)); _, p = U4(E)
        S1, S2, S3, S4 = S_r(E,1), S_r(E,2), S_r(E,3), S_r(E,4)
        ie = 1 - S1 + S2 - S3 + S4; mom = p[0] + p[5] + 5*p[6]
        print(f"  k={k} consec: IE={ie}  moment={mom}  match={ie==mom}  (={float(mom):.5f})")

    print("\n[TEST 1] N determined by score multiset / (scores,c3) of T(x)?")
    for E in ([0,1,2,3,4], [0,1,2,3,4,5], [0,1,3,7,15]):
        bps = breakpoints7(E); smap = {}; cmap = {}
        for i in range(len(bps)-1):
            x0, x1 = bps[i], bps[i+1]
            if x1 <= x0: continue
            xm = (x0 + x1) / 2; A, tf = round_tournament(E, xm)
            if not tf: continue
            if len(phases(E, xm)) < len(set(E)): continue
            s = scores_ms(A); N = N_at(E, xm)
            smap.setdefault(s, set()).add(N)
            cmap.setdefault((s, c3_count(A)), set()).add(N)
        print(f"  E={E}: score-classes={len(smap)} (>1N:{sum(len(v)>1 for v in smap.values())}); "
              f"(score,c3)-classes={len(cmap)} (>1N:{sum(len(v)>1 for v in cmap.values())})")

    print("\n[TEST 3] N determined by sorted GAP SPECTRUM (cut-space data)?")
    for E in ([0,1,2,3,4], [0,1,2,3,4,5]):
        bps = breakpoints7(E); gmap = {}
        for i in range(len(bps)-1):
            x0, x1 = bps[i], bps[i+1]
            if x1 <= x0: continue
            xm = (x0 + x1) / 2; P = phases(E, xm)
            if len(P) < len(set(E)): continue
            g = []
            for t in range(len(P)):
                nxt = P[(t+1) % len(P)] + (F(1) if t+1 == len(P) else F(0))
                g.append(nxt - P[t])
            gmap.setdefault(tuple(sorted(g, reverse=True)), set()).add(N_at(E, xm))
        print(f"  E={E}: gap-spectra={len(gmap)} (>1N:{sum(len(v)>1 for v in gmap.values())})")

    print("\n[TEST 4] N determined by the SCALE-1/7 SEPTILE DIGRAPH (phase-ordered)?")
    for E in ([0,1,2,3,4], [0,1,3,7,15]):
        bps = breakpoints7(E); dmap = {}
        for i in range(len(bps)-1):
            x0, x1 = bps[i], bps[i+1]
            if x1 <= x0: continue
            xm = (x0 + x1) / 2
            if len(phases(E, xm)) < len(set(E)): continue
            dmap.setdefault(septile_digraph(E, xm), set()).add(N_at(E, xm))
        print(f"  E={E}: septile-digraphs={len(dmap)} (>1N:{sum(len(v)>1 for v in dmap.values())})")

    print("\n[SALVAGE] consec U4-extremal? AND does it dominate the maxgap CDF (cut-space)?")
    k = 8; cE = list(range(k)); cU4, cG, _ = U4_and_maxgapCDF(cE)
    beats = 0; checked = 0; exceptions = []
    for rest in itertools.combinations(range(1, 12), k-1):
        E = [0] + list(rest); u, G, _ = U4_and_maxgapCDF(E); checked += 1
        if u > cU4 + F(1, 10**9): beats += 1
        if any(G[c] > cG[c] + F(1, 10**9) for c in cG):
            exceptions.append((E, u, [str(c) for c in cG if G[c] > cG[c]+F(1,10**9)]))
    print(f"  k=8 bank spread<=11: checked={checked}, consec U4={float(cU4):.5f}, beats={beats}")
    print(f"  shapes with maxgap-survival EXCEEDING consec at some threshold: {len(exceptions)}")
    for E, u, thr in exceptions:
        print(f"    {E}: U4={float(u):.5f} (<<consec) but lonely-meas exceeds consec at {thr}")
    print("\n  => consec maximizes U4 but NOT the tournament lonely-measure: the U4-extremum")
    print("     is a SCALE-1/7 WALL-ALIGNMENT effect, NOT a tournament cut/cycle-space effect.")
    print("\nDONE. Route A (U4 = tournament score/cycle moment functional) = DEAD END.")
