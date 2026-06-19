#!/usr/bin/env python3
"""
lrc14_tournament_reframe_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ANGLE E — TOURNAMENT / PARITY HOME-TURF REFRAME of the LRC(14) residual.

SET-UP.  The cluster offset set E (0 in E, |E|=k) acts on the circle R/Z.  At each
phase x in [0,1), the k points p_i(x) = frac(e_i x) sit on the circle.  The
DIFFERENCE-WINDING map (HYP-2576) makes a tournament T(x) on the k cluster members:
    i -> j  iff  frac((e_i - e_j) x) in (0, 1/2).
HYP-2576 PROVED: when tie-free this is a LOCAL/ROUND tournament — i.e. its vertices
sit in CYCLIC ORDER on the circle and each beats the next (k-1)/2 ... actually the
out-set of each vertex is the half-open semicircle CLOCKWISE from it.  A round
tournament's score sequence is determined by the GAPS between consecutive points.

THE DICTIONARY we test here (small k, exact Fraction):
  (D1) T(x) is round for a.e. x  (re-confirm, and give the explicit semicircle rule).
  (D2) The 1/7-GAP event {maxgap(x) > 1/7}  <->  a tournament-DEGREE event.
       In a round tournament the score s_i = #points in the clockwise semicircle from p_i.
       A circular gap of length g between consecutive points p_i (clockwise-next p_j)
       means: NO point in the arc (p_i, p_j).  We test: is "max circular gap > 1/7"
       equivalent to a statement about the SCORE SEQUENCE / a forbidden sub-pattern?
  (D3) meas(S7) (seven-sector cover) and mu_{1/7} as AVERAGES over x of a 0/1 tournament
       statistic.  Compute E_x[ phi(T(x)) ] for several tournament functionals phi and
       see which one reproduces meas(S7) / measN / mu_{1/7}.
  (D4) THE ODD-CYCLE / CONDORCET ANGLE.  A round tournament has a Hamiltonian cycle
       (it is strongly connected iff every open semicircle is nonempty <=> maxgap < 1/2).
       Its number of directed 3-cycles c3(T) and Hamiltonian-path count H(T) are
       gap-determined.  We compute E_x[c3(T(x))], E_x[H(T(x))] and relate to the
       singular-series / sector quantities.  Key Redei fact: H(T) is ALWAYS ODD.
  (D5) AP-EXTREMALITY <-> transitive tournament?  When E is an AP {0,1,..,k-1}, the
       points frac(e_i x) = frac(i*(x)) are EQUALLY SPACED in winding -> the round
       tournament is as TRANSITIVE-like as possible (a "rotation" of the canonical
       round tournament).  Test whether AP gives the EXTREME tournament statistic.

All exact.  Outputs the dictionary table.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

# ---------- circle / tournament primitives ----------
def phases(E, x):
    """clockwise positions frac(e x) on the circle, with the member index."""
    return [((e * x) % 1, idx) for idx, e in enumerate(E)]

def round_tournament(E, x):
    """Build T(x): i->j iff frac((e_i-e_j)x) in (0,1/2). Return adjacency + tiefree."""
    n = len(E); A = [[0]*n for _ in range(n)]; tiefree = True
    for i in range(n):
        for j in range(n):
            if i == j: continue
            rel = ((E[i]-E[j]) * x) % 1
            if 0 < rel < F(1,2): A[i][j] = 1
            elif rel > F(1,2):   A[i][j] = 0
            else:                A[i][j] = 1 if E[i] < E[j] else 0; tiefree = False
    return A, tiefree

def scores(A):
    return [sum(row) for row in A]

def c3(A):
    n = len(A); c = 0
    for a,b,cc in itertools.combinations(range(n),3):
        # directed 3-cycle either orientation
        if (A[a][b] and A[b][cc] and A[cc][a]) or (A[a][cc] and A[cc][b] and A[b][a]):
            c += 1
    return c

def Hpaths(A):
    """number of directed Hamiltonian PATHS (Redei: always odd). n! perms — small n only."""
    n = len(A); h = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): h += 1
    return h

def is_local(A):
    """round/local: in- and out-neighbourhoods are each transitive."""
    n = len(A)
    for v in range(n):
        for nb in ([u for u in range(n) if A[v][u]], [u for u in range(n) if A[u][v]]):
            for a,b,cc in itertools.permutations(nb,3) if len(nb)>=3 else []:
                if A[a][b] and A[b][cc] and A[cc][a]: return False
    return True

# ---------- circle gap statistics ----------
def maxgap(E, x):
    """largest circular gap between consecutive cluster phases."""
    ps = sorted(set((e*x) % 1 for e in E))
    if len(ps) == 1: return F(1)
    g = F(0)
    for i in range(len(ps)):
        nxt = ps[(i+1) % len(ps)] + (F(1) if i+1==len(ps) else F(0))
        g = max(g, nxt - ps[i])
    return g

def sectors_hit(E, x):
    """which of the 7 fixed sectors [j/7,(j+1)/7) are hit."""
    return set(int(((e*x)%1)*7) for e in E)

# ---------- exact piecewise integration over x ----------
def breakpoints(E):
    """all x where any frac(e_i x), any difference winding, or sector membership changes.
    Use the finest mesh: m/(L) for L = lcm-ish; we just take m/(7*e) and m/d for diffs."""
    bps = set([F(0), F(1)])
    Es = sorted(set(E))
    for e in Es:
        if e == 0: continue
        for m in range(0, 7*e+1): bps.add(F(m, 7*e))   # sector edges
        for m in range(0, 2*e+1): bps.add(F(m, 2*e))   # semicircle / phase=1/2 edges
    diffs = set()
    for a in range(len(Es)):
        for b in range(a+1, len(Es)):
            diffs.add(Es[b]-Es[a]); diffs.add(Es[a]+Es[b])
    for d in diffs:
        if d == 0: continue
        for m in range(0, 2*d+1): bps.add(F(m, 2*d))   # difference-winding crossings of 0,1/2
    return sorted(x for x in bps if 0 <= x <= 1)

def integrate_dictionary(E):
    """Return exact measures of several events + averages of tournament functionals."""
    k = len(E)
    bps = breakpoints(E)
    meas = {"S7": F(0), "N(maxgap<=1/7)": F(0), "nonlocal": F(0),
            "strong(maxgap<1/2)": F(0)}
    avg = {"c3": F(0), "H": F(0), "score_var": F(0), "maxscore": F(0)}
    # also: average of indicator-that-some-score==k-1 (a "near-transitive" / Condorcet-ish event)
    meas["has_full_dominator(score=k-1)"] = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        w = x1 - x0
        xm = (x0 + x1) / 2
        A, tf = round_tournament(E, xm)
        s = scores(A)
        mg = maxgap(E, xm)
        sh = sectors_hit(E, xm)
        if len(sh) == 7: meas["S7"] += w
        if mg <= F(1,7): meas["N(maxgap<=1/7)"] += w
        if not is_local(A): meas["nonlocal"] += w
        if mg < F(1,2): meas["strong(maxgap<1/2)"] += w
        if (k-1) in s: meas["has_full_dominator(score=k-1)"] += w
        # functional averages (weighted by w)
        avg["c3"] += w * c3(A)
        if k <= 8:  # H is n! — keep small
            avg["H"] += w * Hpaths(A)
        mean = F(sum(s), k)
        avg["score_var"] += w * sum((F(si)-mean)**2 for si in s) / k
        avg["maxscore"] += w * max(s)
    return meas, avg

# ===================================================================
print("="*92)
print("ANGLE E: TOURNAMENT REFRAME of LRC(14) residual — exact small-k dictionary")
print("="*92)

# --- D1/D2: the round-tournament <-> max-gap dictionary, exact, small k ---
print("\n[D1+D2] round-tournament round-ness, and max circular gap vs SCORE structure.")
print("Fact to test: in a tie-free round tournament on k points, the maxgap(x) (largest")
print("empty circular arc) > 1/7  <=>  some consecutive pair has an empty 1/7-arc between")
print("them.  We check: maxgap>1/7  <=>  NOT every sector hit?  (S7 vs N relationship).")
print("-"*92)

cap = {8: F(2243,5880), 9: None, 10: None, 11: None, 12: None, 13: None}

# small primitive shapes per k for the dictionary
test_shapes = {
  5:  [("consec", [0,1,2,3,4]),
       ("AP-step2", [0,2,4,6,8]),
       ("perforated", [0,1,2,3,5]),
       ("dissoc", [0,1,3,7,15])],
  6:  [("consec", [0,1,2,3,4,5]),
       ("perforated", [0,1,2,3,4,6]),
       ("dissoc", [0,1,3,7,15,31])],
  7:  [("consec", list(range(7))),
       ("perforated", [0,1,2,3,4,5,7]),
       ("dissoc", [0,1,3,7,15,31,63])],
  8:  [("consec", list(range(8))),
       ("perforated", [0,1,2,3,4,5,6,8]),
       ("dissoc", [0,1,3,7,15,31,63,127])],
}

for k in (5,6,7,8):
    print(f"\n=== k = {k} ===")
    print(f"{'shape':<14}{'meas(S7)':>11}{'measN':>11}{'mu_1/7':>11}{'nonloc':>9}"
          f"{'E[c3]':>9}{'E[H]':>9}{'E[maxsc]':>10}{'P[dom]':>9}")
    for name, E in test_shapes[k]:
        meas, avg = integrate_dictionary(E)
        muN = F(1) - meas["N(maxgap<=1/7)"]
        print(f"{name:<14}{float(meas['S7']):>11.5f}"
              f"{float(meas['N(maxgap<=1/7)']):>11.5f}"
              f"{float(muN):>11.5f}"
              f"{float(meas['nonlocal']):>9.5f}"
              f"{float(avg['c3']):>9.4f}"
              f"{float(avg['H']):>9.4f}"
              f"{float(avg['maxscore']):>10.4f}"
              f"{float(meas['has_full_dominator(score=k-1)']):>9.5f}")

print("\n" + "="*92)
print("DICTIONARY CLAIMS UNDER TEST:")
print(" (C1) nonloc == 0 for all shapes  => T(x) is a.e. a ROUND tournament  [HYP-2576].")
print(" (C2) measN (1/7-net) = P[no empty 1/7-arc] = P[T(x) strongly-connected-ish + dense].")
print(" (C3) mu_1/7 = P[maxgap>1/7] = P[ some vertex's clockwise gap to its successor >1/7 ].")
print("      In round-tournament terms: the consecutive-on-circle pair (v, succ(v)) with the")
print("      empty arc is a 'NEAR-DOMINANCE' edge.  mu_1/7 = P[ >=1 such weak edge ].")
print(" (C4) AP {0..k-1}: does it EXTREMISE E[c3] / E[maxscore] / measN among shapes?")
print("="*92)
print("DONE.")
