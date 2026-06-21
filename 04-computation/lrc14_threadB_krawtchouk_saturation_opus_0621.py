#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_krawtchouk_saturation_opus_0621.py  (opus, 2026-06-21, THREAD B)

THE DELSARTE-LP SATURATION of "consec maximizes L_y".  All exact (Fractions).

SETUP (THM-534, HYP-2726):
  N(x) = #missed sectors among {1..6} (sector 0 always hit by e=0).  measS7(E)=P(N=0)=p_0.
  p_t = meas{x: N(x)=t}, t=0..6.  S_r(E)=E[C(N,r)]=sum_t C(t,r) p_t  (factorial moments).
  Dual g(t)=sum_r y_r C(t,r) >= 1[t=0] on t in {0..6}  =>  measS7(E) <= L_y(E)=sum_r y_r S_r(E)
  for EVERY E (PROVED Bonferroni).  g expands in binary Krawtchouk K_j(t;6) with c_j>=0
  (Delsarte-positive):  g(t)=sum_j c_j K_j(t),  L_y(E)=sum_j c_j M_j(E),  M_j(E):=E[K_j(N)].
    k=8:    c=(1/16, 0, 1/40, 0, 3/80, 0, 0)        (EVEN only: j=0,2,4)
    k=9,10: c=(1/12, 1/72, 1/36, 1/48, 0, 0, 0)     (j=0,1,2,3)
    k=11-13:c=(1/8,  1/24, 1/24, 0,   0, 0, 0)       (j=0,1,2)

THE OPEN CLAIM (THREAD B): max_E L_y(E) = L_y(consec) <= cap_k.

THIS SCRIPT ATTACKS IT THROUGH THE KRAWTCHOUK MOMENTS M_j (the clean Delsarte basis):

(1) KRAWTCHOUK-MOMENT EXTREMALITY.  L_y=sum_j c_j M_j with c_j>=0 and M_0=1.  So a
    SUFFICIENT (stronger) statement is: consec maximizes EACH M_j(E) with c_j>0.
    We test, per k, whether consec maximizes every binding M_j over the bounded-spread
    window.  If yes for the dangerous rows, the LP saturation has a CLEAN per-moment proof:
    each binding Krawtchouk moment is consec-extremal (one inequality each), and a
    nonnegative combination preserves it.

(2) COMPLEMENTARY SLACKNESS.  The dual g(t) is TIGHT (g(t)=1[t=0], i.e. g(t)=0 for t!=0)
    exactly at the INTEGER ROOTS of g.  Primal complementary slackness: the LP-extremizing
    occupancy p_t is supported on {t: g(t) binding}.  k=8 g-roots: {1,2,4,5} plus t=0
    (g(0)=1 active).  So the LP optimum's support is {0,1,2,4,5}; t in {3,6} carry slack.
    We compute consec's actual p_t and ask: is consec's p_t SUPPORTED on the binding set,
    i.e. p_3(consec)=p_6(consec)=0?  (If so, consec SATURATES the LP -> L_y(consec)=measS7
    would follow from LP duality.)  We report the actual support and the slack mass.

(3) LP TIGHTNESS / does L_y(consec) == measS7(consec)?  We compute the gap
    L_y(consec)-measS7(consec) exactly and decompose it as sum_t g(t) p_t over the SLACK
    cells (t with g(t)>1[t=0]).  This says EXACTLY why the bound is not tight at consec
    (which slack cells carry the lost mass).

(4) ANTI-MDS.  consec has a support-2 relation (2*e_1 = e_2, i.e. 2*1=2).  We compute the
    relation-code min distance d(E) and correlate it with M_j and with L_y.  anti-MDS
    (d=2, max low-weight shells) <=> max M_j.

OUTPUT: exact tables.  HONEST about which rows the per-moment extremality holds on.
"""
import sys, itertools, random
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
sys.stdout.reconfigure(line_buffering=True)
random.seed(621621)
H = F(1, 14)

# ---------- binary Krawtchouk on 6 inner sectors ----------
def Kraw(j, t, n=6):
    return sum((-1)**i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))
KTAB = [[F(Kraw(j, t)) for t in range(7)] for j in range(7)]

# ---------- occupancy law p_t(E), exact ----------
def occupancy(E):
    """p_t for t=0..6 (t=#missed inner sectors), exact via one breakpoint pass."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)  # includes 0 from e=0
        free = 6 - len([s for s in hit if s != 0])
        p[free] += L
    return p

def Svec_from_p(p):
    return [sum(comb(t, r) * p[t] for t in range(7)) for r in range(7)]

def Mvec_from_p(p):
    """M_j = E[K_j(N)] = sum_t K_j(t) p_t."""
    return [sum(KTAB[j][t] * p[t] for t in range(7)) for j in range(7)]

def measS7_from_p(p): return p[0]

# ---------- caps (gp danger-zone) ----------
def danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def mgmerge(iv):
    iv = sorted(iv); o = []
    for a, b in iv:
        if o and a <= o[-1][1]: o[-1] = (o[-1][0], max(o[-1][1], b))
        else: o.append((a, b))
    return o
def measGP(P):
    if not P: return F(1)
    dz = mgmerge([iv for u in P for iv in danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev: s += a - prev
        prev = max(prev, b)
    if prev < 1: s += 1 - prev
    return s
@lru_cache(None)
def cap(k):
    psz = 13 - k
    if psz == 0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1, 14), psz))

# ---------- duals: y (factorial) and c (Krawtchouk) per k ----------
DUAL_G = {  # g(t) values t=0..6
    8:  [F((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)],
    9:  [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)],
    10: [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)],
    11: [F((t-3)*(t-4), 12) for t in range(7)],
    12: [F((t-3)*(t-4), 12) for t in range(7)],
    13: [F((t-3)*(t-4), 12) for t in range(7)],
}
def kraw_expand(g):
    Aug = [[KTAB[j][t] for j in range(7)] + [g[t]] for t in range(7)]
    for col in range(7):
        piv = next(r for r in range(col, 7) if Aug[r][col] != 0)
        Aug[col], Aug[piv] = Aug[piv], Aug[col]
        pv = Aug[col][col]; Aug[col] = [x / pv for x in Aug[col]]
        for r in range(7):
            if r != col and Aug[r][col] != 0:
                f = Aug[r][col]; Aug[r] = [Aug[r][i] - f * Aug[col][i] for i in range(8)]
    return [Aug[i][7] for i in range(7)]
DUAL_C = {k: kraw_expand(g) for k, g in DUAL_G.items()}
def L_y(p, k):  # = sum_t g(t) p_t = sum_j c_j M_j
    g = DUAL_G[k]
    return sum(g[t] * p[t] for t in range(7))

# ---------- relation code min distance (anti-MDS test) ----------
def relation_min_dist(E, box=3):
    """min support of nonzero integer relation sum n_i e_i = 0, |n_i|<=box (truncated)."""
    Enz = [e for e in E if e != 0]
    best = 99
    # search small-support relations by trying subsets
    idx = list(range(len(Enz)))
    for r in range(2, min(5, len(Enz)) + 1):
        for sub in itertools.combinations(idx, r):
            es = [Enz[i] for i in sub]
            # find nonzero integer combo summing to 0 with |coeff|<=box, full support
            found = False
            for coeffs in itertools.product(range(-box, box + 1), repeat=r):
                if all(c != 0 for c in coeffs) and sum(c * e for c, e in zip(coeffs, es)) == 0:
                    found = True; break
            if found:
                best = min(best, r); break
        if best <= r: break
    return best

def banner(t): print("\n" + "=" * 96 + f"\n{t}\n" + "=" * 96)

# ================================================================================
banner("(0) The duals: g(t), factorial y, Krawtchouk c.  Verify c>=0 (Delsarte) & L_y forms agree.")
for k in [8, 9, 10, 11, 12, 13]:
    g = DUAL_G[k]; c = DUAL_C[k]
    roots = [t for t in range(7) if g[t] == 0]
    binding = [t for t in range(7) if (g[t] == (1 if t == 0 else 0))]  # tight cells
    print(f"  k={k}: g(t)={[str(x) for x in g]}")
    print(f"       c_j={[str(x) for x in c]}  all c>=0:{all(x>=0 for x in c)}  nonzero-j:{[j for j in range(7) if c[j]!=0]}")
    print(f"       g-roots (t with g(t)=0): {roots}   tight cells (g=1[t=0]): {binding}")

# ================================================================================
banner("(1) KRAWTCHOUK-MOMENT EXTREMALITY: does consec maximize each binding M_j(E)?")
print("  For each k, the binding moments are j with c_j>0.  We sweep bounded spread and")
print("  report, per binding j, whether consec is the max of M_j (and the # of beaters).")
WINDOWS = {8: 16, 9: 15, 10: 14, 11: 13}
for k in [8, 9, 10, 11]:
    maxE = WINDOWS[k]
    Ec = list(range(k)); pc = occupancy(Ec); Mc = Mvec_from_p(pc); Lc = L_y(pc, k)
    binding_j = [j for j in range(7) if DUAL_C[k][j] != 0 and j >= 1]
    beat = {j: 0 for j in binding_j}; Ly_beat = 0; nset = 0
    worstM = {j: Mc[j] for j in binding_j}; worstLy = Lc
    argM = {j: 'consec' for j in binding_j}
    for rest in itertools.combinations(range(1, maxE + 1), k - 1):
        E = [0] + list(rest); p = occupancy(E); M = Mvec_from_p(p); Ly = L_y(p, k); nset += 1
        for j in binding_j:
            if M[j] > worstM[j] + F(0): worstM[j] = M[j]; argM[j] = E
            if M[j] > Mc[j]: beat[j] += 1
        if Ly > worstLy: worstLy = Ly
        if Ly > Lc: Ly_beat += 1
    print(f"\n  k={k}  (maxE<={maxE}, {nset} sets):  binding j={binding_j}")
    for j in binding_j:
        print(f"    M_{j}: consec={float(Mc[j]):.5f}  windowmax={float(worstM[j]):.5f}  "
              f"consec_is_max={beat[j]==0}  beaters={beat[j]}  (argmax={argM[j] if beat[j] else 'consec'})")
    print(f"    L_y: consec={float(Lc):.5f}  windowmax={float(worstLy):.5f}  consec_is_max={Ly_beat==0}  "
          f"beaters={Ly_beat}  cap={float(cap(k)):.5f}  L_y(consec)<=cap:{Lc<=cap(k)}")

# ================================================================================
banner("(2) COMPLEMENTARY SLACKNESS: consec's occupancy p_t vs the dual binding set.")
print("  Dual tight cells = {t: g(t)=1[t=0]}.  LP-optimal primal supported there.")
print("  We list consec p_t, mark which cells are slack (g(t)>1[t=0]), and the mass consec")
print("  puts on slack cells (= why L_y(consec) exceeds measS7(consec)).")
for k in [8, 9, 10, 11, 12, 13]:
    Ec = list(range(k)); pc = occupancy(Ec); g = DUAL_G[k]
    tight = [t for t in range(7) if g[t] == (1 if t == 0 else 0)]
    slack = [t for t in range(7) if t not in tight]
    slackmass = sum(pc[t] for t in slack)
    onslack_Lcontrib = sum(g[t] * pc[t] for t in slack)
    print(f"\n  k={k}: consec p_t = {[str(x) for x in pc]}")
    print(f"       float p_t   = {[round(float(x),4) for x in pc]}")
    print(f"       tight cells={tight}  slack cells={slack}")
    print(f"       consec mass on slack cells = {float(slackmass):.5f}  "
          f"(p_t on slack: {[(t, round(float(pc[t]),4)) for t in slack]})")
    print(f"       => consec saturates LP support? (zero mass on slack): {slackmass==0}")

# ================================================================================
banner("(3) LP TIGHTNESS: gap L_y(consec) - measS7(consec), decomposed over slack cells.")
print("  L_y - measS7 = sum_t g(t)p_t - p_0 = sum_{t>=1} g(t)p_t  (since g(0)=1).")
print("  This is the EXACT non-tightness.  measS7=L_y iff p_t=0 for all t>=1 with g(t)>0.")
for k in [8, 9, 10, 11, 12, 13]:
    Ec = list(range(k)); pc = occupancy(Ec); g = DUAL_G[k]
    meas = pc[0]; Ly = L_y(pc, k); gap = Ly - meas
    contrib = [(t, g[t] * pc[t]) for t in range(1, 7) if g[t] * pc[t] != 0]
    print(f"\n  k={k}: measS7={float(meas):.5f}  L_y={float(Ly):.5f}  gap={float(gap):.5f}={gap}")
    print(f"       gap contributions g(t)p_t (t>=1): {[(t, str(v)) for t,v in contrib]}")
    print(f"       tight at consec (gap=0)? {gap==0}")

# ================================================================================
banner("(4) ANTI-MDS: relation min-distance d(E) vs L_y(E) -- does d=2 (anti-MDS) max L_y?")
for k in [8, 9, 10]:
    maxE = WINDOWS.get(k, 14)
    Ec = list(range(k)); dc = relation_min_dist(Ec); pc = occupancy(Ec); Lc = L_y(pc, k)
    print(f"\n  k={k}: consec d={dc}  L_y(consec)={float(Lc):.5f}")
    # sample: best L_y per distance class
    bydist = {}
    cnt = 0
    for rest in itertools.combinations(range(1, maxE + 1), k - 1):
        E = [0] + list(rest); cnt += 1
        if cnt > 4000: break
        d = relation_min_dist(E); p = occupancy(E); Ly = L_y(p, k)
        if d not in bydist or Ly > bydist[d][0]: bydist[d] = (Ly, E)
    for d in sorted(bydist):
        Ly, E = bydist[d]
        print(f"       d={d}: max L_y={float(Ly):.5f} at {E}")

print("\nDONE.")
