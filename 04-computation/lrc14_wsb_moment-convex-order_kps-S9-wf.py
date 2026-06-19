#!/usr/bin/env python3
"""
lrc14_wsb_moment-convex-order_kps-S9-wf.py   (kind-pasteur-2026-06-19-S9, Angle: moment-convex-order)

CRUX (k=8 form): consec {0..7} maximizes  L_y(E) = P(N=0) + (1/10) P(N=3) + P(N=6)
  where N_E(x) = #{ sectors j in {1..6} EMPTY of the orbit {frac(e_i x)} } in {0..6} (sector 0
  always hit by e=0).  meas(S7(E)) = P(N_E=0) = p_0(E) is what we actually want bounded by cap_k.

We attack via CONVEX / STOCHASTIC ORDER on the empty-sector count N.

PLAN
 (P0) Build the EXACT law p_0..p_6 of N_E for any E (exact rationals, breakpoint engine).
      Cross-check: sum p_i = 1, p_0 = meas(S7), S_r = E[C(N,r)] matches Sr() from the canon tools.
 (P1) THM-536 B2 reinterpreted as a COUPLING: for E subset {0..N}, N_E(x) >= N_{AP_{N+1}}(x)
      POINTWISE (E hits a subset of AP's sectors => more empty). => N_E >=_st N_{AP_{N+1}} (usual
      stochastic order, full distribution).  Confirm exactly.  This is the wrong direction for p_0
      vs consec_k (it dominates by the BIGGER AP). Quantify the span gap.
 (P2) The real question: between two SAME-k sets, what order makes consec extremal for the
      g-weighted functional?  g(0)=g(6)=1, g(3)=1/10, else 0 -- NOT monotone, NOT convex, NOT
      concave.  Identify the order class that g respects and test whether N_consec is extremal in it.
      Compute exact law of N for consec_8 and a thorough bank; record (a) usual stoch order pairs,
      (b) convex order pairs, (c) the sign pattern of g vs each candidate order.
 (P3) Exact law of N_{AP} via the Sturmian walk (THM-536): partial sums S_e=floor(e*theta) mod 7;
      N = 6 - #{distinct nonzero residues visited}.  Get p_i(AP_k) exactly and inspect its shape.
"""
import sys, itertools, functools
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------------
# EXACT law of N_E:  p_t = meas{ x in [0,1): exactly t of sectors {1..6} empty }.
# Breakpoint engine: the function x -> set-of-hit-sectors is piecewise constant; breakpoints are
# all x where some frac(e*x) crosses a sector boundary j/7.
# ----------------------------------------------------------------------------
def law_N(E):
    """Return list p[0..6] (Fractions) = meas{ N_E == t }.  N counts EMPTY sectors among {1..6}."""
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for j in range(7):                      # boundaries j/7, j=0..6 (and 1==7/7)
            c = F(j, 7)
            for m in range(e):                  # frac(e x)=c  => x=(c+m)/e, m=0..e-1
                bps.add((c + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    p = [F(0)] * 7
    for i in range(len(bps)):
        x0 = bps[i]
        x1 = bps[i + 1] if i + 1 < len(bps) else F(1)
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        hit = set()
        for e in E:
            s = int((e * xm * 7) % 7)            # exact: e*xm rational; floor(7*frac(e xm))
            # compute sector = floor(7 * frac(e*xm)) exactly
            fr = (e * xm) % 1
            sec = int(fr * 7)                    # fr in [0,1); 7*fr in [0,7); floor via int on Fraction*int? do exact:
            num = fr.numerator * 7
            den = fr.denominator
            sec = num // den
            hit.add(sec)
        empties = 6 - len([s for s in hit if 1 <= s <= 6])  # sector 0 always hit by e=0
        # (sectors {1..6}; count how many are NOT hit)
        nempty = len([j for j in range(1, 7) if j not in hit])
        p[nempty] += x1 - x0
    return p

# canonical Sr / Ly tools (copied from prompt) for cross-check
def J(A, E):
    E = sorted(set(E)); arcs = [(F(j, 7), F(j + 1, 7)) for j in A]; bp = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for (a, b) in arcs:
            for end in (a, b):
                m = 0
                while True:
                    xv = (end + m) / e
                    if xv >= 1: break
                    if xv >= 0: bp.add(xv)
                    m += 1
    bp = sorted(b for b in bp if 0 <= b < 1); tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if all(not any(a < ((e * mid) % 1) < b for (a, b) in arcs) for e in E): tot += hi - lo
    return tot
def Sr(E, r): return sum((J(set(A), E) for A in itertools.combinations(range(1, 7), r)), F(0))
def Sr_from_law(p, r): return sum(F(comb(t, r)) * p[t] for t in range(7))
def Ly_from_law(p, k):
    if k == 8:   return p[0] + F(1, 10) * p[3] + p[6]
    if k in (9, 10):
        g = lambda t: F(-(t - 2) * (t - 3) * (t - 6), 36)
        return sum(g(t) * p[t] for t in range(7))
    g = lambda t: F((t - 3) * (t - 4), 12)
    return sum(g(t) * p[t] for t in range(7))

# cap_k engine
H = F(1, 14)
def danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def measGP(P):
    dz = sorted([iv for u in P for iv in danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev: s += a - prev
        prev = max(prev, b)
    if prev < 1: s += 1 - prev
    return s
@functools.lru_cache(None)
def cap(k): return min(measGP(P) for P in itertools.combinations(range(1, 14), 13 - k))

# ----------------------------------------------------------------------------
print("="*78)
print("(P0) EXACT law of N + cross-checks (consec_k, k=8,9,10)")
print("="*78)
for k in (8, 9, 10):
    E = list(range(k))
    p = law_N(E)
    tot = sum(p)
    p0 = p[0]
    ly = Ly_from_law(p, k)
    capk = cap(k)
    # cross-check Sr
    ok = all(Sr_from_law(p, r) == Sr(E, r) for r in range(5))
    print(f" k={k} consec: law p0..p6 = {[str(x) for x in p]}")
    print(f"        sum={tot}  p0=meas(S7)={float(p0):.5f}  L_y={float(ly):.5f}  cap={float(capk):.5f}"
          f"  {'CLOSES' if ly<=capk else 'OVER!'}  Sr-xcheck={'OK' if ok else 'FAIL'}")

print()
print("="*78)
print("(P1) THM-536 B2 as a coupling:  E subset {0..N}  =>  N_E(x) >= N_{AP_{N+1}}(x) pointwise")
print("="*78)
# verify pointwise domination on a fine common breakpoint grid for a few E
def hit_sectors(E, xm):
    hit = set()
    for e in E:
        fr = (e * xm) % 1
        hit.add((fr.numerator * 7) // fr.denominator)
    return hit
def Nval(E, xm):
    hit = hit_sectors(E, xm)
    return len([j for j in range(1, 7) if j not in hit])
def pointwise_dom(E, F2):
    """check N_E(x) >= N_{F2}(x) at all common cells; return (#cells, #violations)."""
    allE = sorted(set(E) | set(F2))
    bps = set([F(0), F(1)])
    for e in allE:
        if e == 0: continue
        for j in range(7):
            c = F(j, 7)
            for m in range(e): bps.add((c + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    nv = 0; nc = 0
    for i in range(len(bps)):
        x0 = bps[i]; x1 = bps[i+1] if i+1 < len(bps) else F(1)
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; nc += 1
        if Nval(E, xm) < Nval(F2, xm): nv += 1
    return nc, nv
tests = [([0,1,3,7], 7), ([0,2,5], 5), ([0,1,4,9], 9), ([0,3,6,10], 10)]
for E, N in tests:
    AP = list(range(N+1))
    nc, nv = pointwise_dom(E, AP)
    print(f"  E={E} (subset {{0..{N}}})  vs AP_{N+1}: cells={nc} violations of N_E>=N_AP = {nv}"
          f"  {'POINTWISE DOM OK' if nv==0 else 'FAILS'}")

print()
print("="*78)
print("(P2) SAME-k order: consec_8 law vs a bank; usual stoch order? convex order? g-respect?")
print("="*78)
k = 8
pc = law_N(list(range(k)))
lyc = Ly_from_law(pc, k)
def cdf(p):
    c=[F(0)]*7; run=F(0)
    for t in range(7): run+=p[t]; c[t]=run
    return c
def tail(p):  # P(N>=t)
    return [sum(p[t:],F(0)) for t in range(7)]
print(f"  consec_8 law: {[str(x) for x in pc]}  L_y={float(lyc):.5f}")
print(f"  consec_8 P(N=0)={float(pc[0]):.4f} P(N=3)={float(pc[3]):.4f} P(N=6)={float(pc[6]):.4f}")
# bank: all primitive [0]+combos of {1..12} size 7
bank = [[0]+list(r) for r in itertools.combinations(range(1,13), k-1)]
n_over = 0          # L_y(E) > L_y(consec) ?  (would refute crux)
n_usual_dom_c = 0   # N_E >=_st N_consec  (consec stoch-min)?
n_usual_dom_E = 0   # N_consec >=_st N_E ?
n_cx_dom = 0        # convex-order comparable?
worst = None
for E in bank:
    p = law_N(E)
    ly = Ly_from_law(p, k)
    if ly > lyc + F(1,10**15):
        n_over += 1
        if worst is None or ly > worst[1]: worst = (E, ly)
    # usual stoch order: compare tails
    tE = tail(p); tC = tail(pc)
    if all(tE[t] >= tC[t] for t in range(7)): n_usual_dom_c += 1   # N_E >=_st N_consec
    if all(tC[t] >= tE[t] for t in range(7)): n_usual_dom_E += 1   # N_consec >=_st N_E
print(f"  bank size={len(bank)}; L_y(E)>L_y(consec): {n_over}  (0 => consec is L_y-max here)")
print(f"  N_E >=_st N_consec for {n_usual_dom_c} sets;  N_consec >=_st N_E for {n_usual_dom_E} sets")
if worst: print(f"  WORST overshoot E={worst[0]} L_y={float(worst[1]):.6f}")

print()
print("="*78)
print("(P3) exact AP law via Sturmian walk (THM-536) and its shape")
print("="*78)
def law_AP_sturmian(k):
    """meas of N for consec_k via theta=7x reparam: residues floor(e theta) mod 7, e=0..k-1.
       N = #{j in 1..6 : j not visited}.  Breakpoints in theta at e*theta integer."""
    bps = set([F(0), F(7)])
    for e in range(1, k):
        for t in range(0, 7*e + 1):
            bps.add(F(t, e))
    bps = sorted(z for z in bps if F(0) <= z < F(7))
    p = [F(0)]*7
    for i in range(len(bps)):
        t0 = bps[i]; t1 = bps[i+1] if i+1 < len(bps) else F(7)
        if t1 <= t0: continue
        tm = (t0+t1)/2
        vis = set()
        for e in range(k):
            v = e*tm
            vis.add((v.numerator // v.denominator) % 7)   # floor(e*theta) mod 7
        nempty = len([j for j in range(1,7) if j not in vis])
        p[nempty] += (t1 - t0)/7    # /7 measure-preserving
    return p
for k in (7,8,9,10):
    pA = law_AP_sturmian(k)
    pdirect = law_N(list(range(k)))
    match = (pA == pdirect)
    print(f"  k={k}: sturmian law={[str(x) for x in pA]}  matches-direct={match}  p0={float(pA[0]):.5f}")

print("\nDONE (stage 1).")
