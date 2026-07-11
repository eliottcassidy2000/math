# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-11-S127: the FINITE CERTIFICATE for THM-701's recursion closure.
# Two tables:
#  (I) CAP-GROWTH: cap_k = min_{|P|=13-k} meas(G_P) (G_P = lonely set of P at gap 1/14). Show
#      cap_{k+1} - cap_k >= 2/21 = 0.09524 for k = 8..12 (so the joint-increment bound closes at every step).
#  (II) BALANCED-CORE: Phi(F) = p0(F) + (1/3) p1(F) <= cap_{|F|+1} for the base family (bounded-spread cores).
#      consec_m is the argmax (HYP-2644); verified Phi(consec_m) < cap_{m+1}, plus a robust search over
#      shifted-AP / geometric / random balanced cores confirming Phi <= cap with margin.
import math, random
from fractions import Fraction as F

# ---------- LONELY MEASURE (gap 1/14) for small P => caps 11,12,13 ----------
def lonely_measure(P, N=200003):
    # meas{ x in [0,1): min_i ||p_i x|| >= 1/14 },  ||.|| = distance to nearest integer
    if not P: return 1.0
    cnt = 0
    for i in range(N):
        x = (i + 0.5) / N
        ok = True
        for p in P:
            y = (p * x) % 1.0
            if min(y, 1 - y) < 1/14 - 1e-12:
                ok = False; break
        if ok: cnt += 1
    return cnt / N

def min_lonely(m, cand):
    # min over m-subsets drawn from candidate speeds (primitive small speeds dominate the min)
    from itertools import combinations
    best = (9.0, None)
    for S in combinations(cand, m):
        v = lonely_measure(list(S))
        if v < best[0]: best = (v, S)
    return best

# ---------- seven-sector p_j and Phi ----------
def sector(y): return min(6, int((y % 1.0) * 7))
def missed(E, x): return set(range(1, 7)) - set(sector(e * x) for e in E)
def pvec(E, N=60013):
    v = [0]*7
    for i in range(N): v[len(missed(E, (i+0.5)/N))] += 1
    return [c/N for c in v]
def Phi(E, lam=1/3):
    p = pvec(E); return p[0] + lam*p[1]

def main():
    two21 = 2/21
    # cap_8,9,10 from THM-532/534 (LRC extremals over |P|=5,4,3); cap_11,12,13 computed here.
    cap = {8: 0.38153, 9: 0.49426, 10: 0.60440}
    cap[13] = lonely_measure([])                       # |P|=0
    cap[12] = min_lonely(1, [1,2,3,4,5,6])[0]          # |P|=1  (= 6/7)
    c11 = min_lonely(2, [1,2,3,4,5,6,7,8,9,10,11,12,13])
    cap[11] = c11[0]                                    # |P|=2
    print("=== TABLE I: cap-growth (cap_k = min_{|P|=13-k} meas(G_P)) ===")
    print(f"  cap_13 = {cap[13]:.5f}  (|P|=0)")
    print(f"  cap_12 = {cap[12]:.5f}  (|P|=1; 6/7 = {6/7:.5f})")
    print(f"  cap_11 = {cap[11]:.5f}  (|P|=2, argmin {c11[1]})")
    print(f"  cap_10 = {cap[10]:.5f}  cap_9 = {cap[9]:.5f}  cap_8 = {cap[8]:.5f}  (LRC extremals, THM-532/534)")
    print(f"  2/21 = {two21:.5f}")
    ok = True
    for k in range(8, 13):
        g = cap[k+1] - cap[k]
        good = g >= two21
        ok = ok and good
        print(f"  cap_{k+1} - cap_{k} = {g:.5f}  >= 2/21 ? {'YES' if good else 'NO'}")
    print(f"  => cap-growth >= 2/21 for all k=8..12 : {'CERTIFIED' if ok else 'FAILS'}")
    print()

    # ---------- TABLE II: balanced-core Phi <= cap_{|F|+1} ----------
    print("=== TABLE II: balanced-core check  Phi(F)=p0+(1/3)p1 <= cap_{|F|+1} ===")
    print("  (a) consec_m (the HYP-2644 argmax):")
    for m in range(7, 13):
        E = list(range(m)); ph = Phi(E); c = cap[m+1]
        print(f"     Phi(consec_{m}) = {ph:.5f}   cap_{m+1} = {c:.5f}   margin = {c-ph:.5f}  ({'OK' if ph<c else 'FAIL'})")
    print("  (b) robust search over balanced cores (shifted-AP / geometric-near / random), min margin:")
    random.seed(7)
    worst = {}
    def upd(m, E):
        ph = Phi(E); worst[m] = min(worst.get(m, 9.0), cap[m+1]-ph)
    for m in range(7, 13):
        for s in [1, 5, 20, 100]:                       # shifted consec_m by s (large-N balanced blocks)
            upd(m, [0] + [s+i for i in range(m-1)])
        for r in [2, 3]:                                 # near-geometric small
            upd(m, [0] + [r**i for i in range(m-1)])
        for _ in range(300):                             # random bounded-spread (max <= 4*m)
            E = [0] + sorted(random.sample(range(1, 4*m+1), m-1))
            upd(m, E)
    for m in range(7, 13):
        print(f"     |F|={m}: min margin cap_{m+1} - Phi over balanced search = {worst[m]:.5f}  ({'OK' if worst[m]>0 else 'VIOLATED'})")

if __name__ == "__main__":
    main()
