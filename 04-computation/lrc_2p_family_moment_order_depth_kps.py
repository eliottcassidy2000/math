"""
lrc_2p_family_moment_order_depth_kps.py  (kind-pasteur-2026-06-28-S31aq)

VERIFY the HYP-3216 family prediction: for LRC(2p) (apex prime p), the cover-bound
moment-order ladder has DEPTH (p-1)/2 (= the cyclotomic degree), i.e. the minimal
moment order t* that closes the binding (consec) row over the p sectors is
(p+1)/2 [the apex node degree], and lower orders FAIL.

For n=2p the relaxed threshold is 1/p; p sectors [j/p,(j+1)/p); N = #empty inner
sectors (j=1..p-1). The moment-LP degree-t bound:
   U_t(E) = max{ q_0 : sum_i C(i,r) q_i = S_r(E), r=0..t, q_i>=0 }.
t*(E) = minimal t with U_t(E) "tight" (= the actual max coverage, i.e. the bound stops
overshooting). Find t* for the binding consec clusters at p=3,5,7; check t* = (p+1)/2.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb

def sector_of(x, p): return int((x % 1) * p)

def Ndist(E, p):
    """q[t] = meas{N=t}, N = #empty inner sectors of the p-sector partition, exact."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, p * e + 1): b.add(F(m, p * e))
    b = sorted(b); q = [F(0)] * p   # N in 0..p-1
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        cov = set(sector_of(e * ((x0 + x1) / 2), p) for e in E)
        t = (p - 1) - len([s for s in range(1, p) if s in cov])  # empty inner
        if 0 <= t <= p - 1: q[t] += x1 - x0
    return q

def moment_LP_Ut(q, t, p):
    """max q_0 s.t. sum_i C(i,r) q_i = S_r (r<=t), q>=0, i in 0..p-1. Vertex enum (support <= t+1)."""
    pdim = p
    S = [sum(comb(i, r) * q[i] for i in range(pdim)) for r in range(t + 1)]
    best = F(0)
    # optimum supported on a subset of size <= t+1 INCLUDING index 0 (to max q_0)
    others = list(range(1, pdim))
    for extra in itertools.combinations(others, min(t, len(others))):
        supp = [0] + list(extra)
        # solve linear system: for r=0..t, sum_{i in supp} C(i,r) q_i = S_r
        # (|supp| = t+1 unknowns, t+1 equations)
        if len(supp) != t + 1:
            # underdetermined: pad — skip (take t+1 exactly)
            continue
        import numpy as np
        A = [[comb(i, r) for i in supp] for r in range(t + 1)]
        Sb = [float(x) for x in S]
        try:
            sol = np.linalg.solve(np.array(A, float), np.array(Sb, float))
        except Exception:
            continue
        if all(x >= -1e-9 for x in sol):
            q0 = sol[0]
            if q0 > best: best = F(0) + round(q0, 9) if False else max(best, F(0)) ; best = max(float(best), q0) if isinstance(best, float) else q0
    # robust float version
    bestf = 0.0
    Sf = [float(x) for x in S]
    import numpy as np
    for extra in itertools.combinations(others, min(t, len(others))):
        supp = [0] + list(extra)
        if len(supp) != t + 1: continue
        A = np.array([[comb(i, r) for i in supp] for r in range(t + 1)], float)
        try:
            sol = np.linalg.solve(A, np.array(Sf, float))
        except Exception:
            continue
        if all(x >= -1e-9 for x in sol):
            if sol[0] > bestf: bestf = sol[0]
    return bestf

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    print("Verify: LRC(2p) moment-order DEPTH = (p-1)/2; apex-node closing order t* = (p+1)/2")
    for p in (3, 5, 7):
        print(f"\n=== p={p} (n={2*p}), cyclotomic degree (p-1)/2={(p-1)//2}, predicted apex t*=(p+1)/2={(p+1)//2} ===")
        # binding (apex) cluster: the smallest k where coverage is near-tight. Use consec_k for a range of k.
        apex_depths = []
        for k in range(2*p - 6 if p > 3 else 2, 2*p):
            consec = tuple(range(k))
            q = Ndist(consec, p)
            actual_q0 = float(q[0])
            cap_k = float(F(comb(k+1, 2), comb(2*p, 2)))   # pair-Pascal cap = C(k+1,2)/C(2p,2)
            Us = [round(moment_LP_Ut(q, t, p), 4) for t in range(1, p)]
            # DEPTH = minimal t with U_t <= cap_k (the cover bound closes at that order)
            depth = next((t for t in range(1, p) if Us[t-1] <= cap_k + 1e-9), p-1)
            apex_depths.append(depth)
            print(f"  consec_{k}: q0={actual_q0:.4f} cap_k=C({k+1},2)/{comb(2*p,2)}={cap_k:.4f}  "
                  f"U_t={Us}  -> DEPTH(U_t<=cap)={depth}")
        apex = max(apex_depths)
        pred = (p+1)//2
        print(f"  APEX DEPTH (deepest binding row) = {apex}; predicted (p+1)/2 = {pred}  -> {'CONFIRMS' if apex==pred else 'MISMATCH'} the family law")
