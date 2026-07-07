"""
Toward a PROVABLE lower bound on mu_{1/7}(E) (the load-bearing tail lemma (A')).
boxeph-2026-07-07-S1 cont., HYP-4760.

The fleet converged: the one honest open lemma is mu_{1/7}(E) >= mu_{1/7}(AP_k)
(sharp AP-minimality of the theta=1/7 tail; opus-S134 has the AP side via the
Farey roof, death-star hardened it numerically). Reverse-Markov/E[maxgap] is dead.

I test two candidate PROVABLE routes to a universal floor, and locate where the
mu_{1/7} mass lives:

ROUTE A (single-gap tail):  maxgap >= gap@0 pointwise => mu_{1/7}(E) >= P(gap@0>1/7).
   gap@0 = R+L, R=min_i frac(e_i x), L=min_i frac(-e_i x).  If inf_E P(gap@0>1/7)
   is comfortably positive, that's a clean single-gap route (no full max needed).

ROUTE B (universal q<=6 grid window):  at x=p/q the config lies on the q-slot
   grid {j/q}, so maxgap(p/q) >= 1/q > 1/7 for q<=6, for EVERY integer family.
   By continuity maxgap>1/7 on a neighborhood.  Measure = a rigorous universal
   lower bound (family-dependent width).  How big is it, and is it dilation-robust?
"""
from fractions import Fraction as F

# ---------- exact P(gap@0 > thr) ----------
def P_gap0_gt(E, thr=F(1,7)):
    """EXACT meas{x: R(x)+L(x) > thr}, R=min frac(e x), L=min frac(-e x)."""
    E = [e for e in E if e != 0]
    # breakpoints: R,L are piecewise linear; order changes where e_i x - e_j x in Z
    # and wraps e_i x in Z. Plus we need where R+L crosses thr (linear in each piece).
    bps = {F(0), F(1)}
    n = len(E)
    for i in range(n):
        ai = abs(E[i])
        for m in range(0, ai+1):
            bps.add(F(m, ai))
        for j in range(i+1, n):
            d = abs(E[i]-E[j])
            for m in range(0, d+1):
                bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        # R = min_i frac(e_i x): argmin line eR x + bR
        rbest = None; lbest = None
        for e in E:
            flr = (e*mid).__floor__()
            v = e*mid - flr
            if rbest is None or v < rbest[0]:
                rbest = (v, e, flr)   # R = e x - flr
            fll = (-e*mid).__floor__()
            w = -e*mid - fll
            if lbest is None or w < lbest[0]:
                lbest = (w, -e, fll)  # L = -e x - fll
        _, eR, flR = rbest
        _, eL, flL = lbest
        # S(x) = R+L = (eR+eL) x - (flR+flL); linear on (a,b). Where S>thr?
        c = eR + eL; b0 = -(flR + flL)
        # S(x)=c x + b0 > thr on (a,b)
        if c == 0:
            if b0 > thr: total += (b-a)
        elif c > 0:
            x0 = (thr - b0)/c; lo = max(a, x0)
            if lo < b: total += (b - lo)
        else:
            x0 = (thr - b0)/c; hi = min(b, x0)
            if a < hi: total += (hi - a)
    return total

# ---------- exact mu_{1/7} (reuse) ----------
def mu_exact(E, thr=F(1,7)):
    E = list(E)
    denoms = set()
    for i in range(len(E)):
        for j in range(len(E)):
            if E[i] != E[j]:
                denoms.add(abs(E[i]-E[j]))
        denoms.add(abs(E[i]) if E[i] != 0 else 1)
    bps = {F(0), F(1)}
    for d in denoms:
        for m in range(0, d+1):
            bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        fl = {e: (e*mid).__floor__() for e in E}
        order = sorted(E, key=lambda e: e*mid - fl[e])
        m = len(order)
        subs = []
        for s in range(m):
            e1 = order[s]; e2 = order[(s+1) % m]
            if s < m-1:
                c = F(e2-e1); b0 = F(-(fl[e2]-fl[e1]))
            else:
                c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]]) + 1)
            if c == 0:
                if b0 > thr: subs.append((a, b))
            elif c > 0:
                x0 = (thr-b0)/c; lo = max(a, x0)
                if lo < b: subs.append((lo, b))
            else:
                x0 = (thr-b0)/c; hi = min(b, x0)
                if a < hi: subs.append((a, hi))
        subs.sort(); cur = None
        for lo, hi in subs:
            if lo >= hi: continue
            if cur is None: cur = [lo, hi]
            elif lo <= cur[1]: cur[1] = max(cur[1], hi)
            else:
                total += cur[1]-cur[0]; cur = [lo, hi]
        if cur is not None:
            total += cur[1]-cur[0]
    return total

ap13 = list(range(1,14))
gw   = [1,2,3,4,5,6,7,8,9,10,11,13,24]
print("="*66)
print("ROUTE A: single-gap tail  mu_{1/7}(E) >= P(gap@0 > 1/7)")
print("="*66)
for nm,E in [("AP",ap13),("GW",gw)]:
    mu = mu_exact(E); pg = P_gap0_gt(E)
    print(f"  {nm:4s}  mu_1/7={float(mu):.4f}  P(gap0>1/7)={float(pg):.4f}  ratio={float(pg/mu):.3f}")
print("  (AP: gap0=maxgap a.e. => P(gap0>1/7)=mu; other families P<mu)")

# adversarial inf of P(gap@0>1/7)
import random
random.seed(7)
def P_gap0_fast(E, G=40000, thr=1/7):
    c=0
    for aa in range(G):
        x=(aa+0.5)/G
        R=2.0;L=2.0
        for e in E:
            fe=(e*x)%1.0
            if fe<R:R=fe
            g=1.0-fe
            if g<L:L=g
        if R+L>thr:c+=1
    return c/G
worst=(1.0,None); below=0
cands=[ap13,gw,[2*j for j in range(1,14)],[1,2,3,4,5,6,20,21,22,23,24,25,26]]
for _ in range(250):
    cands.append(sorted(random.sample(range(1,60),13)))
# structured spread adversaries (these minimize gap0)
for d in range(2,10):
    for a0 in range(1,d):
        E=[a0+d*j for j in range(13)]
        if len(set(E))==13: cands.append(E)
for E in cands:
    v=P_gap0_fast(E)
    if v<worst[0]: worst=(v,E)
    if v<1/7: below+=1
print(f"  adversarial ({len(cands)} families): inf P(gap0>1/7) = {worst[0]:.4f} at {worst[1]}")
print(f"    {below} below 1/7; ROUTE A {'CLEARS m_P=0.057' if worst[0]>0.057 else 'may fail'} "
      f"(needs >= m_P~0.057)")

print("\n" + "="*66)
print("ROUTE B: universal q<=6 grid window (rigorous, family-indep at rationals)")
print("  maxgap(p/q) >= 1/q > 1/7 for q<=6 for EVERY integer family.")
print("="*66)
# measure of mu_{1/7} within +-w of some p/q, q<=6, vs total mu -- how much mass is 'universal'?
def mass_near_smallq(E, qmax=6, G=60000, thr=1/7):
    # fraction of mu_{1/7} mass with x within the Farey-cell of a q<=6 rational
    smallq_rats=sorted(set(F(p,q) for q in range(1,qmax+1) for p in range(0,q+1)))
    smallq=[float(r) for r in smallq_rats]
    tot=0;near=0
    def near_smallq(x):
        return min(abs(x-r) for r in smallq)
    for aa in range(G):
        x=(aa+0.5)/G
        pts=sorted((e*x)%1.0 for e in E)
        mg=0.0
        for i in range(1,len(pts)):
            g=pts[i]-pts[i-1]
            if g>mg:mg=g
        w=1.0-pts[-1]+pts[0]
        if w>mg:mg=w
        if mg>thr:
            tot+=1
            # 'near' = closest small-q rational within 1/(2*qmax*... ) heuristic band
            if near_smallq(x) < 1.0/(2*7):  # within 1/14 of a q<=6 rational
                near+=1
    return (near/G, tot/G)
for nm,E in [("AP",ap13),("GW",gw),("primes",[2,3,5,7,11,13,17,19,23,29,31,37,41]),
             ("spreadAP{6+5j}",[6+5*j for j in range(13)])]:
    nearf,totf=mass_near_smallq(E)
    print(f"  {nm:15s} mu_1/7={totf:.4f}  mass within 1/14 of q<=6 rational={nearf:.4f} "
          f"({100*nearf/totf:.0f}% of mu)")
print("\n  => if most mu-mass is near q<=6 rationals for the WORST (near-AP) families,")
print("     the universal grid window is the load-bearing part; saturated/spread")
print("     families get extra mass elsewhere (decorrelation).")
