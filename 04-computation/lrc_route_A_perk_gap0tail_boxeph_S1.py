"""
ROUTE A per-k: inf over k-element families of P(gap@0 > 1/7), vs the union-bound
per-k thresholds T_k (monad-explorer HYP-4787).  boxeph-2026-07-07-S1, HYP-4760.

Route A: mu_{1/7}(E) >= P(gap@0 > 1/7)  (pointwise maxgap >= gap@0).  Fed into the
skeleton union bound G2 >= meas(G_P) + mu_{1/7}(E) - 1, Route A discharges leg k
iff  inf_k P(gap@0>1/7) >= T_k := 1 - meas(G_P^{min,k}) + m_P.

monad-explorer T_k (k=8..12): 0.6185, 0.5057, 0.3956, 0.2747, 0.1429; k=13 (P=empty): m_P~0.0565.

If Route A discharges k=11,12,13 and opus-S131 verified AP-minimality exhaustively
for k=8,9,10, the whole (A') ledger is covered by TWO different, simpler arguments:
a single-gap avoidance tail (large k) + a finite AP-min check (small k).
"""
from fractions import Fraction as F
import random

def P_gap0_gt_exact(E, thr=F(1,7)):
    E = [e for e in E if e != 0]
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
        rbest=lbest=None
        for e in E:
            flr=(e*mid).__floor__(); v=e*mid-flr
            if rbest is None or v<rbest[0]: rbest=(v,e,flr)
            fll=(-e*mid).__floor__(); w=-e*mid-fll
            if lbest is None or w<lbest[0]: lbest=(w,-e,fll)
        _,eR,flR=rbest; _,eL,flL=lbest
        c=eR+eL; b0=-(flR+flL)
        if c==0:
            if b0>thr: total+=(b-a)
        elif c>0:
            x0=(thr-b0)/c; lo=max(a,x0)
            if lo<b: total+=(b-lo)
        else:
            x0=(thr-b0)/c; hi=min(b,x0)
            if a<hi: total+=(hi-a)
    return total

def P_gap0_fast(E, G=30000, thr=1/7):
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

Tk = {8:0.6185, 9:0.5057, 10:0.3956, 11:0.2747, 12:0.1429, 13:0.0565}
random.seed(2026)
print(f"{'k':>3} | {'inf P(gap0>1/7)':>15} | {'T_k':>7} | {'Route A discharges?':>20} | minimizer")
print("-"*90)
results={}
for k in range(8,14):
    worst=(1.0,None)
    cands=[list(range(1,k+1))]
    # spread APs {a+dj} (the gap0 minimizers)
    for d in range(2,12):
        for a0 in range(0,d):
            E=[a0+d*j for j in range(k)]
            E=[e for e in E if e!=0]
            if len(set(E))==k: cands.append(sorted(set(E))[:k] if len(set(E))>=k else None)
    cands=[c for c in cands if c and len(c)==k]
    # random spreads
    for _ in range(400):
        cands.append(sorted(random.sample(range(1,70),k)))
    for E in cands:
        v=P_gap0_fast(E)
        if v<worst[0]: worst=(v,E)
    # exact-verify the found minimizer
    exact=float(P_gap0_gt_exact(worst[1]))
    inf_est=min(worst[0], exact)
    ok="YES" if inf_est>=Tk[k] else f"no (short {Tk[k]-inf_est:+.3f})"
    results[k]=(inf_est,worst[1])
    print(f"{k:>3} | {inf_est:>15.4f} | {Tk[k]:>7.4f} | {ok:>20} | {worst[1]}")

print("\nSUMMARY (Route A = single-gap avoidance tail, feeds union bound):")
covered=[k for k in range(8,14) if results[k][0]>=Tk[k]]
uncov=[k for k in range(8,14) if results[k][0]<Tk[k]]
print(f"  Route A discharges legs k = {covered}")
print(f"  remaining (need AP-minimality, opus-S131 exhaustive k<=10): k = {uncov}")
print("  => the (A') ledger splits: single-gap tail (large k) + finite AP-min (small k).")
