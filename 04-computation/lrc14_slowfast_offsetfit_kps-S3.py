#!/usr/bin/env python3
"""
kind-pasteur-2026-06-17-S3 — the SLOW-FAST (two-scale) reduction of criterion C on S3.

For a covering 13-set S = P (small, <=13) U L (large cluster), V0=min(L), take v=max(S).
Cluster' = L\\{max}; offsets d_i = u - V0 for u in cluster'.  Write theta = V0*tau (FAST phase,
cycles ~V0 times as tau sweeps a small-part-safe arc).  At slow time tau, cluster' is jointly
level-1/14-safe iff the fast phase theta lies in  W_th(tau) := INTERSECTION_i (S_mid - d_i*tau),
S_mid=(1/14,13/14) width 6/7.  w_th(tau) = |W_th(tau)| = 6/7 - circ_width({d_i*tau mod 1}) (if>0).
PREDICTION: the widest safe arc of S\\{max} ~ max_{tau in G_P} w_th(tau) / V0.
CRITERION (v=max): W(S\\max) > 1/(7 max).  For V0 large this is ~ [exists tau in G_P: w_th(tau) > 1/7]
i.e. width({d_i*tau mod1}) < 5/7.  RHS depends ONLY on (P, offset-pattern) -> BOUNDED data.

This script: (1) validates the prediction W ~ max w_th / V0; (2) tests the offset-fit condition and
whether it is V0-independent; (3) checks it predicts C(S) via v=max.  Exact rationals.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
random.seed(7)

def safe_components(A, h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc and sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append(sc[0][1]+(1-sc[-1][0]))
    return max(ws)
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def cls(S):
    S=sorted(S); k=sum(1 for v in S if v>13)
    if k<=1: return 'S1'
    if S[-1]<13*S[0]: return 'S2'
    return 'S3'

def circ_width(points):
    """smallest arc (mod 1) containing all points; = 1 - largest gap."""
    pts=sorted(p%1 for p in points)
    if len(pts)<=1: return F(0)
    gaps=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return 1-max(gaps)

def w_theta(offsets, tau):
    """width of valid fast-phase set at slow time tau = 6/7 - circ_width({d*tau})."""
    w=F(6,7)-circ_width([d*tau for d in offsets])
    return w if w>0 else F(0)

def predict_and_actual(S):
    S=sorted(S); V=max(S); P=[u for u in S if u<=13]; L=[u for u in S if u>13]
    V0=min(L); clus2=[u for u in L if u!=V]; offs=[u-V0 for u in clus2]
    Wact=Wwidth([u for u in S if u!=V])
    thr=F(1,7*V)
    # scan tau over G_P safe arcs (sample interior points), maximize w_theta
    GP=safe_components(P) if P else [(F(0),F(1))]
    best=F(0); bestau=None
    for (a,b) in GP:
        # sample several interior taus exactly
        for t in [a+(b-a)*F(k,40) for k in range(1,40)]:
            wt=w_theta(offs,t)
            if wt>best: best=wt; bestau=t
    pred=best/V0
    return dict(V=V,V0=V0,spread=L[-1]-L[0],offs=offs,P=P,Wact=Wact,thr=thr,
                wthmax=best,pred=pred,critmax=(Wact>thr),offsetfit=(best>F(1,7)))

# gather S3 sets
S3sets=[]; tries=0
while len(S3sets)<40 and tries<300000:
    tries+=1
    nsmall=random.choice([2,3,4]); small=random.sample(range(1,14),nsmall)
    base=random.randint(20,900); spread=random.choice([14,20,28,35,45])
    nlarge=13-nsmall
    large=sorted(set(base+random.randint(0,spread) for _ in range(nlarge*3)))[:nlarge]
    if len(large)<nlarge: continue
    S=sorted(set(small+large))
    if len(S)!=13 or reduce(gcd,S)!=1 or not is_cov(S) or cls(S)!='S3': continue
    S3sets.append(S)

print("=== (1) PREDICTION W(S\\max) ~ max_tau w_th(tau)/V0  and  (3) offset-fit predicts C ===")
ok_pred=0; ok_fit=0; bad=[]
ratios=[]
for S in S3sets:
    r=predict_and_actual(S)
    rel=float(r['pred']/r['Wact']) if r['Wact']>0 else 0
    ratios.append(rel)
    fitpred = r['offsetfit']
    cmax = r['critmax']
    if fitpred==cmax: ok_fit+=1
    else: bad.append((S,r))
    # print a few
for S in S3sets[:8]:
    r=predict_and_actual(S)
    print(" V=%4d V0=%4d s=%2d  Wact=%.5f pred=%.5f (pred/act=%.2f)  w_th_max=%.4f  crit(v=max)=%s offsetfit(>1/7)=%s"
          %(r['V'],r['V0'],r['spread'],float(r['Wact']),float(r['pred']),
            float(r['pred']/r['Wact']) if r['Wact']>0 else 0, float(r['wthmax']), r['critmax'], r['offsetfit']))
import statistics
print("prediction pred/actual ratio: min=%.2f med=%.2f max=%.2f (want ~1; <1 = small-part chopping)"
      %(min(ratios),statistics.median(ratios),max(ratios)))
print("offset-fit (w_th_max>1/7) matches crit(v=max):",ok_fit,"/",len(S3sets))
if bad:
    print("MISMATCHES:")
    for S,r in bad[:5]:
        print("   S=",S," wthmax=",float(r['wthmax'])," crit=",r['critmax']," fit=",r['offsetfit'])

print()
print("=== (2) V0-INDEPENDENCE: fix (P, offset-pattern), vary V0; does offset-fit stay constant? ===")
# take one offset pattern + small part, embed at growing V0 (keep covering by choosing base)
P=[1,2,3]; offs_pat=[0,4,6,9,13,15,22,24,27]  # 9 cluster offsets -> |L|=10, +3 small =13
for V0 in [200,400,700,1100,2000]:
    L=[V0+d for d in [0]+offs_pat]  # careful: offsets here from V0, include 0; need covering
    S=sorted(set(P+L))
    if len(S)!=13:
        print("  V0=%d skip (size %d)"%(V0,len(S))); continue
    cov=is_cov(S)
    r=predict_and_actual(S)
    print("  V0=%4d covering=%s  w_th_max=%.4f  offsetfit=%s  crit(v=max)=%s  Wact*7V=%.2f"
          %(V0,cov,float(r['wthmax']),r['offsetfit'],r['critmax'],float(r['Wact']*7*r['V'])))
print("(w_th_max should be ~constant in V0 = the bounded-data reduction; crit should stabilize)")
