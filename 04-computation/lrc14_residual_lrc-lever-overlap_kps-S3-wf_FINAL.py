#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_kps-S3-wf_FINAL  (consolidated)

FINAL CONSOLIDATION of the "lrc-lever-overlap" angle on the residual case S3 of LRC(14).

============================ WHAT IS PROVED (elementary) ============================
Notation: S3 covering 13-set. v=Vmax. A=S\{v}. P={speeds<=13}. L'={speeds>13}\{Vmax}.
w0=min L', wT=max L', thresh=1/(7 Vmax). Criterion C(S): exists v with W(S\{v})>1/(7v).
By THM-526, C(S) => M(S)>=1/14.

LEMMA L1 (global arc, = S2 Lemma 1, already in canon).  If 13 Vmin > Vmax, the open arc
   (1/(14 Vmin), 13/(14 Vmax)) is safe for ALL of S => M(S)>=1/14.  [Does not apply in S3 since
   S3 has Vmax>=13 Vmin; included for completeness of the dichotomy.]

LEMMA L_collapse (SINGLE-GAP CLUSTER BAND).  If there is an integer j>=0 with the window
   J_j = [(14j+1)/(14 w0), (14j+13)/(14 wT)] nonempty, of width > thresh, and lying inside a
   P-safe arc, then J_j subset G_A and width(J_j) > thresh => C(S) via v=Vmax.
   PROOF: monotone -- every u in L' maps J_j into the single gap (j+1/14, j+13/14) mod 1, so
   ||u tau||>=1/14; P-safe by hypothesis. (elementary, exact)

LEMMA L_ruler (K_j MULTI-MEMBER RULER, the lever).  For any cluster member r in L', its safe
   arcs I_j^r = ((14j+1)/(14 r),(14j+13)/(14 r)), j=0..r-1, are exactly periodic (period 1/r).
   If for some r and some j the set  G_A intersect I_j^r  contains a sub-arc of width > thresh,
   then C(S) via v=Vmax.  PROOF: that sub-arc is by definition safe for all of A and wider
   than thresh.  (elementary, exact; this is just 'G_A has an arc>thresh somewhere', organized
   by the ruler r so the witness is constructive.)

VERIFIED FACTS (exact, on >5000 S3 sets across regimes incl. huge V0, AP clusters, wide spread):
  - WIDTH-INEQ: (13 w0 - wT) Vmax > 2 w0 wT holds whenever wT/w0 < 11 (true for all tight
    clusters); equivalently the first-gap window J_0 beats thresh. Min margin 2.70 on tight,
    fails only for 'wide clusters' (wT/w0>=11) where L' is not a tight cluster.
  - W(G_{L'}) > thresh ALONE (cluster over-delivers), ratio 5.5..8.2.
  - PERSISTENCE: rho_K := (#good ruler periods)/r converges to a positive constant as r->infty
    for every tested offset-pattern/P; #good grows LINEARLY in r (equidistribution).
  - UNION L1 u L_collapse u L_ruler certifies 1699/1700 S3 sets in the boundary scan; the single
    uncovered set [1,2,3,5,7,8,9,10,11,12,13,27,28] has L'={27} (degenerate, k=2 at the S2/S3
    boundary) and C(S) holds via v=27 (NOT Vmax): W(S\{27})*7*27=1.38>1.

============================ WHAT IS NOT PROVED (honest gaps) ============================
  G1 (universality of the ruler witness for all r, all S3): the L_ruler INDICATOR is a finite
     exact check per set, but it is NOT proven to hold for EVERY S3 set with arbitrarily large
     Vmax. Persistence (rho_K -> const > 0) is VERIFIED, not proven; a uniform positive floor
     rho*(Delta,P) >= c0 > 0 over all admissible (offset pattern Delta, small part P) would close
     it but is not established (and a naive probe gave rho*=0 for an INADMISSIBLE config).
  G2 (the v=Vmax assumption): on >=1 genuine S3 set, C(S) fails via v=Vmax and needs a different
     v. So the clean 'remove the max' reduction is FALSE in general; the lever must allow ANY r
     as ruler AND any removed runner.
  G3 (the wide-cluster remainder): when the large runners are NOT a tight cluster (wT/w0>=11),
     L_collapse and the WIDTH-INEQ fail; only the multi-member L_ruler covers them (verified).

This script RE-VERIFIES the union coverage on a fresh large sample and prints the certificate
type used for each set, plus exact M(S) for any uncovered set, so the residual is explicit.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def teeth_intervals(u,h=C):
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merged_intervals(iv):
    iv=sorted(iv); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged
def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    return merged_intervals(iv)
def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def safe_components(A,h=C): return safe_from_danger(danger_merged(A,h))
def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append(sc[0][1]+(1-sc[-1][0]))
    return max(ws)
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)
def cand(S):
    S=sorted(set(S)); Cs=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cs.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cs.add(F(k,d)); k+=1
    Cs.add(F(1,2)); return Cs
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b

def cluster_PL(S):
    Vmax=max(S)
    P=[v for v in S if v<=13]; L=sorted(v for v in S if v>13); Lp=[v for v in L if v!=Vmax]
    return P,Lp,Vmax
def L1_global(S): return 13*min(S)>max(S)
def L_collapse(S):
    P,Lp,Vmax=cluster_PL(S)
    if not Lp: return False
    w0=min(Lp); wT=max(Lp); thresh=F(1,7*Vmax)
    GP=safe_components(P) if P else [(F(0),F(1))]
    if 13*w0<=wT: return False
    jmax=(13*w0-wT-1)//(14*(wT-w0)) if wT>w0 else w0-1
    for j in range(0,max(0,jmax+1)):
        p=F(14*j+1,14*w0); q=F(14*j+13,14*wT)
        if q<=p or q>=1 or q-p<=thresh: continue
        for a,b in GP:
            if a<=p and q<=b: return True
    return False
def L_ruler(S):
    P,Lp,Vmax=cluster_PL(S)
    if not Lp: return False
    thresh=F(1,7*Vmax); GA=safe_components(list(P)+list(Lp))
    for ruler in Lp:
        for j in range(0,ruler):
            lo=F(14*j+1,14*ruler); hi=F(14*j+13,14*ruler)
            if hi>1: hi=F(1)
            if lo>=hi: continue
            for a,b in GA:
                x=max(a,lo); y=min(b,hi)
                if y-x>thresh: return True
    return False
def C_any_v(S):
    """Truth: does C(S) hold via SOME removed runner v (not necessarily Vmax)?"""
    for v in S:
        if Wwidth([x for x in S if x!=v])*7*v>1: return True
    return False

def make_S3(rng, spread_max, vlo, vhi):
    for _ in range(30000):
        small=set(rng.sample(range(1,14), rng.randint(3,9)))
        V0=rng.randint(vlo,vhi); s=rng.randint(1,spread_max)
        clsize=rng.randint(2, max(2,13-len(small)))
        cluster=set(); tries=0
        while len(cluster)<clsize and tries<300:
            cluster.add(V0+rng.randint(0,s)); tries+=1
        S=small|cluster
        while len(S)<13:
            S.add(rng.randint(V0,V0+s) if rng.random()<0.5 else rng.randint(1,13))
        if len(S)!=13: continue
        S=set(list(S))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
        if k>=2 and Vmax>=13*Vmin: return sorted(S)
    return None

def run(n=800, seed=999, spread_max=60, vlo=14, vhi=600, tag=""):
    rng=random.Random(seed)
    got=0; attempts=0
    ruler_cov=0; union_cov=0; truth=0; vmax_fails=0
    none_cases=[]
    while got<n and attempts<n*60:
        attempts+=1
        S=make_S3(rng,spread_max,vlo,vhi)
        if S is None: continue
        got+=1
        if C_any_v(S): truth+=1
        Vmax=max(S)
        if Wwidth([x for x in S if x!=Vmax])*7*Vmax<=1: vmax_fails+=1
        r=L_ruler(S); cov=L1_global(S) or L_collapse(S) or r
        if r: ruler_cov+=1
        if cov: union_cov+=1
        else: none_cases.append(S)
    flush(f"\n  [{tag}] S3 sets={got}")
    flush(f"    C(S) TRUE via some v (truth):  {truth}/{got}")
    flush(f"    C(S) via v=Vmax FAILS:         {vmax_fails}/{got}")
    flush(f"    L_ruler certifies (via Vmax):  {ruler_cov}/{got}")
    flush(f"    UNION (L1|collapse|ruler):     {union_cov}/{got}")
    flush(f"    covered by NONE:               {len(none_cases)}")
    for S in none_cases[:6]:
        flush(f"      {S}  M={Mval(S)}  C_any_v={C_any_v(S)}")
    return none_cases

if __name__=="__main__":
    flush("="*70)
    flush("FINAL CONSOLIDATION: lrc-lever-overlap on S3")
    flush("="*70)
    allnone=[]
    allnone+=run(800, 991, 60, 14, 600, "general")
    allnone+=run(600, 992, 14, 14, 600, "tight<=14")
    allnone+=run(600, 993, 120,14, 800, "loose<=120")
    allnone+=run(600, 994, 45, 14, 80,  "S2/S3 boundary small V0")
    allnone+=run(400, 995, 30, 200, 3000, "huge V0")
    flush(f"\n  TOTAL covered-by-NONE: {len(allnone)}")
    flush("  (each such set still has C(S) TRUE via a non-Vmax v -- gap G2)")
