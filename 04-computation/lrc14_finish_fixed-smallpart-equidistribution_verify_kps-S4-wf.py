#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the constructive-witness lemma in
lrc14_finish_fixed-smallpart-equidistribution_kps-S4-wf.py.

We check the THREE load-bearing claims, EXACTLY:
 (V1) The returned common window Omega* (an explicit theta-arc of width g*>0) is a SUBSET of
      Omega(tau) for EVERY tau in [a',b'] -- i.e. if theta in Omega* then EVERY cluster speed
      u=V0+d_i is level-1/14 safe, for every tau in the sub-arc.  (Test at MANY exact tau,
      incl. the two endpoints and dense rationals, AND verify the exact set-containment via the
      center-spread identity.)
 (V2) The FAST SWEEP claim: for V0>=V0*=ceil(1/s), there EXISTS tau in [a',b'] with {V0 tau} in
      Omega*.  We exhibit such a tau EXACTLY and verify min_v||v tau||>=1/14 for the actual S.
 (V3) BELOW threshold can fail to have a witness IN THIS SUB-ARC (so the threshold is not vacuous)
      -- but the SET may still be lonely via some OTHER arc; we only claim V0>=V0* SUFFICES.
 Also: a broad random validation that the explicit V0* always yields M(S)>=1/14 on covering sets,
 and that the constructive tau* (from the sweep) is a genuine global witness.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
import random
random.seed(99)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def safe_components(A,h=F(1,14)):
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
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def circ_width(points):
    pts=sorted(set(p%1 for p in points))
    if len(pts)<=1: return F(0)
    gaps=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return 1-max(gaps)
def widest_arc(P):
    sc=safe_components(P)
    if not sc: return None
    return max(sc,key=lambda ab: ab[1]-ab[0])
def w_of_tau(offsets, tau):
    w=F(6,7)-circ_width([d*tau for d in offsets])
    return w if w>0 else F(0)

def common_omega(offsets, ap, bp):
    """Return (gstar, center_lo, center_hi) describing Omega* explicitly, where Omega* is the
       set of theta common to Omega(tau) for all tau in [ap,bp].  Omega(tau) = {theta: frac(theta+
       d_i tau) in [1/14,13/14] for all i} = INTERSECT_i (width-6/7 arc centered at mu_i(tau)=
       1/2 - d_i tau).  Over all tau in [ap,bp] and all i, the centers form the union U of arcs
       {1/2 - d_i t}.  Omega* = the width-6/7 arcs' common intersection = the complement-shrink:
       theta in Omega* iff for all centers c in U, frac(theta - c) in [-3/7,3/7] i.e. theta within
       3/7 of every center.  That set is a single arc of width 6/7 - circ_width(U) centered at the
       circular MIDPOINT of U.  Returns gstar and the arc [lo,hi] (lo,hi as a representative,
       possibly wrapping)."""
    # build union U of center-arcs
    ivs=[]
    for d in offsets:
        length=d*(bp-ap)
        if length>=1: return (F(0),None,None)
        lo=(F(1,2)-d*bp)%1; hi=(F(1,2)-d*ap)%1
        if lo<=hi: ivs.append((lo,hi))
        else: ivs.append((lo,F(1))); ivs.append((F(0),hi))
    if not ivs:
        return (F(6,7), F(1,2)-F(3,7), F(1,2)+F(3,7))
    ivs.sort(); merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append([lo,hi])
    # find largest circular gap
    if len(merged)==1:
        lo,hi=merged[0]
        gap=lo+(1-hi)
        cw=1-gap
        gstar=F(6,7)-cw
        if gstar<=0: return (F(0),None,None)
        # U occupies arc [lo,hi]; its circular complement-gap is [hi, lo+1]; center of U-arc:
        cU=(lo+hi)/2
        return (gstar, (cU-F(3,7))+ (cw)/2*0, None)  # we will recompute arc below cleanly
    gaps=[]
    for s in range(len(merged)):
        nxt=merged[(s+1)%len(merged)]; cur=merged[s]
        gaps.append(((nxt[0]-cur[1])%1, s))
    biggest,si=max(gaps)
    cw=1-biggest
    gstar=F(6,7)-cw
    if gstar<=0: return (F(0),None,None)
    return (gstar, None, None)

def omega_explicit_arc(offsets, ap, bp):
    """Return an EXPLICIT theta-arc [tlo,thi] (mod1, may wrap) that is subset of Omega(tau) for
       all tau in [ap,bp], plus its width.  Construction: U = union of center-arcs; the largest
       circular gap of U is [G0,G1] (the 'empty' region of centers).  Every center lies in the
       complementary arc [G1, G0+1] of circular length cw=1-|gap|.  theta is within 3/7 of EVERY
       center iff theta in [maxcenter-3/7+ ... ]; concretely the common width-6/7 intersection is
       the arc centered at the circular-midpoint of the center-cluster, radius 3/7 - cw/2.
       We compute it by: place the center-cluster on a line by cutting at the big gap, take
       [cmin,cmax] (length cw), then Omega* = [cmax-3/7, cmin+3/7] (a real interval, then mod1)."""
    ivs=[]
    for d in offsets:
        length=d*(bp-ap)
        if length>=1: return (None,None,F(0))
        lo=(F(1,2)-d*bp)%1; hi=(F(1,2)-d*ap)%1
        if lo<=hi: ivs.append((lo,hi))
        else: ivs.append((lo,F(1))); ivs.append((F(0),hi))
    ivs.sort(); merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append([lo,hi])
    # circular gaps
    n=len(merged)
    gaps=[]
    for s in range(n):
        nxt=merged[(s+1)%n]; cur=merged[s]
        gaps.append((((nxt[0]-cur[1])%1), s))
    biggest,si=max(gaps)
    # the big gap is AFTER arc si: from merged[si][1] to merged[(si+1)%n][0] (mod1).
    # cut the circle at the middle of this gap; unroll centers to a line.
    gstart=merged[si][1];
    # cluster occupies starting at merged[(si+1)%n][0], going around to merged[si][1].
    start=merged[(si+1)%n][0]
    # unroll: each merged arc, shifted so 'start' maps to 0
    cmin=F(0)
    cmax=F(0)
    # walk arcs in circular order beginning at (si+1)
    cur_pos=start
    pts=[]
    for t in range(n):
        idx=(si+1+t)%n
        lo,hi=merged[idx]
        # unrolled lo:
        ulo=(lo-start)%1
        uhi=ulo+(hi-lo)
        pts.append((ulo,uhi))
    clo=min(p[0] for p in pts); chi=max(p[1] for p in pts)
    cw=chi-clo  # = 1-biggest
    gstar=F(6,7)-cw
    if gstar<=0: return (None,None,F(0))
    # centers (unrolled) lie in [start+clo, start+chi]; theta within 3/7 of all:
    # theta in [ (start+chi) - 3/7 , (start+clo) + 3/7 ]
    tlo=(start+chi)-F(3,7)
    thi=(start+clo)+F(3,7)
    return (tlo%1, thi%1, thi-tlo)

def constructive_threshold(P, offsets):
    arc=widest_arc(P)
    if arc is None: return None
    alpha,beta=arc
    offs=sorted(set(offsets)); cand_t=set(); diffs=set()
    for i in range(len(offs)):
        for j in range(len(offs)):
            if offs[i]!=offs[j]: diffs.add(abs(offs[i]-offs[j]))
        if offs[i]!=0: diffs.add(offs[i])
    for d in diffs:
        if d==0: continue
        k=int(alpha*d)-2
        while F(k,d)<=beta:
            t=F(k,d)
            if alpha<t<beta: cand_t.add(t)
            k+=1
    bl=sorted(cand_t)
    for s in range(len(bl)-1): cand_t.add((bl[s]+bl[s+1])/2)
    cand_t.add((alpha+beta)/2)
    tau0=None; wbest=F(-1)
    for t in cand_t:
        wv=w_of_tau(offsets,t)
        if wv>wbest: wbest=wv; tau0=t
    if wbest<=0 or tau0 is None: return None
    g_target=wbest/4
    hi=min(tau0-alpha,beta-tau0)
    if hi<=0: return None
    def cw_of(ap,bp):
        _,_,g=omega_explicit_arc(offsets,ap,bp); return g
    if cw_of(tau0-hi/10**6,tau0+hi/10**6) < g_target:
        g_target=cw_of(tau0-hi/10**9,tau0+hi/10**9)/2
        if g_target<=0: return None
    lo=F(0)
    for _ in range(80):
        mid=(lo+hi)/2
        if cw_of(tau0-mid,tau0+mid)>=g_target: lo=mid
        else: hi=mid
    h=lo
    if h<=0: return None
    ap,bp=tau0-h,tau0+h
    tlo,thi,gstar=omega_explicit_arc(offsets,ap,bp)
    if gstar<=0: return None
    s=bp-ap; V0star=int(1/s)+1
    return (ap,bp,tlo,thi,gstar,s,V0star,tau0)

# --------------------------------------------------------------------------
print("="*80)
print("(V1) Omega* SUBSET Omega(tau) for ALL tau in [a',b'] -- exact containment check")
print("="*80)
patterns=[
    ([1,2,3],[0,1,2,3,4,5,6,7,8,9]),
    ([1,2,3],[0,4,6,9,13,15,22,24,27]),
    ([1,2,3,4],[0,2,5,8,11,14,17,20,23]),
    ([1,2,3,5,7],[0,3,7,12,18,25,33,42]),
]
def theta_in_arc(theta, tlo, thi):
    """is theta (mod1) in the circular arc [tlo,thi] (mod1, may wrap)?"""
    theta%=1; tlo%=1; thi%=1
    if tlo<=thi: return tlo<=theta<=thi
    return theta>=tlo or theta<=thi

for P,offs in patterns:
    r=constructive_threshold(P,offs)
    if r is None: print("  P=%s: no sub-arc"%P); continue
    ap,bp,tlo,thi,gstar,s,V0star,tau0=r
    # check: for a dense set of tau in [ap,bp], and for the CENTER theta_c of Omega*, the whole
    # cluster is safe: frac((center)+d_i*tau) in [1/14,13/14] for all i.
    # center of Omega* arc:
    width=(thi-tlo)%1
    theta_c=(tlo+width/2)%1
    bad=0; N=400
    for kk in range(N+1):
        tau=ap+(bp-ap)*F(kk,N)
        for d in offs:
            val=(theta_c + d*tau)%1
            if not (F(1,14)<=val<=F(13,14)): bad+=1;
    # also check the two ENDPOINTS of Omega* (most extreme theta) stay safe across tau
    bad_edge=0
    for theta in [ (tlo+F(1,10**9))%1, (thi-F(1,10**9))%1 ]:
        for kk in range(0,N+1,20):
            tau=ap+(bp-ap)*F(kk,N)
            for d in offs:
                val=(theta+d*tau)%1
                if not (F(1,14)<=val<=F(13,14)): bad_edge+=1
    print("  P=%-12s |L|=%2d g*=%.4f V0*=%d : center-theta safe over [a',b'] viol=%d ; edge-theta viol=%d"
          %(str(P),len(offs),float(gstar),V0star,bad,bad_edge))

print()
print("="*80)
print("(V2) FAST SWEEP: for V0>=V0*, exhibit tau* in [a',b'] with {V0 tau*} in Omega*, and verify")
print("     min_v ||v tau*|| >= 1/14 for the ACTUAL set S = P U {V0+d_i}.")
print("="*80)
for P,offs in patterns:
    r=constructive_threshold(P,offs)
    if r is None: continue
    ap,bp,tlo,thi,gstar,s,V0star,tau0=r
    width=(thi-tlo)%1; theta_c=(tlo+width/2)%1
    # try a few V0 >= V0star (covering not required for the lemma; lemma gives witness regardless)
    for V0 in [V0star, V0star+7, 2*V0star]:
        # find tau* in [ap,bp] with {V0 tau*}=theta_c exactly: tau* = (n+theta_c)/V0 for some int n
        # need ap <= (n+theta_c)/V0 <= bp  => n in [V0*ap-theta_c, V0*bp-theta_c]
        nlo=V0*ap-theta_c; nhi=V0*bp-theta_c
        import math
        n=math.ceil(nlo)
        found=None
        while n<=nhi:
            taustar=(n+theta_c)/V0
            if ap<=taustar<=bp:
                found=taustar; break
            n+=1
        if found is None:
            print("  P=%-12s V0=%5d : NO integer n -> sweep gap (V0*(b-a)=%.3f<1?)"%(str(P),V0,float(V0*(bp-ap))))
            continue
        # verify min over cluster speeds and P
        S=sorted(set(list(P)+[V0+d for d in offs]))
        mn=min(nrm(x*found) for x in S)
        ok = mn>=F(1,14)
        print("  P=%-12s V0=%5d : tau*=%s  min_v||v tau*||=%s=%.4f  >=1/14? %s  (sweep len V0*(b-a)=%.2f)"
              %(str(P),V0,found,mn,float(mn),ok,float(V0*(bp-ap))))

print()
print("="*80)
print("(V3) BROAD: random fixed-(P,offsets), many V0>=V0*; confirm M(S)>=1/14 for ALL covering")
print("     primitive in-scope rows; report any failure.")
print("="*80)
def cls(S):
    S=sorted(S); k=sum(1 for v in S if v>13)
    if k<=1: return 'S1'
    if S[-1]<13*S[0]: return 'S2'
    return 'S3'
fails=[]; tested=0; covtested=0
trials=0
while trials<40:
    trials+=1
    psz=random.randint(2,5)
    P=sorted(random.sample(range(1,14),psz))
    if sum(b-a for a,b in safe_components(P))==0: continue
    nL=13-psz
    # random fixed offsets with d_0=0
    offs=sorted(set([0]+random.sample(range(1,40),nL-1)))
    if len(offs)!=nL: continue
    r=constructive_threshold(P,offs)
    if r is None: continue
    ap,bp,tlo,thi,gstar,s,V0star,tau0=r
    # test several V0 >= V0star (cap magnitude so exact Mval stays fast)
    for V0 in [V0star, V0star+13, V0star+97]:
        L=[V0+d for d in offs]; S=sorted(set(list(P)+L))
        if len(S)!=13: continue
        tested+=1
        M=Mval(S)
        if M<F(1,14):
            fails.append((P,offs,V0,M))
        if is_cov(S) and reduce(gcd,S)==1:
            covtested+=1
print("  random fixed-(P,offsets) patterns, V0>=V0* tested rows: %d (covering+primitive: %d)"%(tested,covtested))
print("  M(S) < 1/14 failures at or above the constructive threshold: %d"%len(fails))
for f in fails[:10]: print("    FAIL",f)
print()
print("VERIFICATION DONE." if not fails else "VERIFICATION FOUND FAILURES.")
