#!/usr/bin/env python3
"""
lrc14_sector_fast_Jall_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

FAST joint-occupancy engine: compute ALL S_r(E) and the certificate functional L_y(E)
in ONE breakpoint pass over x in [0,1).  At each elementary interval the set of OCCUPIED
sectors occ subset {0..6} is constant; that interval contributes its length to J(A,E)
for every A subset {1..6} with A disjoint from occ.  Summing C(.,.) we get S_r directly:
  S_r += len * (number of r-subsets of ({1..6} \ occ))  = len * C(6 - |occ cap {1..6}|, r),
but careful: sector 0 always occupied (e=0) so occ always contains 0; A subset {1..6}.
  free = {1..6} \ occ ;  contributes to J(A) iff A subset free  =>  S_r += len*C(|free|, r).
So per interval we only need |free| = 6 - (#occupied among 1..6).  VERY fast & exact.

Then L_y(E)=sum_r y_r S_r and meas(S7)=sum over intervals with |free|=0 of len (since
S7 = all sectors hit = no free sector).  All exact rational.

This replaces the per-A breakpoint engine (63x speedup).  We re-verify against the slow
engine on a few shapes, then run the heavy maximizer sweeps that were too slow before:
  - k=8,9,10,11: does consec maximize L_y over a LARGE bounded-spread window?
  - any E (bounded spread or random wide) with L_y(E) > cap_k?
"""
import sys, itertools, random
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(771717)
H=F(1,14)

def Sfree_and_S7(E):
    """one pass: returns (S_0..S_6 list, meas S7) exactly."""
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    S=[F(0)]*7; s7=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        L=x1-x0; xm=(x0+x1)/2
        occ=set(int(((e*xm)%1)*7) for e in E)   # includes 0
        free=6-len(occ & {1,2,3,4,5,6})
        for r in range(7):
            S[r]+=L*F(comb(free,r))
        if free==0: s7+=L
    return S, s7

# ---- slow reference for cross-check ----
def measJ_slow(A,E):
    Aset=set(A); E=sorted(set(E))
    if 0 in Aset: return F(0)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(int(((e*xm)%1)*7) not in Aset for e in E): total+=x1-x0
    return total
def Svec_slow(E):
    S=[F(0)]*7;S[0]=F(1)
    for r in range(1,7): S[r]=sum(measJ_slow(A,E) for A in itertools.combinations(range(1,7),r))
    return S

def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgmerge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    if not P: return F(1)
    dz=mgmerge([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k):
    psz=13-k
    if psz==0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1,14),psz))

CERT={8:[F(1),F(-1),F(1),F(-9,10),F(3,5)],9:[F(1),F(-13,18),F(4,9),F(-1,6)],
 10:[F(1),F(-13,18),F(4,9),F(-1,6)],11:[F(1),F(-1,2),F(1,6)],
 12:[F(1),F(-1,2),F(1,6)],13:[F(1),F(-1,2),F(1,6)]}
def Ly(y,S): return sum(y[r]*S[r] for r in range(len(y)))

# cross-check fast vs slow
print("cross-check fast Sfree engine vs slow per-A engine:")
for E in [list(range(8)),[0,1,3,7,15,31,63,127],[0,2,3,4,5,6,7,9],list(range(11))]:
    Sf,s7=Sfree_and_S7(E); Ss=Svec_slow(E)
    print(f"  E(k={len(E)}): S match={Sf==Ss}  s7={float(s7):.5f}")
print()

print("="*94)
print("CONSEC maximizes L_y?  LARGE bounded-spread sweep + does any E exceed cap_k?")
print("="*94)
windows={8:16,9:15,10:14,11:14,12:14,13:13}
for k in range(8,14):
    y=CERT[k]; ck=cap(k); Lc=Ly(y,Sfree_and_S7(list(range(k)))[0])
    maxL=Lc; argmax='consec'; over_cap=[]; over_consec=0; nchk=0
    me=windows[k]
    for rest in itertools.combinations(range(1,me+1), k-1):
        E=[0]+list(rest); S,_=Sfree_and_S7(E); nchk+=1
        L=Ly(y,S)
        if L>Lc: over_consec+=1
        if L>maxL: maxL=L; argmax=E
        if L>ck: over_cap.append((E,float(L)))
    status = "consec is MAX" if argmax=='consec' else f"NOT consec (max at {argmax})"
    print(f"  k={k}: {nchk} sets spread<= {me}.  L_y(consec)={float(Lc):.5f} cap={float(ck):.5f}")
    print(f"        max L_y={float(maxL):.5f}  [{status}]  #over consec={over_consec}  #over cap={len(over_cap)}")
    for E,L in over_cap[:4]:
        print(f"          *** L_y({E})={L:.5f} > cap_{k}={float(ck):.5f}  (meas S7={float(Sfree_and_S7(E)[1]):.5f}) ***")

print("\n"+"="*94)
print("random WIDE-spread stress (each k, 500 sets): max L_y vs cap_k")
print("="*94)
for k in range(8,14):
    y=CERT[k]; ck=cap(k); Lc=Ly(y,Sfree_and_S7(list(range(k)))[0])
    worst=Lc; worstE='consec'; nover=0
    for _ in range(500):
        sp=random.randint(k-1, 120); pool=random.sample(range(1,sp), k-2)
        E=sorted(set([0,sp]+pool))
        if len(E)!=k: continue
        L=Ly(y,Sfree_and_S7(E)[0])
        if L>worst: worst=L; worstE=E
        if L>ck: nover+=1
    print(f"  k={k}: worst L_y={float(worst):.5f} (cap={float(ck):.5f}, consec={float(Lc):.5f})  "
          f"#over cap={nover}  {'at '+str(worstE) if worstE!='consec' else ''}")
print("\nDONE.")
