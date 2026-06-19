#!/usr/bin/env python3
"""
lrc14_sector_dual_verify_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

VERIFY + STRESS-TEST the degree-R dual certificates found in lrc14_sector_lp_moments.

For each k the certificate is a vector y=(y_0..y_R) with
   meas(S7(E)) = p_0  <=  L_y(E) := sum_{r=0}^R y_r S_r(E)
PROVABLY for EVERY E, provided the dual is feasible:
   g(t) := sum_{r=0}^R y_r C(t,r)  >=  1[t=0]   for all integers t in {0,...,6}.   (*)
This is the Bonferroni/linear-programming sign condition; it is a FINITE check (7 ineqs).

So the chain "meas(S7(E)) <= cap_k for all E" reduces to:
   PROVE:  L_y(E) <= cap_k  for all valid offset sets E (0 in E, |E|=k).
We test whether CONSEC maximizes L_y(E) (then L_y(E) <= L_y(consec) <= cap_k closes it).
Since L_y is a fixed linear combination of the moments S_r, and each S_r(E) is a sum
of J(A,E) measures, "consec maximizes L_y" is the SAME shape as the seven-sector
AP-maximizer hypothesis HYP-2604 but for the SINGLE scalar L_y --- a much cleaner target.

OUTPUT:
 (1) verify (*) for each k's y over t=0..6  (PROVES the per-E bound).
 (2) verify L_y(consec)=moment-LP value and <= cap_k (PROVES consec is certified).
 (3) brute-force max of L_y(E) over bounded-spread E (the dangerous k=8,9,10,11):
     is the max == L_y(consec)?  If yes for all spreads in the finite window, the
     certificate closes that k modulo the (finite) bounded-spread reduction THM-527.
 (4) report any E with L_y(E) > cap_k  (a genuine threat) vs L_y(E) > L_y(consec)
     (an AP-beater for the SCALAR functional, harmless if still <= cap_k).
"""
import sys, itertools, random
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(7177)
H=F(1,14)

def measJ(A, E):
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
def Svec(E):
    S=[F(0)]*7; S[0]=F(1)
    for r in range(1,7):
        S[r]=sum(measJ(A,E) for A in itertools.combinations(range(1,7), r))
    return S
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(int(((e*xm)%1)*7) for e in E)
        if len(secs)==7: total+=x1-x0
    return total
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

# the certificates found (exact):
CERT = {
 8:  [F(1),F(-1),F(1),F(-9,10),F(3,5)],
 9:  [F(1),F(-13,18),F(4,9),F(-1,6)],
 10: [F(1),F(-13,18),F(4,9),F(-1,6)],
 11: [F(1),F(-1,2),F(1,6)],
 12: [F(1),F(-1,2),F(1,6)],
 13: [F(1),F(-1,2),F(1,6)],
}
def g_of_t(y,t): return sum(y[r]*F(comb(t,r)) for r in range(len(y)))
def Ly(y,S):     return sum(y[r]*S[r] for r in range(len(y)))

print("="*92)
print("(1) DUAL FEASIBILITY  g(t)=sum_r y_r C(t,r) >= 1[t=0]  for t=0..6   (PROVES per-E bound)")
print("="*92)
for k,y in CERT.items():
    gs=[g_of_t(y,t) for t in range(7)]
    feas=all(gs[t] >= (F(1) if t==0 else F(0)) for t in range(7))
    print(f"  k={k}: g(t) for t=0..6 = {[str(x) for x in gs]}   feasible={feas}")
    assert feas, f"DUAL INFEASIBLE k={k}"
print("  -> All duals feasible: meas(S7(E)) <= L_y(E) holds for EVERY offset set E. PROVED.")

print("\n"+"="*92)
print("(2) CONSEC is certified:  L_y(consec_k) <= cap_k  (exact)")
print("="*92)
for k,y in CERT.items():
    E=list(range(k)); S=Svec(E); ck=cap(k); val=Ly(y,S); s7=measS7(E)
    assert val>=s7, (k,val,s7)
    print(f"  k={k}: meas(S7)={float(s7):.5f} <= L_y(consec)={float(val):.5f} "
          f"= {val} <= cap_k={float(ck):.5f}  {'OK' if val<=ck else '*** FAIL ***'}")

print("\n"+"="*92)
print("(3) Does CONSEC maximize the scalar functional L_y(E)?  brute force, bounded spread.")
print("    If max L_y over the window == L_y(consec), certificate closes k on that window.")
print("="*92)
windows={8:13, 9:13, 10:13, 11:13}  # max element budget (spread): consec uses k-1, allow a few more
for k in [8,9,10,11]:
    y=CERT[k]; Lc=Ly(y,Svec(list(range(k)))); ck=cap(k)
    maxL=Lc; argmax='consec'; beaters_over_cap=[]; beaters_over_consec=0; nchk=0
    maxE=windows[k]
    for rest in itertools.combinations(range(1,maxE+1), k-1):
        E=[0]+list(rest); nchk+=1
        L=Ly(y,Svec(E))
        if L>Lc+F(0): beaters_over_consec+=1
        if L>maxL: maxL=L; argmax=E
        if L>ck: beaters_over_cap.append((E,L))
    print(f"  k={k}: checked {nchk} sets (spread<= {maxE}).  L_y(consec)={float(Lc):.5f}")
    print(f"        max L_y = {float(maxL):.5f} at {argmax}")
    print(f"        # E with L_y>L_y(consec): {beaters_over_consec}   "
          f"# E with L_y>cap_k (REAL threats): {len(beaters_over_cap)}")
    if beaters_over_cap:
        for E,L in beaters_over_cap[:5]:
            print(f"           THREAT E={E}  L_y={float(L):.5f} > cap={float(ck):.5f}  "
                  f"(actual meas S7={float(measS7(E)):.5f})")

print("\n"+"="*92)
print("(4) random wide-spread stress test (k=8): any E with L_y(E)>cap_8?")
print("="*92)
y=CERT[8]; ck=cap(8); worst=F(0); worstE=None; nover=0
for _ in range(400):
    sp=random.randint(8,80); E=sorted(set([0,sp]+random.sample(range(1,sp),6)))
    if len(E)!=8: continue
    L=Ly(y,Svec(E))
    if L>worst: worst=L; worstE=E
    if L>ck: nover+=1
print(f"  k=8 random: worst L_y={float(worst):.5f} at {worstE}; #over cap_8={nover}")
print(f"  cap_8={float(ck):.5f}; L_y(consec_8)={float(Ly(y,Svec(list(range(8))))):.5f}")
print("\nDONE.")
