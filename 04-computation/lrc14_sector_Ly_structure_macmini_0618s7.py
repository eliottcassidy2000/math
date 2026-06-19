#!/usr/bin/env python3
"""
lrc14_sector_Ly_structure_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

Analyze the certificate functional  L_y(E) = sum_r y_r S_r(E)  structurally, to turn
"consec maximizes L_y" into a provable statement (or find its limits).

For the degree-2 cert (k>=11):  L_y = 1 - (1/2)S_1 + (1/6)S_2.
  S_1(E) = sum_{j=1}^6 meas{ all frac(e x) avoid sector j }  = sum_j J({j},E)
  S_2(E) = sum_{1<=j<l<=6} meas{ all avoid sectors j AND l } = sum_{j<l} J({j,l},E)

KEY simplification via the e=0 point: 0 always sits in sector 0, so sectors 1..6 are
symmetric ONLY up to the offsets.  But J({j},E) = meas{ all frac(e x) avoid arc of
length 1/7 } -- this is a SINGLE missing-arc probability.  Define for an arc I of
length L:  q_I(E) = meas{ x : frac(e x) notin I for all e in E }.  Then
  S_1 = sum over the 6 unit sectors j=1..6 of q_{sector j}.

We compute, for consec and many E:
  (a) the per-sector single-miss q_j and the pair-miss q_{jl},
  (b) whether the functional L_y is dominated by consec, and by HOW MUCH margin,
  (c) the analytic candidate bound: is S_1(E) <= S_1(consec) and S_2(E) <= S_2(consec)
      SEPARATELY? (a 2-sided structural fact) or only the combination L_y?
  (d) decompose L_y - cap_k into a positive certificate margin and report the minimum
      over a large bank (the empirical safety margin of the whole scheme).

This tells us exactly which lemma ("S_1 max at consec", "S_2 max at consec", or the
joint "L_y max at consec") is the remaining provable target.
"""
import sys, itertools, random
from math import comb, gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(70717)
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

CERT = {
 8:[F(1),F(-1),F(1),F(-9,10),F(3,5)], 9:[F(1),F(-13,18),F(4,9),F(-1,6)],
 10:[F(1),F(-13,18),F(4,9),F(-1,6)], 11:[F(1),F(-1,2),F(1,6)],
 12:[F(1),F(-1,2),F(1,6)], 13:[F(1),F(-1,2),F(1,6)],
}
def Ly(y,S): return sum(y[r]*S[r] for r in range(len(y)))

print("="*90)
print("Per-sector single-miss S_1 and pair-miss S_2 at consec, and separability test")
print("="*90)
for k in [8,11,13]:
    E=list(range(k)); S=Svec(E)
    print(f"  k={k} consec: S_1={S[1]}={float(S[1]):.5f}  S_2={S[2]}={float(S[2]):.5f}")
    # individual sector miss probs
    qj=[measJ((j,),E) for j in range(1,7)]
    print(f"         per-sector single-miss q_j (j=1..6): {[float(x) for x in qj]}")

print("\n"+"="*90)
print("Separability: is S_1(E)<=S_1(consec) and S_2(E)<=S_2(consec) for all E? (bounded spread)")
print("="*90)
for k in [11]:
    Ec=list(range(k)); Sc=Svec(Ec)
    s1c, s2c = Sc[1], Sc[2]
    s1_beats=0; s2_beats=0; Ly_beats=0; nchk=0
    worst_s1=s1c; worst_s2=s2c; worst_Ly=Ly(CERT[k],Sc)
    for rest in itertools.combinations(range(1,k+3), k-1):
        E=[0]+list(rest); S=Svec(E); nchk+=1
        if S[1]>s1c: s1_beats+=1; worst_s1=max(worst_s1,S[1])
        if S[2]>s2c: s2_beats+=1; worst_s2=max(worst_s2,S[2])
        L=Ly(CERT[k],S)
        if L>Ly(CERT[k],Sc): Ly_beats+=1; worst_Ly=max(worst_Ly,L)
    print(f"  k={k} (spread<={k+2}, {nchk} sets):")
    print(f"    S_1 beaters: {s1_beats} (worst {float(worst_s1):.5f} vs consec {float(s1c):.5f})")
    print(f"    S_2 beaters: {s2_beats} (worst {float(worst_s2):.5f} vs consec {float(s2c):.5f})")
    print(f"    L_y beaters: {Ly_beats} (worst {float(worst_Ly):.5f} vs consec {float(Ly(CERT[k],Sc)):.5f}, cap {float(cap(k)):.5f})")

print("\n"+"="*90)
print("Whole-scheme empirical safety margin: min over bank of (cap_k - L_y(E))")
print("  (if always >0, the dual bound never exceeds cap on the bank)")
print("="*90)
for k in [8,9,10,11]:
    y=CERT[k]; ck=cap(k); minmarg=ck-Ly(y,Svec(list(range(k)))); worstE='consec'
    nover=0
    # structured bank
    bank=[]
    for rest in itertools.combinations(range(1,min(k+4,16)), k-1):
        bank.append([0]+list(rest))
    # random wide
    for _ in range(150):
        sp=random.randint(k-1, 70); E=sorted(set([0,sp]+random.sample(range(1,sp),k-2)))
        if len(E)==k: bank.append(E)
    for E in bank:
        marg=ck-Ly(y,Svec(E))
        if marg<minmarg: minmarg=marg; worstE=E
        if marg<0: nover+=1
    print(f"  k={k}: min (cap - L_y) over {len(bank)} sets = {float(minmarg):.5f} at {worstE}; "
          f"# violations (L_y>cap)={nover}")
print("\nDONE.")
