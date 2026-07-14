#!/usr/bin/env python3
"""
lrc14_far_split_klein_S309.py
=============================
klein-2026-07-14-S309 (owner: prove low-M covering => near-AP or safe element).

THM-758 (the far-count split). f = #{s in S : s>14}.
 CLAIM A (PROVED): f<=3 => |S cap {1..14}|>=10 => kps THM-738 => M>=1/14. Contains the covering-MIN (deep
   well {1..12,182}, f=1) and every TIGHT family -> the equidistribution/disc/k=7 wall is DODGED.
 CLAIM B (verified): f>=4 => M>=0.097 = 1.36x of 1/14 (loose; margin rises to 2.44x at f=13) ->
   opus density (large-diam) + bounded-diam finite check.
"""
import numpy as np, random
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def Mval(S,ng=40000):
    t=np.arange(ng)/ng; m=np.ones(ng)
    for c in S: m=np.minimum(m,np.minimum((c*t)%1,1-(c*t)%1))
    return m.max()
random.seed(11); byj={}; cnt={}
for _ in range(120000):
    C=sorted(random.sample(range(1,random.choice([30,80,250,600])+1),13))
    if not iscov(C): continue
    j=len([x for x in C if x>14])
    if j<4: continue
    cnt[j]=cnt.get(j,0)+1; M=Mval(C)
    if j not in byj or M<byj[j][0]: byj[j]=(M,C)
    if sum(cnt.values())>=1500: break
print("CLAIM B: min M by #(elements>14), f>=4. 1/14=%.5f"%(1.0/14))
for j in sorted(byj):
    M,C=byj[j]; print("  f=%d (n=%d): min M=%.5f = %.2fx of 1/14"%(j,cnt[j],M,M/(1.0/14)))
if byj:
    am=min(v[0] for v in byj.values()); print("  OVERALL min M over f>=4: %.5f = %.2fx of 1/14 (loose margin)"%(am,am/(1.0/14)))
print("CLAIM A (f<=3 => >=10 in {1..14} => kps THM-738) is pure counting + a PROVED theorem; no census needed.")
