"""
UNIFORM BOUND chase: drop j binding units from {1..13}, add killers (mults of 14 covering the dropped q).
Does M -> 1/14 (displacement->0, razor-thin, NO margin) or stay bounded (margin)?
"""
from math import lcm
import numpy as np
def M_grid(S,Q):
    t=np.arange(1,Q)/Q; f=np.ones(Q-1)
    for v in S: f=np.minimum(f,np.abs(((v*t+0.5)%1)-0.5))
    return f.max()
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
ap=list(range(1,14)); units=[13,11,9,5,3,1]  # drop in this order (13 first = tightest)
print(f"1/14 = {1/14:.6f}\n")
print("drop j units (13,11,9,...), replace each missing q by a killer (mult of 14):")
for j in range(0,6):
    dropped=units[:j]
    rem=[x for x in ap if x not in dropped]
    # need killers for q in dropped (each dropped unit q needs a mult of lcm(q,14)); use distinct killers
    killers=[]
    for q in dropped:
        w=lcm(q,14)
        while w in rem or w in killers: w+=lcm(q,14)
        killers.append(w)
    S=rem+killers
    if len(set(S))!=13: 
        print(f"  j={j}: size {len(set(S))} (skip)"); continue
    Q=max(400000, max(S)*15)
    M=M_grid(S,Q); cov=is_cov(S)
    print(f"  j={j}: drop {dropped}, killers {killers} -> M={M:.6f} (delta={M-1/14:+.6f}) covering={cov}")
print()
# also: the divisor-loaded the OTHER way (drop slack 12, killer 84*m) for contrast
print("contrast -- drop slack 12, killer 84*m (divisor-loaded, M INCREASES away from 1/14):")
for m in [1,2,5]:
    S=[1,2,3,4,5,6,7,8,9,10,11,13,84*m]
    print(f"  m={m}: M={M_grid(S,max(400000,84*m*15)):.6f}")
print()
# the tightest single-drop: is there a drop+killer giving M even closer to 1/14?
print("scan ALL single drops k=2..13 with smallest killer, find min M:")
best=(1.0,None)
for k in range(2,14):
    w=lcm(k,14)
    S=[x for x in ap if x!=k]+[w]
    if len(set(S))==13 and is_cov(S):
        M=M_grid(S,max(400000,w*15))
        if M<best[0]: best=(M,(k,w))
print(f"  min single-drop M = {best[0]:.6f} at drop {best[1][0]}, killer {best[1][1]} (delta={best[0]-1/14:+.6f})")
