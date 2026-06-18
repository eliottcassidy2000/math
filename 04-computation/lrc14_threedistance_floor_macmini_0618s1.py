#!/usr/bin/env python3
"""
lrc14_threedistance_floor_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

The three-distance structure of the LRC(14) S3 residual floor.
good x  <=>  the k cluster phases {frac(e_i x)} fit in an arc of length < 5/7
        <=>  maxgap{frac(e_i x)} > 2/7.
Worst cluster = consecutive offsets E={0,1,...,k-1} => phases {0,x,2x,...,(k-1)x} (Steinhaus).

COMPUTE:
 (A) pure measure mu_k = meas{ x in [0,1): maxgap{jx: j<k} > 2/7 }, k=3..13 (NO G_P).
     Is it bounded below? Where do the good x live (rational neighborhoods a/b)?
 (B) consecutive vs random k-subset offsets at fixed k: is consecutive the WORST (min measure)?
 (C) with G_P: floor_k = min over P (|P|=13-k, P subset {1..13}) of rho*(consecutive,P).
     The true uniform floor c0.
 (D) the three-gap structure near x=j/7 (the tight points): how the integer cluster opens a gap.
"""
import sys, random, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(6181)
H=F(1,14)

def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A,h=H):
    dz=merge([iv for u in A for iv in danger_arcs(u,h)]); safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(arcs): return float(sum(b-a for a,b in arcs))

def maxgap(points):
    pts=sorted(p%1.0 for p in points); g=0.0
    for i in range(len(pts)-1): g=max(g,pts[i+1]-pts[i])
    return max(g,(pts[0]+1.0)-pts[-1])

def pure_measure(E, N=300000):
    cnt=0
    for t in range(N):
        x=(t+0.5)/N
        if maxgap([(e*x)%1.0 for e in E])>2.0/7.0: cnt+=1
    return cnt/N

def rho_with_GP(P,E,N=120000):
    GPf=[(float(a),float(b)) for a,b in safe_set(P)]
    cnt=0
    for t in range(N):
        x=(t+0.5)/N
        ins=any(a<=x<b for a,b in GPf)
        if not ins: continue
        if maxgap([(e*x)%1.0 for e in E])>2.0/7.0: cnt+=1
    return cnt/N, meas(GPf)

print("="*84)
print("(A) PURE measure mu_k = meas{x: maxgap{jx:j<k} > 2/7}, consecutive cluster")
print("="*84)
for k in range(3,14):
    E=list(range(k))
    mk=pure_measure(E)
    print(f"   k={k:2d}: mu_k = {mk:.5f}   (1/k={1/k:.4f})")

print("\n"+"="*84)
print("(B) is CONSECUTIVE the WORST shape? (min pure measure over k-subsets), k=6,8,10")
print("="*84)
for k in (6,8,10):
    cons=pure_measure(list(range(k)), N=120000)
    worst=cons; worstE=list(range(k)); samples=[]
    for _ in range(40):
        E=sorted(random.sample(range(0,3*k), k))
        if 0 not in E: E[0]=0
        m=pure_measure(E, N=80000); samples.append(m)
        if m<worst: worst=m; worstE=E
    print(f"   k={k}: consecutive mu={cons:.4f} | min over 40 random subsets={min(samples):.4f} | "
          f"overall-min={worst:.4f} at E={worstE if worstE!=list(range(k)) else 'consecutive'}")

print("\n"+"="*84)
print("(C) FLOOR with G_P: min over P (|P|=13-k) of rho*(consecutive cluster, P)")
print("="*84)
overall=(9.9,None,None)
for k in range(3,14):
    psz=13-k
    E=list(range(k))
    best=(9.9,None)
    if psz==0:
        r,_=rho_with_GP([],E);
        print(f"   k={k:2d} (P=empty): rho*={r:.5f}")
        if r<overall[0]: overall=(r,[],E)
        continue
    # try all P of size psz if small count, else sample
    Ps=list(itertools.combinations(range(1,14),psz))
    if len(Ps)>60: Ps=random.sample(Ps,60)
    for P in Ps:
        r,gp=rho_with_GP(list(P),E,N=60000)
        if r<best[0]: best=(r,P)
    print(f"   k={k:2d}: min rho* over P = {best[0]:.5f}  at P={best[1]}")
    if best[0]<overall[0]: overall=(best[0],best[1],E)
print(f"\n   OVERALL FLOOR c0 ~ {overall[0]:.5f}  at P={overall[1]}, consecutive cluster k={len(overall[2])}")
print("\nDONE.")
