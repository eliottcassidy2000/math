#!/usr/bin/env python3
"""
lrc14_rho_star_scan_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

Scan the limiting density  rho*(Delta,P) = meas{ x in G_P : maxgap{frac(e_i x)} > 2/7 }
(validated reformulation of kind-pasteur's good-period density; THM-526/OPEN-Q-108).
e_i = co-offsets of the cluster (s - d_i, s=spread); k=|cluster|; P subset {1..13}, |P|=13-k.

QUESTIONS:
 (1) What is inf rho* over admissible (P,Delta)? Bounded away from 0, or ->0 (tight)?
 (2) The canonical worst case (canon): P={1,2,3}, k=10 dense cluster, tau*=j/7.
 (3) Which (P,Delta) give rho*=0 (cluster fully covers G_P) -- are they NON-covering?
 (4) Anatomy of the minimizer: where is x bad, what does the cluster do there.
"""
import sys, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(618)
H = F(1,14)

def danger_arcs(u, h=H):
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
def safe_set(A, h=H):
    dz=merge([iv for u in A for iv in danger_arcs(u,h)])
    safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(arcs): return float(sum(b-a for a,b in arcs))

def maxgap_float(points):
    pts=sorted(p%1.0 for p in points)
    g=0.0
    for i in range(len(pts)-1): g=max(g, pts[i+1]-pts[i])
    g=max(g, (pts[0]+1.0)-pts[-1])
    return g

def rho_star(P, Delta, N=60000):
    s=max(Delta); e=[s-d for d in Delta]
    GPf=[(float(a),float(b)) for a,b in safe_set(P)]
    cnt=0; goodm=0
    thr=2.0/7.0
    for t in range(N):
        x=(t+0.5)/N
        ins=False
        for a,b in GPf:
            if a<=x<b: ins=True; break
        if not ins: continue
        if maxgap_float([(ei*x)%1.0 for ei in e]) > thr: cnt+=1
    return cnt/N, meas(GPf)

print("="*86)
print("(2) CANONICAL WORST: P={1,2,3}, k=10 dense clusters")
print("="*86)
P3=[1,2,3]
for Delta in [list(range(10)), [0,1,2,3,4,5,6,7,8,9], [0,1,2,3,4,5,6,7,8,18],
              sorted(random.sample(range(0,20),10)) , sorted(random.sample(range(0,40),10))]:
    Delta=sorted(set(Delta))
    if len(Delta)<2: continue
    r,gp=rho_star(P3,Delta)
    print(f"   P={P3} Delta(size {len(Delta)})={Delta}: rho*={r:.5f}  meas(G_P)={gp:.4f}")

print("\n"+"="*86)
print("(1)(3) BROAD SCAN: min rho* over admissible (P, dense Delta); flag rho*~0")
print("="*86)
best=(9.9,None,None); zeros=[]
results=[]
for trial in range(700):
    k=random.randint(3,11)          # cluster size; |P|=13-k in [2,10]
    psz=13-k
    if psz<1 or psz>13: continue
    P=sorted(random.sample(range(1,14), psz))
    # dense-ish cluster offsets (spread modest so it's a tight cluster)
    spread=random.choice([k-1,k,k+2,2*k,3*k,5*k])
    Delta=sorted(set([0]+random.sample(range(1,spread+1), min(k-1,spread))))
    if len(Delta)!=k:
        # pad
        extra=[d for d in range(1,spread+2) if d not in Delta]
        random.shuffle(extra)
        Delta=sorted(set(Delta+extra[:k-len(Delta)]))
    if len(Delta)<3: continue
    r,gp=rho_star(P,Delta, N=40000)
    results.append((r,gp,k,P,Delta))
    if r<best[0]: best=(r,P,Delta)
    if r<1e-4: zeros.append((P,Delta,gp))
results.sort()
print("  10 SMALLEST rho* found:")
for r,gp,k,P,Delta in results[:10]:
    print(f"    rho*={r:.5f} (k={k}, meas G_P={gp:.4f}) P={P} Delta={Delta}")
print(f"\n  MIN rho* = {best[0]:.5f}  at P={best[1]} Delta={best[2]}")
print(f"  #(rho* < 1e-4): {len(zeros)}")
for P,Delta,gp in zeros[:6]:
    print(f"     ZERO: P={P} Delta={Delta} meas(G_P)={gp:.4f}")
print("\nDONE.")
