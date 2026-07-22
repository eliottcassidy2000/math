import numpy as np
from math import log
# r2(k) = # of (a,b) in Z^2 with a^2+b^2=k  (ordered, signed) = 4*(d1-d3)
def r2(k):
    if k==0: return 1
    c=0; a=0
    while a*a<=k:
        b2=k-a*a; b=int(round(b2**0.5))
        if b*b==b2:
            # count (±a,±b) and (±b,±a) properly -- just brute over b too
            pass
        a+=1
    # brute (small k)
    cnt=0; s=int(k**0.5)+1
    for x in range(-s,s+1):
        for y in range(-s,s+1):
            if x*x+y*y==k: cnt+=1
    return cnt

print("="*70)
print("ANCHOR 1: r2 maximizer over k<=N  =  product of primes = 1 (mod 4)")
def factor(k):
    f={}; d=2
    while d*d<=k:
        while k%d==0: f[d]=f.get(d,0)+1; k//=d
        d+=1
    if k>1: f[k]=f.get(k,0)+1
    return f
for N in (50,200,1000,5000,20000):
    best_k,best=max(((k,r2(k)) for k in range(1,N+1)), key=lambda t:t[1])
    fac=factor(best_k)
    primes_mod4 = {p:(p%4) for p in fac}
    print(f"  N={N:6d}: argmax_k r2 = {best_k:6d} = {fac}  r2={best:3d}  primes mod4={primes_mod4}")

print("\n"+"="*70)
print("ANCHOR 2: grid lower bound  u(m x m grid) / n  grows (n^{c/loglog n})")
print("  best squared-distance k*, unit-distance count U, U/n, and k* factorization")
for m in (5,10,20,40,70,100):
    n=m*m
    # count unordered pairs at each squared distance via displacement multiplicities
    from collections import defaultdict
    cnt=defaultdict(float)
    for dx in range(0,m):
        for dy in range(0,m):
            if dx==0 and dy==0: continue
            k=dx*dx+dy*dy
            mult = (m-dx)*(m-dy)
            # displacement (dx,dy) with signs: (dx,dy),(dx,-dy),(-dx,dy),(-dx,-dy) but unordered pairs:
            # number of unordered pairs with |Δx|=dx,|Δy|=dy:
            if dx>0 and dy>0: pairs = 2*mult      # (dx,dy) and (dx,-dy) give distinct pair-sets
            else: pairs = mult                     # axis-aligned: one orientation
            cnt[k]+=pairs
    kbest=max(cnt, key=lambda k:cnt[k]); U=cnt[kbest]
    print(f"  m={m:3d} n={n:5d}: k*={kbest:5d}={factor(kbest)}  U={int(U):7d}  U/n={U/n:6.2f}  r2(k*)={r2(kbest)}")

print("\n"+"="*70)
print("ANCHOR 3: incidence structure -- unit-distance bipartite graph is K_{2,3}-free")
print("  (two points share <=2 common unit-neighbors: two unit circles meet in <=2 pts)")
np.random.seed(0)
def unit_graph_random(n, scale):
    P=np.random.rand(n,2)*scale
    # count unit distances (tolerance) -- use exact rational-ish: place on fine grid
    U=0; commonmax=0
    D=np.sqrt(((P[:,None,:]-P[None,:,:])**2).sum(-1))
    adj=(np.abs(D-1.0)<1e-9)
    return adj
# geometric K_{2,3}-free check on a structured set (grid scaled so many unit dists)
m=15; pts=np.array([(i,j) for i in range(m) for j in range(m)],float)
k=5  # squared unit distance 5 (dist sqrt5), r2(5)=8
ud=1.0*(np.abs(((pts[:,None]-pts[None])**2).sum(-1)-k)<1e-9)
# common neighbors of each pair
worst=0
for a in range(len(pts)):
    na=np.where(ud[a]>0)[0]
    for b in range(a+1,len(pts)):
        common=np.intersect1d(na, np.where(ud[b]>0)[0])
        worst=max(worst,len(common))
print(f"  grid {m}x{m}, dist=sqrt{k}: max common unit-neighbors over all pairs = {worst}  (<=2 confirms K_2,3-free)")
print("  => Kovari-Sos-Turan gives O(n^3/2); Szemeredi-Trotter refines to O(n^4/3).")
