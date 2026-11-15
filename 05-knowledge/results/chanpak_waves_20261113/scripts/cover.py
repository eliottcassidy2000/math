"""Wave 25 (P66, P67, P69; with argv[2] = 10 also P70, P71): odd N <= X attained by the cyclic join T = A -> B -> c -> A,
H(T) = (f_A - h_A)(f_B - h_B) + 2 h_A h_B, with cores A, B from the profile census of orders <= 9 and padding
(h, f) -> (h, 2^k f). r(N) counts ordered pairs of distinct profile VALUES (A, B).
'bulk' marks unpadded profiles of order 6..9 within 2 residual sd of their order's ridge (band.py fit)."""
import numpy as np, collections, sys
X=int(sys.argv[1]) if len(sys.argv)>1 else 200000
prof={}   # (h,f) -> bulk flag (True if some realization is an unpadded bulk profile)
ridge={}
MAXN=int(sys.argv[2]) if len(sys.argv)>2 else 9
for n in range(1,MAXN+1):
    pts=collections.Counter()
    if n<=9:
        for line in open(f'prof_{n}.txt'):
            h,f,_=map(int,line.split()); pts[(h,f)]+=1
    else:
        for line in open(f'profhf_{n}.txt'):
            h,f,c=map(int,line.split()); pts[(h,f)]+=c
    if n>=6:
        H=np.array([p[0] for p in pts],float); F=np.array([p[1] for p in pts],float); W=np.array(list(pts.values()),float)
        mh=np.average(H,weights=W); sh=np.sqrt(np.average((H-mh)**2,weights=W))
        beta=np.average((H-mh)*(F-np.average(F,weights=W)),weights=W)/sh**2; alpha=np.average(F,weights=W)-beta*mh
        s=np.sqrt(np.average((F-alpha-beta*H)**2,weights=W)); ridge[n]=(alpha,beta,s)
    for (h,f) in pts:
        bulk = n>=6 and abs(f-(ridge[n][0]+ridge[n][1]*h))<=2*ridge[n][2]
        prof[(h,f)] = prof.get((h,f),False) or bulk
        k=1; g=2*f
        while g-h<=X:
            prof.setdefault((h,g),False); k+=1; g*=2
P=list(prof.items())
U=np.array([f-h for (h,f),b in P],dtype=np.int64); V=np.array([h for (h,f),b in P],dtype=np.int64); Bk=np.array([b for p,b in P])
o=np.argsort(U); U,V,Bk=U[o],V[o],Bk[o]
print("profile values used:",len(P)," of which bulk:",int(Bk.sum()))
r=np.zeros(X//2+1,dtype=np.int64); rb=np.zeros(X//2+1,dtype=np.int64)
for i in range(len(U)):
    x,y=int(U[i]),int(V[i])
    j=np.searchsorted(U,X//x,side='right')
    if j==0: break
    Ns=x*U[:j]+2*y*V[:j]; m=Ns<=X
    np.add.at(r,Ns[m]//2,1)
    mb=m&Bk[:j]
    np.add.at(rb,Ns[mb]//2,1)
odd=np.arange(1,X+1,2)
miss=odd[r[odd//2]==0]; missb=odd[rb[odd//2]==0]
print(f"X={X}: missed odd N: {len(miss)}; the first 40: {miss[:40].tolist()}; largest {miss.max() if len(miss) else None}")
print(f"  misses above 1000: {miss[miss>1000].tolist()[:30]} (count {int((miss>1000).sum())})")
print(f"  bulk-B only: misses above 1000: {int((missb>1000).sum())}; largest {missb.max() if len(missb) else None}")
from sympy import primerange
pr=np.array(list(primerange(100000,X+1)))
if len(pr): print(f"  primes in [1e5,{X}]: min r = {int(r[pr//2].min())} (at {int(pr[np.argmin(r[pr//2])])}), mean r = {r[pr//2].mean():.1f}; all odd N in that window: mean r = {r[odd[odd>=100000]//2].mean():.1f}")
np.save(f'r_cyclic_{MAXN}.npy',r)
for k in range(10,24):
    a,b=1<<k,min(1<<(k+1),X+1); w=odd[(odd>=a)&(odd<b)]; rr=r[w//2]
    if len(w)==0: continue
    print(f'  [2^{k},2^{k+1}): missed {int((rr==0).sum())} of {len(w)}; mean r {rr.mean():.2f}; min r {int(rr.min())}; missed by N mod 3: '+' '.join(str(int(((rr==0)&(w%3==c)).sum())) for c in range(3)))
