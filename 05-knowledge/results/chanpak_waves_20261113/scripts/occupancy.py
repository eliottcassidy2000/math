"""Wave 25: bulk occupancy of the (h, f) profile set of order n, counting only ADMISSIBLE lattice points:
h odd, f = 0 mod 4 (mod 8 for even n), f = (-1)^n h mod 3. Input: lines 'h f count' (profhf) or 'h f F3' (prof)."""
import sys, collections, numpy as np
def load(path,n,counted):
    cnt=collections.Counter()
    for line in open(path):
        a=line.split(); h,f=int(a[0]),int(a[1]); cnt[(h,f)]+= int(a[2]) if counted else 1
    return cnt
for path,n,counted in [(a.split(':')[0],int(a.split(':')[1]),a.split(':')[2]=='c') for a in sys.argv[1:]]:
    cnt=load(path,n,counted); P=set(cnt)
    H=np.array([p[0] for p in P],float); F=np.array([p[1] for p in P],float); Wt=np.array([cnt[p] for p in P],float)
    mh=np.average(H,weights=Wt); mf=np.average(F,weights=Wt)
    sh=np.sqrt(np.average((H-mh)**2,weights=Wt)); sf=np.sqrt(np.average((F-mf)**2,weights=Wt))
    step=8 if n%2==0 else 4; eps=1 if n%2==0 else -1
    bad=sum(1 for h,f in P if h%2==0 or f%step or (f-eps*h)%3)
    out=[]
    for r in (0.1,0.25,0.5):
        tot=occ=0
        for h in range(int(mh-r*sh)|1,int(mh+r*sh)+1,2):
            for f in range((int(mf-r*sf)//step+1)*step,int(mf+r*sf)+1,step):
                if (f-eps*h)%3: continue
                tot+=1; occ+=((h,f) in P)
        out.append(f"+-{r}sd: {occ}/{tot} = {occ/max(tot,1):.3f}")
    print(f"n={n}: {len(P)} profiles ({bad} violate the congruences); center ({mh:.0f},{mf:.0f}), sd ({sh:.0f},{sf:.0f}); occupancy "+"; ".join(out))
