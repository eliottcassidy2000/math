"""Wave 25: occupancy of admissible lattice points in a band around the ridge of the (h, f) cloud of order n.
Ridge: f ~ alpha + beta h (class-weighted least squares); residual sd s. Band: |h - mh| <= r_h sd_h and
|f - ridge(h)| <= r s. Admissible: h odd, f = 0 mod 4 (8 if n even), f = (-1)^n h mod 3."""
import sys, collections, numpy as np
for arg in sys.argv[1:]:
    path,n,mode=arg.split(':'); n=int(n)
    cnt=collections.Counter()
    for line in open(path):
        a=line.split(); cnt[(int(a[0]),int(a[1]))]+= int(a[2]) if mode=='c' else 1
    P=set(cnt); H=np.array([p[0] for p in P],float); F=np.array([p[1] for p in P],float); W=np.array([cnt[p] for p in P],float)
    mh=np.average(H,weights=W); sh=np.sqrt(np.average((H-mh)**2,weights=W))
    beta=np.average((H-mh)*(F-np.average(F,weights=W)),weights=W)/sh**2; alpha=np.average(F,weights=W)-beta*mh
    s=np.sqrt(np.average((F-alpha-beta*H)**2,weights=W))
    step=8 if n%2==0 else 4; eps=1 if n%2==0 else -1
    res=[]
    for rh,rs in ((0.25,0.25),(0.5,0.5),(1,1)):
        tot=occ=0
        for h in range(int(mh-rh*sh)|1,int(mh+rh*sh)+1,2):
            c=alpha+beta*h
            for f in range((int(c-rs*s)//step+1)*step,int(c+rs*s)+1,step):
                if (f-eps*h)%3: continue
                tot+=1; occ+=((h,f) in P)
        res.append(f"band({rh},{rs}): {occ}/{tot} = {occ/max(tot,1):.3f}")
    print(f"n={n}: ridge f = {alpha:.0f} + {beta:.2f} h, residual sd {s:.0f} (sd_f|h / sd_f = {s/np.sqrt(np.average((F-np.average(F,weights=W))**2,weights=W)):.3f}); "+"; ".join(res))
