"""Witnesses for the lower bound q_R(T) >= c_R T^(-1/2): a residue class y = rho mod 2^t0 on which the partners
y+1..y+R have all merged with each other (equal Terras time, uniformly on the class) by time t0 while y has merged with
none of them.  Merges are checked on three random lifts (uniformity) of the class; then the y-vs-cluster chain state
at t0 is printed.  Uses random 2-adic y (600-bit lifts, t0 <= 400)."""
import random, sys
def T(x): return x//2 if x%2==0 else (3*x+1)//2
def first_merge_times(y, R, L):
    vals=[y+r for r in range(R+1)]
    # cluster union-find over partners 1..R; record y's first merge
    parent=list(range(R+1))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
        return a
    ncl=R; ty=None
    for t in range(1,L+1):
        vals=[T(v) for v in vals]
        seen={}
        for i,v in enumerate(vals):
            if v in seen:
                a,b=find(i),find(seen[v])
                if a!=b:
                    if 0 in (i,seen[v]) or find(0) in (a,b):
                        if ty is None: ty=t
                    parent[a]=b
            else: seen[v]=i
        # number of partner clusters not containing y
        roots=set(find(i) for i in range(1,R+1))
        if ty is None and len(roots)==1: return t, ty, vals
        if ty is not None: return None, ty, vals
    return None, ty, vals
random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
Rmax=int(sys.argv[1])
for R in range(int(sys.argv[3]) if len(sys.argv)>3 else 2, Rmax+1):
    found=None
    for trial in range(20000):
        rho=random.getrandbits(400)
        t0,ty,_=first_merge_times(rho, R, 400)
        if t0 is None: continue
        # uniformity: same t0 behaviour on other lifts rho + 2^t0 * s
        ok=True
        for s in range(3):
            y=(rho % (1<<t0)) + (random.getrandbits(200)<<t0) + (1 << (t0+300))
            t1,ty1,_=first_merge_times(y, R, t0)
            if t1!=t0: ok=False;break
        if ok: found=(trial,t0,rho % (1<<t0)); break
    if found:
        trial,t0,cls=found
        print(f"R={R}: witness after {trial+1} draws: class y = {cls} mod 2^{t0} (partners coalesced at Terras time {t0}, y unmerged)", flush=True)
    else:
        print(f"R={R}: no witness in 20000 draws", flush=True)
