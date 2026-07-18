from fractions import Fraction as F
from math import gcd
def norm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return min(r,1-r)
def minnorm(fam,t): return min(norm(F(v)*t) for v in fam)
def Pmax(P):
    best=F(0); t0=None
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            v=minnorm(P,F(a,q))
            if v>best: best=v; t0=F(a,q)
    return best,t0
def bad_union_frac(P,K,N=6000):
    mu,t0=Pmax(P); smax=(mu-F(1,14))/max(P)
    if smax<=0: return None
    bad=0
    for i in range(N):
        s=smax*F(2*i-N,N); t=t0+s
        if minnorm(K,t)<F(1,14): bad+=1
    return mu,t0,smax,bad/N

# UNIFORM-IN-j test: NARROW near-equal fast clusters of growing j over shrinking core
# cluster centered ~ C, consecutive (spread = j-1, narrow), fast
print("Narrow near-equal fast clusters, growing j -- does bad-union stay << 1 (collapse), NOT ~ j/7?")
for j in range(2,11):
    c=13-j; P=list(range(1,c+1))
    C=400
    K=list(range(C, C+j))          # consecutive: spread j-1
    r=bad_union_frac(P,K)
    if r is None: print(f"  j={j}: no room"); continue
    mu,t0,smax,bf=r
    Delta=K[-1]-K[0]
    Fmeas_pred=F(1,7)+Delta*(t0+smax)   # predicted meas(F): 1/7 + spread
    fast_ok = K[0]*2*smax>=1
    print(f"  j={j:2d} core|{c}| K=[{K[0]}..{K[-1]}] Δ={Delta}: bad-union={bf:.3f}  "
          f"(union-bd j/7={j/7:.2f}, moiré-pred 1/7+Δ(t0+smax)={float(Fmeas_pred):.3f})  "
          f"fast:{fast_ok} good-kick:{bf<1.0}")
