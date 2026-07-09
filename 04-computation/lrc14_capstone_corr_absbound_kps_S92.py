# Capstone inequality: is the A-PRIORI ABSOLUTE bound Sum|W-hat(n)|*min(N,1/2||n.e/Vmax||) < N(6/7)^k ?
# Compare to the SIGNED Corr_N = S_N - N(6/7)^k (r_N). If absolute < 1, clean a-priori proof.
import math, random
from itertools import combinations, product
TH=1.0/7.0
def W_of(E,Vmax,j):
    pts=sorted(((e*j)%Vmax)/Vmax for e in E); n=len(pts); unc=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>TH: unc+=g-TH
    return unc
def frac_dist(x):  # ||x|| distance to nearest integer
    x=x-math.floor(x); return min(x,1-x)
def babs(m):  # |b0(m)| = |sin(pi m/7)|/(pi|m|); 0 if 7|m
    if m%7==0: return 0.0
    return abs(math.sin(math.pi*m/7))/(math.pi*abs(m))
def cabs(s):  # |c(sigma)|; sigma=0 -> use |1-c(0)|=6/7 handled separately
    if s==0: return None
    if s%7==0: return 0.0
    return abs(math.sin(math.pi*s/7))/(math.pi*abs(s))
def What_abs(nvec):  # |W-hat(n)| via LEM-011, n over coords e_1..e_{k-1}
    active=[m for m in nvec if m!=0]; r=len(active)
    if any(m%7==0 for m in active): return 0.0
    sigma=sum(nvec)
    prodb=1.0
    for m in active: prodb*=babs(m)
    if sigma==0: sfac=6.0/7.0
    else:
        if sigma%7==0: return 0.0
        sfac=cabs(sigma)
    kfac=(6.0/7.0)**( (K-1) - r )   # (6/7)^{k-1-r}
    return kfac*prodb*sfac
def abs_bound(E,Vmax,N,supmax=3,Mmax=6):
    # E = [e_0=0, e_1,...,e_{k-1}]; n over e_1..e_{k-1} (K-1 coords)
    coords=list(range(1,K))  # indices into E for e_1..e_{k-1}
    tot=0.0
    rng=[m for m in range(-Mmax,Mmax+1) if m!=0 and m%7!=0]
    for s in range(1,supmax+1):
        for combo in combinations(coords, s):
            for vals in product(rng, repeat=s):
                nvec=[0]*(K-1)
                for idx,ci in enumerate(combo): nvec[ci-1]=vals[idx]
                wa=What_abs(nvec)
                if wa==0: continue
                # n.e = sum vals[idx]*E[combo[idx]]
                ne=sum(vals[idx]*E[combo[idx]] for idx in range(s))
                theta=frac_dist(ne/Vmax)
                gN = N if theta<1e-12 else min(float(N), 1.0/(2*theta))
                tot+=wa*gN
    return tot
def isprim(E):
    from math import gcd; from functools import reduce
    E=sorted(E); return reduce(gcd,[e-E[0] for e in E])==1
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for jj in range(i+1,len(E)):
            d=E[jj]-E[i]; L=2; nx=E[jj]+d
            while nx in S: L+=1; nx+=d
            pv=E[i]-d
            while pv in S: L+=1; pv-=d
            if L>best: best=L
    return best
random.seed(165); K=13; N=math.ceil(7*(K-1)/6); lead=(6.0/7.0)**K
print(f"K={K} N={N} (6/7)^K={lead:.5f}; TARGET: |Corr_N| < N*(6/7)^K = {N*lead:.4f}")
print(f"{'config':>16}{'longAP':>7}{'r_N(signed)':>12}{'absBound/N(6/7)^k':>19}{'<1?':>5}")
Vmax=1001; lo=int(6*Vmax/7)
for label,gen in [("dissoc",0),("near-AP",1),("random",2)]:
    for _ in range(4):
        if gen==1:
            d=random.randint(lo//(K-1)+1,Vmax//(K-1)); L=random.randint(K-3,K-1)
            ap=[d*i for i in range(L)]
            if max(ap)>=Vmax: continue
            E=sorted(set(([0] if 0 not in ap else [])+ap+random.sample([x for x in range(1,Vmax) if x not in ap],K-L)))
        else:
            sp=random.randint(lo,Vmax-1); E=sorted(set([0]+random.sample(range(1,sp),K-2)+[sp]))
        if len(E)!=K: continue
        SN=sum(W_of(E,Vmax,j) for j in range(1,N+1)); corr=SN-N*lead; rN=abs(corr)/(N*lead)
        ab=abs_bound(E,Vmax,N); abr=ab/(N*lead)
        print(f"{label:>16}{longest_ap(E):>7}{rN:>12.3f}{abr:>19.3f}{'YES' if abr<1 else 'NO':>5}")
        break
print("\n=> if absBound/N(6/7)^k < 1 uniformly => CLEAN a-priori proof of |Corr_N|<N(6/7)^k => capstone.")
print("   if >1, the signed r_N<1 needs cancellation (abs too lossy) -- characterize the gap.")
