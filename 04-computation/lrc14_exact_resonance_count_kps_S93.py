# EXACT-RESONANCE COUNT (kps-S93). The Sidon/AP dichotomy lives at n.e = 0 (EXACT additive
# relations), NOT the 1/(2N) shell. Count small-support n with:
#   (Z)  n.e = 0 exactly            <- Sidon-discriminating; these are the (1,-2,1)-type relations
#   (M)  n.e = 0 mod Vmax, n.e != 0 <- wraparound resonances (present for all)
# and the D_N-weighted resonant mass N * |sum_{n.e=0 mod Vmax} What(n)| vs target N(6/7)^k.
import math
from itertools import combinations, product
def babs(m):
    if m%7==0: return 0.0
    return abs(math.sin(math.pi*m/7))/(math.pi*abs(m))
def What_signed(nvec,K):
    active=[m for m in nvec if m!=0]; r=len(active)
    if any(m%7==0 for m in active): return 0.0
    sigma=sum(nvec)
    prodb=1.0
    for m in active: prodb*=babs(m)
    if sigma==0: sfac=6.0/7.0
    else:
        if sigma%7==0: return 0.0
        sfac=babs(sigma)
    return ((-1)**r)*((6.0/7.0)**((K-1)-r))*prodb*sfac
def counts(E,Vmax,K,supmax=3,Mmax=6):
    coords=list(range(1,K)); rng=[m for m in range(-Mmax,Mmax+1) if m!=0 and m%7!=0]
    zc=0; mc=0; z_mass=0.0; m_mass=0.0; res_signed=0.0
    for s in range(1,supmax+1):
        for combo in combinations(coords,s):
            for vals in product(rng,repeat=s):
                nvec=[0]*(K-1)
                for idx,ci in enumerate(combo): nvec[ci-1]=vals[idx]
                ws=What_signed(nvec,K)
                if ws==0.0: continue
                ne=sum(vals[idx]*E[combo[idx]] for idx in range(s))
                if ne==0:
                    zc+=1; z_mass+=abs(ws); res_signed+=ws
                elif ne%Vmax==0:
                    mc+=1; m_mass+=abs(ws); res_signed+=ws
    return zc,mc,z_mass,m_mass,res_signed
def W_of(E,Vmax,j,TH=1.0/7.0):
    pts=sorted(((e*j)%Vmax)/Vmax for e in E); n=len(pts); unc=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>TH: unc+=g-TH
    return unc
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
# Fixed families (deterministic) so the dichotomy is clean.
def make_AP(K,Vmax): d=Vmax//(K+1); return [d*i for i in range(K)]
def make_sidon(K,Vmax):
    # greedy Sidon set mod Vmax starting 0,1
    E=[0,1]; diffs={1,Vmax-1}
    x=2
    while len(E)<K and x<Vmax:
        ok=True; nd=set()
        for e in E:
            d1=(x-e)%Vmax; d2=(e-x)%Vmax
            if d1 in diffs or d2 in diffs or d1 in nd or d2 in nd: ok=False; break
            nd.add(d1); nd.add(d2)
        if ok: E.append(x); diffs|=nd
        x+=1
    return sorted(E)
print("EXACT-RESONANCE COUNT: Sidon vs AP. Z = #{n.e=0 exact}, M = #{n.e=0 mod V, !=0}.")
print("resMass = N*|sum_{n.e=0 mod V} What(n)| (the resonant Dirichlet mass, D_N=N terms).")
print(f"{'K':>3}{'family':>9}{'longAP':>7}{'r_N':>7}{'Z(exact)':>9}{'M(wrap)':>8}{'Zmass':>8}{'resMass/t':>10}")
for K in (11,12,13):
    N=math.ceil(7*(K-1)/6); lead=(6.0/7.0)**K; targ=N*lead; Vmax=1009  # prime ruler
    for name,E in [("AP",make_AP(K,Vmax)),("Sidon",make_sidon(K,Vmax))]:
        if len(set(E))!=K: 
            print(f"{K:>3}{name:>9}  (bad size {len(set(E))})"); continue
        SN=sum(W_of(E,Vmax,j) for j in range(1,N+1)); rN=abs(SN-targ)/targ
        zc,mc,zm,mm,rs=counts(E,Vmax,K)
        resmass=N*abs(rs)
        print(f"{K:>3}{name:>9}{longest_ap(E):>7}{rN:>7.3f}{zc:>9}{mc:>8}{zm:>8.3f}{resmass/targ:>10.3f}")
print()
print("EXPECT: AP has Z >> Sidon-Z (the (1,-2,1) 3-AP relations force n.e=0). If Sidon-Z is small")
print("and Sidon resMass/t < 1 => the COUNT route closes the dissociated branch a-priori.")
