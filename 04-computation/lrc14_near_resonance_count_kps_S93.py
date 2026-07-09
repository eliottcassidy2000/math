# NEAR-RESONANCE COUNT (kps-2026-07-09-S93). Work the last-mile capstone core with the
# MERTENS lens: split Corr_N = sum_n What(n) G_N(theta_n) into NEAR (||theta||<1/(2N)) and
# FAR (>=1/(2N)) resonances. Question: does the 20x absolute-bound bloat live in the NEAR
# part (a COUNT problem -- Mertens-SAFE) or the FAR part (a CANCELLATION problem -- the
# Mertens TRAP, where a signed sum's sqrt-cancellation heuristic can FAIL, cf. Mertens conj)?
import math, random
from itertools import combinations, product
TH=1.0/7.0
def frac(x): x=x-math.floor(x); return x
def frac_dist(x): x=x-math.floor(x); return min(x,1-x)
def babs(m):
    if m%7==0: return 0.0
    return abs(math.sin(math.pi*m/7))/(math.pi*abs(m))
def bsigned(m):  # b0(m) real part sign via sin; use the actual complex b0? We track sign via (-1)^r and c(sigma).
    return babs(m)
def What_signed(nvec,K):
    # LEM-011: What(n) = (-1)^r (6/7)^{k-1-r} [prod b0(n_i)] (1[sigma=0]-c(sigma))
    # For the SIGNED real magnitude we use: sign (-1)^r, prod|b0|, and (1-c) factor.
    # c(sigma)= (6/7)-analog; here we use the ABS structure but keep the (-1)^r sign and the
    # (1[sigma=0]-c) sign. We approximate the phase-real-part by its magnitude*sign((-1)^r).
    active=[m for m in nvec if m!=0]; r=len(active)
    if any(m%7==0 for m in active): return 0.0
    sigma=sum(nvec)
    prodb=1.0
    for m in active: prodb*=babs(m)
    if sigma==0: sfac=6.0/7.0
    else:
        if sigma%7==0: return 0.0
        sfac=babs(sigma)  # |c(sigma)|
    kfac=(6.0/7.0)**((K-1)-r)
    return ((-1)**r)*kfac*prodb*sfac
def W_of(E,Vmax,j):
    pts=sorted(frac((e*j)/Vmax) for e in E); n=len(pts); unc=0.0
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
def resonance_split(E,Vmax,N,K,supmax=3,Mmax=6):
    coords=list(range(1,K)); rng=[m for m in range(-Mmax,Mmax+1) if m!=0 and m%7!=0]
    near_abs=far_abs=0.0; near_cnt=0; far_cnt=0
    near_signed=far_signed=0.0
    cut=1.0/(2*N)
    for s in range(1,supmax+1):
        for combo in combinations(coords,s):
            for vals in product(rng,repeat=s):
                nvec=[0]*(K-1)
                for idx,ci in enumerate(combo): nvec[ci-1]=vals[idx]
                ws=What_signed(nvec,K)
                if ws==0.0: continue
                ne=sum(vals[idx]*E[combo[idx]] for idx in range(s))
                theta=frac_dist(ne/Vmax)
                # G_N(theta) real magnitude bound: min(N, 1/(2 sin(pi theta))) ~ min(N,1/(2 theta_dist))
                gmag = N if theta<1e-12 else min(float(N), 1.0/(2*math.sin(math.pi*theta)))
                # signed real contribution ~ ws * (Dirichlet kernel real part); approx sign by cos(pi(N+1)theta)?
                # We track |contribution| for abs, and ws*gmag*cos-phase is too detailed; use ABS split + COUNT.
                contrib_abs=abs(ws)*gmag
                if theta<cut:
                    near_abs+=contrib_abs; near_cnt+=1
                else:
                    far_abs+=contrib_abs; far_cnt+=1
    return near_abs,far_abs,near_cnt,far_cnt
random.seed(93)
print("NEAR-RESONANCE COUNT + near/far split of the a-priori absolute bound.")
print("cut = 1/(2N); NEAR = ||n.e/Vmax|| < cut (the resonances), FAR = the rest.")
print(f"{'K':>3}{'cfg':>10}{'longAP':>7}{'r_N':>7}{'targ':>7}{'nearAbs/t':>10}{'farAbs/t':>10}{'#near':>7}{'#far':>7}")
for K in (11,12,13):
    N=math.ceil(7*(K-1)/6); lead=(6.0/7.0)**K; targ=N*lead; Vmax=1001; lo=int(6*Vmax/7)
    for label,gen in [("dissoc",0),("AP",1),("random",2)]:
        got=False
        for _try in range(40):
            if gen==1:
                d=random.randint(max(2,lo//(K-1)+1),max(3,Vmax//(K-1))); L=K
                ap=[d*i for i in range(L)]
                if max(ap)>=Vmax: continue
                E=sorted(set(ap))
            else:
                sp=random.randint(lo,Vmax-1); E=sorted(set([0]+random.sample(range(1,sp),K-2)+[sp]))
            if len(E)!=K: continue
            got=True; break
        if not got: continue
        SN=sum(W_of(E,Vmax,j) for j in range(1,N+1)); rN=abs(SN-targ)/targ
        na,fa,nc,fc=resonance_split(E,Vmax,N,K)
        print(f"{K:>3}{label:>10}{longest_ap(E):>7}{rN:>7.3f}{targ:>7.3f}{na/targ:>10.3f}{fa/targ:>10.3f}{nc:>7}{fc:>7}")
print()
print("READ: if the 20x lives in FAR (far_abs/t >> near_abs/t) => the bloat is CANCELLATION among")
print("many far resonances -- the MERTENS TRAP (cannot assume sqrt-cancellation a-priori).")
print("If NEAR dominates and #near is SMALL for dissoc => the COUNT route works (Mertens-SAFE).")
