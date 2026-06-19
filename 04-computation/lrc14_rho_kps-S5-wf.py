"""
Verify the claim's ρ*(P,E) = meas(G_P ∩ Good(E)) computations EXACTLY.
G_P = {x : ||p x|| >= 1/14 for all p in P}, ||.|| = distance to nearest integer.
Good(E) = {x : circular max-gap of {frac(e x)} > 2/7}.

Claimed:
  rho*=0 for P={1,2,3}, E=(0,1,2,3,4,5,6,7,8,10): mu(E)=121/490, meas(G_P)=29/42, rho*=0
  second: P={1,2,3,4}, E=(0,2,3,4,5,6,7,8,10): mu(E)=27/98, meas(G_P)=13/21, rho*=0
We recompute meas(G_P), mu(E), and rho* EXACTLY via shared breakpoints.
kps-S5-wf
"""
from fractions import Fraction as F

def frac_norm_breakpoints(P):
    # ||p x|| = 1/14 boundaries: p x = integer +- 1/14 => x = (n +- 1/14)/p
    bps=set([F(0),F(1)])
    for p in P:
        for n in range(0,p+1):
            for s in (F(1,14),-F(1,14)):
                x=(F(n)+s)/p
                if 0<=x<=1: bps.add(x)
        # also where p x integer (order anchors) not needed but harmless
        for n in range(0,p+1): bps.add(F(n,p))
    return sorted(bps)

def in_GP(x, P):
    for p in P:
        v=(p*x)%1
        d=min(v,1-v)
        if d < F(1,14):
            return False
    return True

def meas_GP(P):
    bps=frac_norm_breakpoints(P)
    tot=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        if in_GP(mid,P): tot+=(b-a)
    return tot

def good_breakpoints(E):
    E=sorted(set(E)); k=len(E)
    TH=F(2,7); bp=set([F(0),F(1)]); diffs=set()
    for i in range(k):
        for j in range(i+1,k): diffs.add(E[j]-E[i])
    for d in diffs:
        for m in range(0,d+1): bp.add(F(m,d))
    obp=sorted(b for b in bp if F(0)<=b<=F(1)); refined=set(obp)
    for a,b in zip(obp,obp[1:]):
        mid=(a+b)/2; floors={e:(e*mid).__floor__() for e in E}
        order=sorted(E,key=lambda e:e*mid-floors[e])
        for t in range(k):
            if t==k-1:
                ef,el=order[0],order[-1]; slope=ef-el; const=F(1)-floors[ef]+floors[el]
            else:
                eh,elo=order[t+1],order[t]; slope=eh-elo; const=-floors[eh]+floors[elo]
            if slope!=0:
                xb=(TH-const)/slope
                if a<xb<b: refined.add(xb)
    return sorted(refined)

def in_good(x,E):
    pts=sorted(set((e*x)%1 for e in E))
    if len(pts)==1: return True
    g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1])
    return max(g)>F(2,7)

def mu_exact(E):
    bps=good_breakpoints(E); tot=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        if in_good(mid,E): tot+=(b-a)
    return tot

def rho_star(P,E):
    # merge breakpoints of G_P and Good(E)
    bps=set(frac_norm_breakpoints(P)) | set(good_breakpoints(E))
    bps=sorted(b for b in bps if 0<=b<=1)
    tot=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        if in_GP(mid,P) and in_good(mid,E):
            tot+=(b-a)
    return tot

if __name__=="__main__":
    cases=[
        ("P={1,2,3}, E=(0..8,10)", [1,2,3], [0,1,2,3,4,5,6,7,8,10], F(121,490), F(29,42)),
        ("P={1,2,3,4}, E=(0,2..8,10)", [1,2,3,4], [0,2,3,4,5,6,7,8,10], F(27,98), F(13,21)),
    ]
    for name,P,E,mu_claim,gp_claim in cases:
        mu=mu_exact(E); gp=meas_GP(P); rs=rho_star(P,E)
        print(f"{name}:", flush=True)
        print(f"   mu(E)={mu}={float(mu):.6f}  claim {mu_claim}  {'OK' if mu==mu_claim else 'MISMATCH'}", flush=True)
        print(f"   meas(G_P)={gp}={float(gp):.6f}  claim {gp_claim}  {'OK' if gp==gp_claim else 'MISMATCH'}", flush=True)
        print(f"   rho*={rs}={float(rs):.6f}  claim 0  {'OK (exactly disjoint)' if rs==0 else 'MISMATCH (NONZERO!)'}", flush=True)
