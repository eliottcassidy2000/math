"""Signed signature of the WORRY-SET (tight, M=1/n). Do tight configs carry shell-partners
(v_i+v_j≡0 mod 2n-1, = zero signed pair-clocks)? Is the signed structure a finer invariant that
separates tight configs sharing the same M? Build n=4..7 exhaustively (window), report each tight
config's shell-partners and sign-orbit. opus-2026-06-06-S699d."""
from itertools import combinations, product
from fractions import Fraction as F
from math import gcd
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for m in range(d):
            t=F(m,d); v=min(min((x*t)%1,1-(x*t)%1) for x in V)
            if v>best: best=v
    return best
def shell_partners(V,C): return [(V[i],V[j]) for i in range(len(V)) for j in range(i+1,len(V)) if (V[i]+V[j])%C==0]
def sign_orbit(V,C):
    orbit=set()
    for eps in product((1,-1),repeat=len(V)):
        cl=tuple(sorted(min((eps[i]*V[i]-eps[j]*V[j])%C,(eps[j]*V[j]-eps[i]*V[i])%C) for i in range(len(V)) for j in range(i+1,len(V))))
        orbit.add(cl)
    return len(orbit)
def main():
    print("WORRY-SET signed signature: tight configs (M=1/n), their shell-partners & sign-orbit")
    for n in range(4,8):
        C=2*n-1; thr=F(1,n); B=2*n; tight=[]
        for V in combinations(range(1,B+1),n-1):
            g=0
            for v in V: g=gcd(g,v)
            if g!=1: continue
            if M_exact(V)==thr: tight.append(V)
        print(f" n={n} (C={C}): {len(tight)} tight configs in [1,{B}]")
        sp_yes=sp_no=0
        for V in tight:
            sp=shell_partners(V,C); orb=sign_orbit(V,C)
            if sp: sp_yes+=1
            else: sp_no+=1
            print(f"    {V}: shell-partners={sp if sp else 'NONE'}; sign-orbit={orb}")
        print(f"    => tight with shell-partner: {sp_yes}, without: {sp_no}")
    print("\nReading: if ALL tight configs have NO shell-partner, that's a worry-set theorem candidate")
    print("(tight ⟹ no pair sums to 2n-1). If some do, the signed structure splits the worry-set.")
if __name__=='__main__': main()
