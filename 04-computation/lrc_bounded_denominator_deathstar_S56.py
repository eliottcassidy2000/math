from math import gcd
n=14
def loneliness_times(fam, Dmax):
    # all reduced t=k/D (D<=Dmax) with min_v ||v t|| == 1/n exactly (=tight loneliness time)
    out=[]
    for D in range(2,Dmax+1):
        for k in range(1,D):
            if gcd(k,D)!=1: continue
            m=min(min((v*k)%D, D-(v*k)%D) for v in fam)
            if n*m==D:   # ||.|| == 1/n exactly
                out.append((k,D))
    return out
def active(fam,k,D):
    return [v for v in fam if min((v*k)%D,D-(v*k)%D)*n==D]

for name,fam in [("AP {1..13}",list(range(1,14))),("GW {1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    Vmax=max(fam); L=loneliness_times(fam,120)
    denoms=sorted(set(D for _,D in L))
    print(f"{name}: Vmax={Vmax}, 2*Vmax={2*Vmax}")
    print(f"   all loneliness times (denom<=120): distinct denominators = {denoms}")
    print(f"   max denom = {max(denoms)}  <=  2*Vmax = {2*Vmax} ?  {max(denoms)<=2*Vmax}")
    print(f"   any denom > n=14 ?  {any(D>14 for D in denoms)}   => R holds: {not any(D>14 for D in denoms)}")
    # show the v+ + v- = denom structure at t=1/14
    act=active(fam,1,14); print(f"   active runners at t=1/14: {act}, sum={sum(act)} (=denom 14)")
    print()
# The bound in action: R for AP needs NO check (2*Vmax=26<28); for GW only D in {28,42} to rule out.
print("R is a FINITE check: denom<=2Vmax and denom divisible by 14 => AP: none in (14,26]; GW: only {28,42}, both absent.")
