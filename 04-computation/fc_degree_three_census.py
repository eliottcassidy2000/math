import sys; sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation')
from fc_factorial_conjecture_census import L_of_power
from itertools import product
print("DEGREE-3 FC CENSUS (the regime where the determinacy argument fails)")
print()
# homogeneous cubics in 2 variables
mons=[(3,0),(2,1),(1,2),(0,3)]
for R,M in ((6,5),(9,4)):
    surv=[]
    for v in product(range(-R,R+1),repeat=4):
        if not any(v): continue
        f={mons[i]:v[i] for i in range(4) if v[i]}
        if all(L_of_power(f,m,2)==0 for m in range(1,M+1)): surv.append(f)
    print(f"  homogeneous cubic, n=2, coeffs[-{R},{R}], m=1..{M}: "
          f"{len(surv)} survivors / {(2*R+1)**4-1}")
    for f in surv[:5]: print(f"      {f}")
# cubic binomials: only two monomials, wide range
print()
best=[]
for i in range(4):
    for j in range(i+1,4):
        for A in range(-40,41):
            for B in range(-40,41):
                if A==0 or B==0: continue
                f={mons[i]:A,mons[j]:B}
                if L_of_power(f,1,2)==0 and L_of_power(f,2,2)==0:
                    best.append((mons[i],A,mons[j],B))
print(f"  two-term cubics with L(f)=L(f^2)=0, coeffs in [-40,40]: {len(best)}")
for b in best[:6]: print(f"      {b}")
