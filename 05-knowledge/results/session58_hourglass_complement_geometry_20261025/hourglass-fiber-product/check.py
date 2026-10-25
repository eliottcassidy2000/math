from itertools import product
from fractions import Fraction
import json
from pathlib import Path

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

try:
    require(False, "intentional guard probe")
except RuntimeError as exc:
    require(str(exc) == "intentional guard probe", "guard probe wrong exception")
else:
    raise RuntimeError("guard probe failed to raise")

def test(p):
    roots={a:[t for t in range(p) if t*t%p==a] for a in range(p)}
    r=lambda a:len(roots[a%p])
    sums=[0,0,0,0]; entries=[0,0,0,0]; even_energy=Fraction(0); odd_energy=Fraction(0)
    for E,x,z in product(range(p),repeat=3):
        common=[E,E+x,E-x,E+z,E-z]
        plus=[E+x+z,E-x-z]; minus=[E+x-z,E-x+z]
        w=1
        for a in common:w*=r(a)
        hp=r(plus[0])*r(plus[1]); hm=r(minus[0])*r(minus[1])
        even_energy+=w*Fraction(hp+hm,2)**2
        odd_energy+=w*Fraction(hp-hm,2)**2
        for i,v in enumerate([w,w*hp,w*hm,w*hp*hm]):sums[i]+=v
        for i,forms in enumerate([common,common+plus,common+minus,common+plus+minus]):
            entries[i]+=all(r(a)>0 for a in forms)
    direct=[0,0,0,0]
    if p<=7:
        for c,a,b,d,e in product(range(p),repeat=5):
            if (a*a+b*b-2*c*c)%p or (d*d+e*e-2*c*c)%p:continue
            hp=r(a*a+d*d-c*c)*r(3*c*c-a*a-d*d)
            hm=r(a*a-d*d+c*c)*r(c*c-a*a+d*d)
            for i,v in enumerate([1,hp,hm,hp*hm]): direct[i]+=v
        require(sums==direct, "independent shared-root enumeration")
    proj=[(v-1)//(p-1) for v in sums]
    require(all((v-1)%(p-1)==0 for v in sums), "projective divisibility")
    even_energy=(even_energy-1)/(p-1)
    odd_energy=odd_energy/(p-1)
    require(even_energy-odd_energy==proj[3], "spectral energy identity")
    cov=Fraction(proj[3],proj[0])-Fraction(proj[1]*proj[2],proj[0]**2)
    # The nine label coefficient vectors remain distinct, nonzero, and span rank 3.
    forms=[(1,i%p,j%p) for i,j in [(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]]
    require(len(set(forms))==9, "distinct parameter forms")
    return dict(p=p,cone_roots=sums,projective_roots=proj,affine_square_entry_parameters=entries,projective_covariance=str(cov),even_energy=str(even_energy),odd_energy=str(odd_energy),direct_shared_root_check=p<=7)

# Symbolic linear reconstruction and six independent equations.
import sympy as s
E,x,z=s.symbols('E x z')
labels=s.Matrix([E,E+x,E-x,E+z,E-z,E+x+z,E-x-z,E+x-z,E-x+z])
C=s.Matrix([[-2,1,1,0,0,0,0,0,0],[-2,0,0,1,1,0,0,0,0],[1,-1,0,-1,0,1,0,0,0],[-3,1,0,1,0,0,1,0,0],[-1,-1,0,1,0,0,0,1,0],[-1,1,0,-1,0,0,0,0,1]])
require(C*labels==s.zeros(6,1), "symbolic reconstruction")
require(C[:,[2,4,5,6,7,8]].det()==1, "unit pivot minor")
results={'symbolic_rank_all_characteristics':6,'monic_pivot_determinant':1,'finite_fields':[test(p) for p in [3,5,7,11,13,17,19,23,29,31]]}
open(Path(__file__).with_name('results.json'),'w').write(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
