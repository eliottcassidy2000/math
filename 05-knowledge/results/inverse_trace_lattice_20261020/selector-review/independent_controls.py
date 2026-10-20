#!/usr/bin/env python3
"""Small independent exact controls; not a test of an actual noninvertible Keller map."""
from pathlib import Path
import json
import sympy as s
u,v,z,t,p,q=s.symbols('u v z t p q')
checks=[]
def check(condition,label):
    if not bool(condition):
        raise RuntimeError(label)
    checks.append(label)

# Specialized inverse polynomials, independently built with arbitrary low coefficients.
for d in (3,5,6,9):
    low=sum((j+2)*z**(d-j) for j in range(1,d-1))
    R=u*(z**d+low)+(v-1)*z+(2+u*v)
    for uv,expected in [((1,0),d),((0,0),1),((0,1),0)]:
        specialized=s.Poly(R.subs({u:uv[0],v:uv[1]}),z)
        check(specialized.degree()==expected,f'specialization d={d},p={uv}')
    linear=s.expand(R.subs(u,0))
    check(s.cancel(linear.subs(z,-2/(v-1)))==0,f'finite selector d={d}')
    zero_fiber=s.expand((u*(z**d+low)+v*z+u+v).subs({u:0,v:0}))
    check(zero_fiber==0,f'triple-zero fiber d={d}')

# Unit statement remains valid for a nonreduced finite intersection algebra.
I=s.groebner([u**2,v**2],u,v)
B=1+u+v
inverse=1-u-v+2*u*v
check(I.reduce(s.expand(B*inverse-1))[1]==0,'unit inverse in C[u,v]/(u^2,v^2)')
check(B.subs({u:0,v:0})!=0,'unit nonzero on reduced support')
check((u+v).subs({u:0,v:0})==0,'triple-zero nonunit control')
H=u**2*(u-1)**3
h=u*(u-1)
A=(v**2-1)*(v-2)
check(s.gcd(H,A)==1,'no common curve component')
check(s.Poly(h,u,v).total_degree()*s.Poly(A,u,v).total_degree()==6,'squarefree Bezout budget')
points=[(i,j) for i in (0,1) for j in (-1,1,2)]
check(all(H.subs(u,i)==A.subs(v,j)==0 for i,j in points),'six distinct affine intersections')
check(len(points)==6,'reduced cardinality is six')
check(s.Poly(H,u,v).total_degree()*s.Poly(A,u,v).total_degree()==15,'unreduced degree budget differs')

# Actual polynomial graph hostile: all finite-source algebra works, global etaleness fails.
for d in (3,6):
    F1=p
    F2=p*q**d-d*q
    R=u*z**d-d*z-v
    check(s.expand(R.subs({u:F1,v:F2,z:q}))==0,f'actual polynomial graph d={d}')
    check(s.diff(R,v)==-1,f'graph is a plane d={d}')
    check(s.expand(R.subs(u,0))==-d*z-v,f'boundary singleton root d={d}')
    J=s.det(s.Matrix([[s.diff(F1,p),s.diff(F1,q)],[s.diff(F2,p),s.diff(F2,q)]]))
    check(J.subs({p:1,q:1})==0,f'global Keller condition fails d={d}')
    critical=s.Poly(R.subs({u:1,v:1-d}),z)
    check(critical.eval(1)==critical.diff().eval(1)==0,f'off-boundary cover ramifies d={d}')

# Correct nonmonicity control. Chosen primitive coordinate q has minimal Z-u.
# The graph closure contains (u,v,z)=(0,1,0), but no affine source point maps there.
chosen_minimal=z-u
check(chosen_minimal.subs({u:0,v:1,z:0})==0,'bounded primitive root on nonfinite graph')
companion_relation=u*p-v
check(companion_relation.subs({u:0,v:1})==-1,'companion has no affine value')
check(s.cancel(v/u).subs(v,1)==1/u,'companion pole is not detected by chosen coordinate')

# Generic bijectivity alone does not imply singleton special fibers: nodal normalization.
nx=t**2-1
ny=t*(t**2-1)
check(s.expand(ny**2-nx**2*(nx+1))==0,'nodal normalization equation')
check(nx.subs(t,1)==nx.subs(t,-1)==ny.subs(t,1)==ny.subs(t,-1)==0,'two preimages at node')
check(s.gcd(s.diff(nx,t),s.diff(ny,t))==1,'normalization parametrization is an immersion')
check(s.cancel(ny/nx)==t,'normalization generically one-to-one')

# Euler accounting from branch-identification data, not delta invariants.
models=[('line',1,[],1),('cusp',1,[1],1),('node',1,[2],0),
        ('disjoint lines',2,[],2),('triangle',3,[2,2,2],0),('triple point',3,[3],1)]
for name,components,branches,wanted in models:
    chi=components-sum(b-1 for b in branches)
    check(chi==wanted,f'normalization Euler {name}')
    for degree in (3,6):
        omitted=(degree-1)*(1-chi)
        check(degree*(1-chi)+chi-omitted==1,f'Euler rearrangement {name},d={degree}')
check((6-1)*(1-2)<0,'disjoint-line model excluded by omitted-set nonnegativity')

out={'status':'PASS','checks':len(checks),'scope':'Independent finite controls only; universal proofs reviewed separately. No Keller realization constructed.','labels':checks}
Path(__file__).with_name('independent-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'status':'PASS'},sort_keys=True))
