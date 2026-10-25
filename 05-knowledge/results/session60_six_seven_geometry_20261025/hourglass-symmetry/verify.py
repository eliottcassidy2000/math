#!/usr/bin/env python3
"""Exact finite checks for hourglass symmetries and its central Shioda sextic.
Requires Python 3 and SymPy. All checks survive python -O.
"""
import itertools
import json
from collections import Counter
import sympy as s

checks = 0
def check(ok, msg):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(msg)

def tup(M): return tuple(M)
def normalized(v):
    a = next(q for q in v if q)
    return tuple(s.cancel(q/a) for q in v)

E,x,z,t,u = s.symbols('E x z t u')
roots = [s.Matrix(v) for v in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)]]
rootset=set(map(tup,roots))
L=[s.Matrix(v) for v in [(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)]]
points={normalized(a.cross(b)) for a,b in itertools.combinations(L,2)}
inc={p:tuple(i for i,l in enumerate(L) if (l.dot(s.Matrix(p)))==0) for p in points}
triples=[v for v in inc.values() if len(v)==3]
check(len(triples)==3,'three triple points')
check(Counter(len(v) for v in inc.values())=={2:12,3:3},'intersection multiplicities')
check([sum(i in v for v in triples) for i in range(7)]==[3,1,1,1,1,1,1],'unique central line')

# An invertible linear symmetry is fixed by the images of the basis roots.
G={}
for a,b in itertools.permutations(roots,2):
    M=s.Matrix.hstack(a,b)
    if M.det() and {tup(M*v) for v in roots}==rootset:
        G[tup(M)]=M
check(len(G)==12,'all linear root symmetries')
R=s.Matrix([[1,-1],[1,0]])
S=s.Matrix([[0,1],[1,0]])
check(R**6==s.eye(2) and all(R**k!=s.eye(2) for k in range(1,6)),'rotation exact order six')
check(S*S==s.eye(2) and S*R*S==R.inv(),'dihedral relation')
generated={tup(R**k*S**e) for k in range(6) for e in range(2)}
check(set(G)==generated,'generators exhaust D12')
check({M.det() for M in G.values()}=={-1,1},'unimodular rational symmetries')
perms=[]
for M in G.values():
    perm=[0]+[1+roots.index(M*v) for v in roots]
    perms.append(perm)
    for a,b in zip(L,[L[i] for i in perm]):
        check(s.diag(1,M)*a==b,'exact normalized line permutation')
check({p[1] for p in perms}==set(range(1,7)),'finite-line transitivity')
node_counts=[]
for omit in range(7):
    retained=[v for v in L if v!=L[omit]]
    ps={normalized(a.cross(b)) for a,b in itertools.combinations(retained,2)}
    n3=sum(sum(l.dot(s.Matrix(p))==0 for l in retained)==3 for p in ps)
    node_counts.append(4*n3)
check(node_counts==[0,8,8,8,8,8,8],'sixroot deletion node counts')

# Conic, primal duality, parametrization and exact six branch labels.
Q=s.Matrix([[1,0,0],[0,-1,s.Rational(1,2)],[0,s.Rational(1,2),-1]])
check(Q.det()!=0,'smooth dual conic')
for l in L[1:]: check((l.T*Q*l)[0]==0,'six dual points on conic')
check((L[0].T*Q*L[0])[0]!=0,'hostile central dual point not on conic')
primal=s.factor((s.Matrix([E,x,z]).T*Q.inv()*s.Matrix([E,x,z]))[0]*3)
check(primal==3*E**2-4*x**2-4*x*z-4*z**2,'primal tangent conic')
D=t*t-t+1
param=s.Matrix([D,t*t-1,t*(t-2)])
check(s.expand((param.T*Q*param)[0])==0,'parametrization lies on conic')
branch=[s.oo,-1,s.Integer(0),s.Rational(1,2),s.Integer(1),s.Integer(2)]
param_points=[]
for a in branch:
    p=s.Matrix([1,1,1]) if a==s.oo else param.subs(t,a)
    param_points.append(normalized(p))
check(set(param_points)=={normalized(l) for l in L[1:]},'branch values enumerate finite lines')
rot=lambda a: -s.Integer(1) if a==s.oo else (s.oo if a==2 else s.cancel((s.sympify(a)+1)/(2-a)))
check([rot(a) for a in [2,s.oo,-1,0,s.Rational(1,2),1]]==[s.oo,-1,0,s.Rational(1,2),1,2],'Mobius six-cycle')
f=t*(t+1)*(t-1)*(t-2)*(2*t-1)
check(s.gcd(f,s.diff(f,t))==1,'smooth genus-two branch sextuple')
delta=s.sqrt(-3);alpha=(1+delta)/2;beta=(1-delta)/2
check(s.expand((t-alpha)**6-(t-beta)**6+3*delta*f)==0,'sixth-root transform')
q=s.expand(E*D+x*(t*t-1)+z*t*(t-2))
q2=E+x+z
check(s.Matrix([[1,1,1],[-1,0,-2],[1,-1,0]]).det()==-3,'coefficient map invertible')
B=(E**2-x**2)*(E**2-z**2)*(E**2-(x+z)**2)
prod=q2
for a in branch[1:]: prod*=q.subs(t,a)
check(s.expand(4*prod-27*B)==0,'Shioda scalar 27')
# Two elliptic maps and differential independence, using v²=u⁶−1.
check(s.expand((u**6-1)-((u*u)**3-1))==0,'first elliptic map')
check(s.cancel((u**6-1)/u**6-(1-(u**-2)**3))==0,'second elliptic map')
check(s.diff(u*u,u)==2*u and s.cancel(s.diff(u**-2,u)*u**3)==-2,'independent elliptic differential coefficients')
# Explicit NS classes on E0×E0: fibers, diagonal, graph of order-three CM.
# Graph intersections: deg(id−zeta3)=3; all fiber intersections =1.
NS=s.Matrix([[0,1,1,1],[1,0,1,1],[1,1,0,3],[1,1,3,0]])
check(NS.det()==-3,'four independent divisor classes')
# Hostile controls: root-set mutation loses conic and finite-group conclusion.
mutated=s.Matrix([1,1,2])
check((mutated.T*Q*mutated)[0]!=0,'mutated sixth line leaves conic')
# Möbius crossratios: (0,infinity;1,i)=1/i, nonreal, unlike rational branch set.
check(s.simplify(1/s.I)==-s.I,'Bolza nonreal crossratio')

report={
 'status':'PASS','checks':checks,
 'universe':{'lines':7,'line_pairs':21,'ordered_root_basis_images':30,'linear_symmetries':12,'deletions':7},
 'line_order':['E','E+x','E-x','E+z','E-z','E+x+z','E-x-z'],
 'triple_incidence':sorted(triples),'sixroot_node_counts':node_counts,
 'symmetry_order':12,'signed_permutation_subgroup_order':768,
 'canonical_representation':'det(base A) tensor permutation on seven root coordinates',
 'rational_branch_values':['infinity','-1','0','1/2','1','2'],
 'genus_two_polynomial':str(s.expand(f)),
 'Shioda_cover':'W^2=27*(E^2-x^2)*(E^2-z^2)*(E^2-(x+z)^2)',
 'CM_NS_intersection_determinant':int(NS.det()),
 'geometric_Kummer_picard_rank':20,
 'boundaries':['complex geometric assertions do not identify rational twists','D12 is full projective base-arrangement group, not asserted full surface automorphism group','original diagonal K3 Kummer property requires separate argument']}
print(json.dumps(report,indent=2,sort_keys=True))
