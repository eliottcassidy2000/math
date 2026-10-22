import sympy as s
import itertools, json
from pathlib import Path
out={}
# Two explicitly scissors-congruent boxes: cut unit cube into y-halves,
# translate one half by (1,-1/2,0), then center resulting 2 x 1/2 x 1 box.
a,b=s.symbols('a b', positive=True)
pbox=a**4/s.Integer(80)+b**4/s.Integer(80)-a*a*b*b/s.Integer(24)
assert pbox.subs({a:1,b:1})==-s.Rational(1,60)
assert pbox.subs({a:2,b:s.Rational(1,2)})==s.Rational(611,3840)
out['equal_volume_Dehn_zero_boxes_harmonic_integrals']=[str(pbox.subs({a:1,b:1})),str(pbox.subs({a:2,b:s.Rational(1,2)}))]
x,y,z=s.symbols('x y z'); p=x**4-6*x*x*y*y+y**4
assert s.expand(sum(s.diff(p,v,2) for v in (x,y,z)))==0
# Planar harmonic frame hostile: Q rotates vectors by 45 degrees and flips A.
H=[s.diag(1,-1)/s.sqrt(2),s.Matrix([[0,1],[1,0]])/s.sqrt(2)]
X=[s.Matrix(v) for v in [(1,0),(-1,0),(0,1),(0,-1)]]
Q=s.Matrix([[1,-1],[1,1]])/s.sqrt(2)
def op(points,basis,n,m):
    M=s.zeros(len(basis))
    for v in points:
        f=s.Matrix([(v.T*h*v)[0] for h in basis]); M+=f*f.T
    c=s.Rational(2*len(points),n*(n+2))*m*m
    return s.simplify(M-c*s.eye(len(basis)))
A=op(X,H,2,1); B=op([Q*v for v in X],H,2,1)
assert A==s.diag(1,-1) and B==-A
assert A.eigenvals()==B.eigenvals()
assert (A+B)==s.zeros(2) and (2*A).det()==-4
out['spectral_addition_hostile']={'A':str(A),'rotated_A':str(B),'same_spectrum':True,'aligned_sum_trace_square':str(s.trace((A+A)**2)),'misaligned_sum_trace_square':str(s.trace((A+B)**2))}
# Exact Gram-cycle formula tested through k=5 on the 3D cross polytope.
X=[sign*s.eye(3)[:,i] for i in range(3) for sign in (-1,1)]
H=[s.diag(1,-1,0)/s.sqrt(2),s.diag(1,1,-2)/s.sqrt(6)]
for i,j in [(0,1),(0,2),(1,2)]:
    h=s.zeros(3); h[i,j]=h[j,i]=1/s.sqrt(2); H.append(h)
A=op(X,H,3,1); c=s.Rational(4,5); M=A+c*s.eye(5)
G=s.Matrix([[(v.dot(w))**2-s.Rational(1,3) for w in X] for v in X])
assert s.trace(A)==0
for k in range(1,6):
    cyc=sum(s.prod(G[t[i],t[(i+1)%k]] for i in range(k)) for t in itertools.product(range(6),repeat=k))
    assert cyc==s.trace(M**k)==s.trace(G**k)
V=sum(v.dot(w)**4 for v in X for w in X)-s.Rational(3*6**2,3*5)
assert s.trace(A*A)==V==s.Rational(24,5)
out['gram_cycles_k_1_to_5']=[str(s.trace(M**k)) for k in range(1,6)]
out['cross_polytope_angular_variance']=str(V)
out['status']='PASS; exact SymPy rational/algebraic arithmetic'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
# Regular tetrahedral edge-direction harmonic sum, exact.
tetra=[s.Matrix(t) for t in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
dirs=[]
for i,j in itertools.combinations(range(4),2):
    v=tetra[i]-tetra[j]; dirs.append(v/s.sqrt(v.dot(v)))
assert s.simplify(sum(p.subs({x:v[0],y:v[1],z:v[2]}) for v in dirs))==-1
out['tetrahedron_harmonic_edge_sum']='-1'
# Irrational-angle finite configuration nevertheless lies in square lattice.
u=s.Matrix([1,2]); v=s.Matrix([2,1])
assert u.dot(v)/s.sqrt(u.dot(u)*v.dot(v))==s.Rational(4,5)
out['square_lattice_equal_norm_pair_cosine']='4/5'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
# Directed 3-cycle skew spectrum is not a real root-of-unity argument.
S=s.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
l=s.symbols('l')
assert s.expand(S.charpoly(l).as_expr()-l*(l*l+3))==0
out['directed_three_cycle_skew_charpoly']='l*(l^2+3)'
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
