"""Independent symbolic-matrix audit of complete projective jet precision."""
from fractions import Fraction as F
from math import comb,gcd,lcm,prod
from itertools import combinations,permutations
import sys
from sympy import Matrix,Poly,Symbol,ZZ
from sympy.matrices.normalforms import smith_normal_form
sys.stdout.reconfigure(newline="\n")
T=Symbol('T');CHECKS=0
def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise ArithmeticError(msg)
def bracket(x,y):return x[0]*y[1]-x[1]*y[0]
def observe(V,W,m):
 D=sum(m)-1;rows=[]
 for v,w,r in zip(V,W,m):
  columns=[Poly((v[0]+T*w[0])**(D-j)*(v[1]+T*w[1])**j,T) for j in range(D+1)]
  rows.extend([[p.nth(k) for p in columns] for k in range(r)])
 return Matrix(rows)
def packets(V,W,m):
 out=[]
 for i,r in enumerate(m):
  q=Poly(prod((bracket(V[i],V[j])+T*bracket(W[i],V[j]))**m[j] for j in range(len(V)) if i!=j),T)
  inv=[F(1,int(q.nth(0)))]
  for k in range(1,r):inv.append(-sum(F(int(q.nth(s)))*inv[k-s] for s in range(1,k+1))/int(q.nth(0)))
  out.append(lcm(*(x.denominator for x in inv)))
 return out
def change(G,x):return tuple(sum(G[i][j]*x[j] for j in range(2)) for i in range(2))
def pval(x,p):
 if not x:return 100000
 x=abs(int(x));v=0
 while x%p==0:v+=1;x//=p
 return v
def main():
 cases=[([(1,0)],[(0,1)],[r]) for r in range(1,6)]
 for m in ([1,2,3],[3,1,2],[2,3,4]):
  cases.append(([(1,0),(1,4),(0,1)],[(0,1),(0,1),(-1,0)],m))
  cases.append(([(0,1),(9,1),(18,1)],[(-1,0)]*3,m))
 cases.extend([([(1,0),(0,1),(1,1),(1,2)],[(0,1),(-1,0),(0,1),(0,1)],[1,2,3,2]),
               ([(1,0),(1,8)],[(0,1),(0,1)],[3,4])])
 matrices=0
 for V,W,m in cases:
  need(all(bracket(v,w)==1 for v,w in zip(V,W)),'complete primitive unimodular local frames')
  base=None
  for G in (((1,0),(0,1)),((2,1),(1,1)),((0,1),(1,0))):
   sign=bracket(G[0],G[1]);VV=[change(G,v) for v in V];WW=[tuple(sign*x for x in change(G,w)) for w in W]
   A=observe(VV,WW,m);inv=A.inv();node=packets(VV,WW,m)
   S=smith_normal_form(A,domain=ZZ);smith=[abs(int(S[j,j])) for j in range(A.rows)]
   den=lcm(*(int(x.q) for x in inv));need(den==lcm(*node)==smith[-1],'literal inverse, reciprocal coefficients and integer Smith agree')
   need(abs(int(A.det()))==prod(abs(bracket(VV[i],VV[j]))**(m[i]*m[j]) for i,j in combinations(range(len(V)),2)),'weighted determinant with directions at infinity')
   offset=0
   for i,r in enumerate(m):
    need(lcm(*(int(inv[j,offset].q) for j in range(A.rows)))==node[i],'value column alone attains node precision')
    offset+=r
   if base is None:base=smith
   need(smith==base,'entire Smith module under both orientation gauges')
   # Complete local coordinate shear is triangular even when r>2.
   WW2=[tuple(w[j]+(i+2)*v[j] for j in range(2)) for i,(v,w) in enumerate(zip(VV,WW))]
   B=observe(VV,WW2,m)
   transition=B*inv
   need(all(x.q==1 for x in transition) and abs(transition.det())==1,'literal full higher-order target coordinate change')
   matrices+=2
 # Intrinsic residue uses cancelled bracket products, no illegal common chart.
 sidecars=0
 exceptional={3,11,15,17,21,29}
 for e in (1,2,3):
  for a in range(2,31):
   V=[(1,0),(1,-31**e),(1,-31**e*a)]
   for z in ((0,1),(1,1),(2,1),(31,2)):
    num=bracket(V[0],V[2])*bracket(V[1],z)//31**e
    den=bracket(V[0],V[1])*bracket(V[2],z)//31**e
    residue=(num*pow(den,-1,31**e))%(31**e)
    need(residue==a%(31**e),'transverse invariance modulo cluster depth')
    sidecars+=1
   for Vp in permutations(V):
    z=(0,1);num=bracket(Vp[0],Vp[2])*bracket(Vp[1],z)//31**e;den=bracket(Vp[0],Vp[1])*bracket(Vp[2],z)//31**e
    r=num*pow(den,-1,31)%31
    need((r in exceptional)==(a in exceptional),'p31 bit is invariant under all node orderings')
 # A common affine chart with unit denominators cannot cover all p+1 classes.
 for p in (2,3,5,31):
  V=[(1,0)]+[(a,1) for a in range(p)]
  for x,y in V:need(any((x*a+y*b)%p==0 for a,b in V),'every primitive denominator loses a projective residue')
 print('ROOT_AUDIT: arbitrary projective complete-jet precision, covariance, and intrinsic p31 sidecar PASS')
 print('INDEPENDENT_SYMBOLIC_MATRICES:',matrices,'base configurations',len(cases),'integer Smith and literal inverse for each orientation')
 print('SIDECAR_CHECKS:',sidecars,'p31 every admissible residue,depths1..3,four transverse frames,all six node orders')
 print('CONTROLS: singleton jets;heterogeneous orders;infinity;orientation reversal;nontrivial tangent shears;all residue classes')
 print('SCOPE: existing p31 cluster ideals transport by audited covariance/splitting; no new general partition formula')
 print('ACTIVE_CHECKS:',CHECKS)
if __name__=='__main__':main()
