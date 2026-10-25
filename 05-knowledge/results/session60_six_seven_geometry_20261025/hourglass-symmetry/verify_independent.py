#!/usr/bin/env python3
"""Standard-library independent Möbius audit of the six branch values."""
from fractions import Fraction as F
from itertools import permutations
import json

checks=0
def require(v,why):
    global checks
    checks+=1
    if not v: raise RuntimeError(why)

def det(p,q): return p[0]*q[1]-p[1]*q[0]
def norm(v):
    a=next(x for x in v if x)
    return tuple(F(x,a) for x in v)
def apply(m,p): return norm((m[0]*p[0]+m[1]*p[1],m[2]*p[0]+m[3]*p[1]))
def mul(m,n):
    return norm((m[0]*n[0]+m[1]*n[2],m[0]*n[1]+m[1]*n[3],m[2]*n[0]+m[3]*n[2],m[2]*n[1]+m[3]*n[3]))
def power(m,k):
    a=(F(1),F(0),F(0),F(1))
    for _ in range(k):a=mul(a,m)
    return a
B=[norm(p) for p in [(1,0),(-1,1),(0,1),(1,2),(1,1),(2,1)]]
S=set(B); group={}
# Map infinity,-1,0 to P,Q,R. Matrix columns aP,bR with a=det(R,Q), b=det(P,Q).
for P,Q,R in permutations(B,3):
    a=det(R,Q);b=det(P,Q)
    M=norm((a*P[0],b*R[0],a*P[1],b*R[1]))
    require(M[0]*M[3]-M[1]*M[2]!=0,'unique map invertible')
    require([apply(M,p) for p in B[:3]]==[P,Q,R],'prescribed triple')
    images=[apply(M,p) for p in B]
    if set(images)==S: group[M]=[B.index(p) for p in images]
require(len(group)==12,'full reduced automorphism group order12')
R=norm((1,1,-1,2)); Sref=norm((-1,1,0,1)); I=norm((1,0,0,1))
require(power(R,6)==I and all(power(R,k)!=I for k in range(1,6)),'rotation order6')
require(power(Sref,2)==I,'reflection order2')
require(mul(mul(Sref,R),Sref)==power(R,5),'dihedral relation')
generated={mul(power(R,k),power(Sref,e)) for k in range(6) for e in range(2)}
require(set(group)==generated,'all stabilizers generated')
# Crossratios are independently all rational; retain their complete finite palette.
crs=set()
for a,b,c,d in permutations(B,4):
    cr=F(det(a,c)*det(b,d),det(a,d)*det(b,c))
    crs.add(cr)
require(all(cr not in (0,1) for cr in crs),'crossratios avoid zero and one')
# Distinct roots, rotation cycle and generic quadratic reconstruction.
require(len(set(B))==6,'six branch points distinct')
# Coefficient matrix inverse: E=(2q0+2q2+q1)/3, x=(-q0+2q2+q1)/3, z=(-q0-q2-2q1)/3.
for E in range(-2,3):
 for x in range(-2,3):
  for z in range(-2,3):
   q2=E+x+z;q1=-E-2*z;q0=E-x
   require((2*q0+2*q2+q1, -q0+2*q2+q1,-q0-q2-2*q1)==(3*E,3*x,3*z),'independent coefficient inverse')
report={'status':'PASS','checks':checks,'ordered_mobius_image_triples':120,'full_reduced_branch_group_order':len(group),'crossratio_palette':sorted(map(str,crs)), 'crossratio_count':len(crs),'group_permutations':sorted(group.values()),'uses_sympy':False}
print(json.dumps(report,indent=2,sort_keys=True))
