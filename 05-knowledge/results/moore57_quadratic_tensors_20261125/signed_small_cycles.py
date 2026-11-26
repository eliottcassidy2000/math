"""Find physical nogoods in the compressed phase equations, without singletons."""
from collections import Counter
from itertools import combinations
import sys
from antipodal_options import compressed
from slot_interface_check import BRANCHES

_,cap,bank=compressed(sys.argv[1])
C=[]
for h,Q in zip(BRANCHES[1:],bank):
    for y,*xyz,m in Q:
        if m==3 and all(cap[j][s%28]==3 for j,s in enumerate(xyz)):
            C.append((h,y,xyz))
N=[{} for _ in C]
two=[]
for i,(h,a,X) in enumerate(C):
 for j in range(i+1,len(C)):
    k,b,Y=C[j]
    if h==k:continue
    e=[(g,x%28,1^(x//28)^(y//28)) for g,(x,y) in enumerate(zip(X,Y)) if x%28==y%28]
    if len(e)>=2:
        assert {z for _,_,z in e}=={0,1}
        two.append((i,j,e))
    if len(e)==1:
        N[i][j]=e[0]; N[j][i]=e[0]
tri=[]
for i in range(len(C)):
 for j in (j for j in N[i] if j>i):
    for k in (k for k in N[i].keys()&N[j].keys() if k>j):
        F=[N[i][j],N[j][k],N[i][k]]
        if len({g for g,r,z in F})!=3:continue
        if sum(z for g,r,z in F)%2:
            tri.append(((i,j,k),F))
print(len(C),"unary-free options;",len(two),"bad parallel pairs;",len(tri),"bad triangles")
q,F=tri[0]
print("Example:",[C[i] for i in q],F)
assert len(C)==677 and len(two)==607 and len(tri)==9449
assert [C[i] for i in q]==[(3,2,[2,48,39]),(4,2,[2,35,32]),(54,18,[18,20,4])]
