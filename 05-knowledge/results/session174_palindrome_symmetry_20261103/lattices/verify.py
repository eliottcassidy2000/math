from itertools import combinations, product
from fractions import Fraction as F
import json

def check(c):
    if not c: raise RuntimeError('Exact lattice check failed')
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def shell(n):return {tuple(1 if i in S else -1 for i in range(n)) for S in combinations(range(n),n//2)}
def reflect(x,r):
    c=F(2*dot(x,r),dot(r,r));return tuple(F(a)-c*b for a,b in zip(x,r))
def det(A):
    A=[[F(x) for x in row] for row in A];d=F(1)
    for j in range(len(A)):
        k=next(k for k in range(j,len(A)) if A[k][j])
        if k!=j:A[j],A[k]=A[k],A[j];d=-d
        p=A[j][j];d*=p
        for k in range(j+1,len(A)):
            c=A[k][j]/p;A[k]=[a-c*b for a,b in zip(A[k],A[j])]
    return d
H=shell(8)
I={tuple(2*(k==i)-2*(k==j) for k in range(8)) for i in range(8) for j in range(8) if i!=j}
R=H|I;check((len(H),len(I),len(R))==(70,56,126))
E=set()
for i,j in combinations(range(8),2):
    for a,b in product((-2,2),repeat=2):E.add(tuple(a*(k==i)+b*(k==j) for k in range(8)))
E|={s for s in product((-1,1),repeat=8) if sum(x<0 for x in s)%2==0}
check(len(E)==240 and {r for r in E if sum(r)==0}==R)
check(all(reflect(x,r) in R for x in R for r in R))
one=H|{reflect(x,r) for x in H for r in H};check(one==R)
basis=[tuple(2*(k==i)-2*(k==i+1) for k in range(8)) for i in range(6)]+[(-1,-1,-1,1,1,1,1,-1)]
G=[[F(dot(a,b),4) for b in basis] for a in basis];check(det(G)==2)
check({(i,j) for i in range(7) for j in range(i+1,7) if G[i][j]==-1}=={(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)})
for n,count in ((2,2),(4,6)):
    S=shell(n);check(len(S)==count and all(reflect(x,r) in S for x in S for r in S))
check(all(F(2)-F(8,n) not in (-2,-1,0,1,2) for n in range(6,42,2) if n!=8))
print(json.dumps({'status':'PASS','H8':len(H),'integer_completion':len(I),'E7_closure':len(R),'E8_ambient':len(E),'reflection_tests':len(R)**2,'one_step_completion':len(one),'simple_gram':[[int(x) for x in row] for row in G],'gram_determinant':int(det(G)),'adjacent_Cartan_coefficients':{str(n):str(F(2)-F(8,n)) for n in range(2,42,2)}},indent=2))
