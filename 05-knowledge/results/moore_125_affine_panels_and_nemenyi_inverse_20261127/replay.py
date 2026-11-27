from itertools import product
from collections import Counter
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
for p in (2,3,5):
 inds=list(product(range(p),repeat=3))
 C=Matrix([[int((z-a*x-b*y-c)%p==0) for a,b,c in inds] for x,y,z in inds])
 S=smith_normal_form(C,domain=ZZ)
 print(p,Counter(abs(S[i,i]) for i in range(len(inds))))
 assert sum(1 for i in range(len(inds)) if S[i,i])==1+p*p*(p-1)
 assert sum(1 for i in range(len(inds)) if int(S[i,i])%p)==(p+2)*(p+1)*p//6
 D=Matrix([[int((x,y)==(a,b)) for a,b,c in inds] for x,y,z in inds])
 J=Matrix.ones(len(inds))
 assert C*C.T==p*p*Matrix.eye(len(inds))+p*(J-D)
 if p==5:assert Counter(abs(S[i,i]) for i in range(len(inds)))=={0:24,1:35,5:56,25:10}
from itertools import product
allowed=[]
for c in (30,120,210):
 hits=[]
 for a,b in product(range(25),repeat=2):
  TA,TB=15*a-160,15*b-160
  cA,cB=(TA-50)//2,(TB-50)//2
  if (0<=TA<=150 and 0<=TB<=150 and
      (TA-50)%2==(TB-50)%2==0 and min(cA,cB)>=0 and cA+cB==c):
   hits.append((cA,cB,TA,TB))
 print('cyclic total',c,'allowed panels and trace pair',hits)
 if hits: allowed.append(c)
assert allowed==[30]
samples=Matrix([[y*y,y,1,0,0] for y in [0,1,2]]+[[-x*x,0,0,x,1] for x in [0,1]])
assert abs(samples.det())==2
print('five marked stress samples determinant',samples.det())
