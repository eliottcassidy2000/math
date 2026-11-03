def require(x):
 if not x: raise AssertionError("lift check failed")
from collections import Counter
A=((11,18),(4,7));I=((1,0),(0,1))
def mul(a,b,m):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))%m for j in range(2)) for i in range(2))
def pw(a,n,m):
 r=I
 while n:
  if n&1:r=mul(r,a,m)
  a=mul(a,a,m);n//=2
 return r
def row(v,a,m):return tuple(sum(v[i]*a[i][j] for i in range(2))%m for j in range(2))
require(pw(A,56,2197)==((1327,2028),(1183,144)))
require(pw(A,56,13)==I and all(pw(A,k,13)!=I for k in range(1,56)))
M=pw(A,56,2197);v=(25%2197,34%2197);sing=[]
for k in range(56):
 if (v[0]+6)%13==6:
  delta=(row(v,M,2197)[0]-v[0])%2197
  require(delta%169==0)
  sing.append((k,(v[0]+6)%169,delta//169))
 v=row(v,A,2197)
require(sing==[(13,45,8),(27,32,1),(41,136,5),(55,149,12)])
missing={6,19,58,71,84,97,110,123,162}
for m in (2,3,4):
 q=13**m;v=(25%q,34%q);image=set()
 for k in range(56*13**(m-1)):
  image.add((v[0]+6)%q);v=row(v,A,q)
 require(len(image)==160*13**(m-2))
 require(image=={r for r in range(q) if r%169 not in missing})
 print('13^',m,'image size',len(image))
p=19;T=171;P=pw(A,T,p*p);C=tuple(tuple(((P[i][j]-(i==j))//p)%p for j in range(2)) for i in range(2));v=(25%p,34%p);h=Counter()
for k in range(T):
 d=row(v,C,p)[0];require(d!=0)
 h[(v[0]+6)%p]+=1;v=row(v,A,p)
require(h==Counter({i:9 for i in range(19)}))
print('19-phase counts all 9; PASS')
