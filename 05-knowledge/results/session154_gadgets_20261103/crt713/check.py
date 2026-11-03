def require(x):
 if not x: raise AssertionError("joint image check failed")
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
for m,T in ((7,48),(169,728)):
 require(pw(A,T,m)==I and all(pw(A,k,m)!=I for k in range(1,T)))
require(pw(A,4368,1183)==I)
v=(25%1183,34%1183);fib=Counter()
for k in range(4368):
 h=(v[0]+6)%1183;fib[h%7,h%169]+=1;v=row(v,A,1183)
require(len(fib)==934)
def allowed(u,v):
 if u==6:return v%13 in (4,9,11,12)
 return v%13!=12 and (v%13!=6 or v in (32,45,136,149))
require(all(((u,v) in fib)==allowed(u,v) for u in range(7) for v in range(169)))
require(Counter(fib.values())=={2:234,3:78,4:156,5:156,6:247,12:39,13:24})
for c in (2,4):
 image={(3*u+c)%7 for u,v in fib} # simple marginal sanity check
 require(image==set(range(7)))
print('mod 7*169 images:',len(fib),'of',7*160,'independent possibilities')
print('row sizes:',[sum(1 for v in range(169) if (u,v) in fib) for u in range(7)])
print('fibre histogram:',sorted(Counter(fib.values()).items()))
print('PASS')
