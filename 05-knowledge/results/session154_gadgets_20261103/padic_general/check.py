def require(x):
 if not x: raise AssertionError("p-adic check failed")
A=((11,18),(4,7)); I=((1,0),(0,1))
def mul(a,b,m):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))%m for j in range(2)) for i in range(2))
def pw(a,n,m):
 r=I
 while n:
  if n&1:r=mul(r,a,m)
  a=mul(a,a,m);n//=2
 return r
def row(v,a,m):return tuple(sum(v[i]*a[i][j] for i in range(2))%m for j in range(2))
def order(p):
 r=I
 for t in range(1,p*(p*p-1)*(p*p-p)+1):
  r=mul(r,A,p)
  if r==I:return t
 raise AssertionError
for p,T in ((7,48),(13,56),(19,171)):
 require(order(p)==T)
 M=pw(A,T,p*p)
 C=tuple(tuple(((M[i][j]-(i==j))//p)%p for j in range(2)) for i in range(2))
 v=(25%p,34%p); hist=[]
 for k in range(T):
  hist.append(((v[0]+6)%p,row(v,C,p)[0]));v=row(v,A,p)
 print('p,T,C,covered,regular:',p,T,C,len({h for h,d in hist}),len({h for h,d in hist if d}))
 if p in (7,19):require({h for h,d in hist if d}==set(range(p)))
 if p==13:
  require({h for h,d in hist}==set(range(p)))
  require({h for h,d in hist if d}==set(range(p))-{6})
  v=(25%169,34%169); hits=[]
  for k in range(T):
   if (v[0]+6)%13==6:hits.append((k,(v[0]+6)%169))
   v=row(v,A,169)
  require(hits==[(13,45),(27,32),(41,136),(55,149)])
# p=5 singular transient and exact mod-25 reachable image
v=(0,4); observed=[]
for k in range(9):observed.append((v[0]+6)%5);v=row(v,A,5)
require(observed==[1,2,4,0,3,2,4,0,3])
v=(0,9);seen=set();h=set()
while v not in seen:
 seen.add(v);h.add((v[0]+6)%25);v=row(v,A,25)
require(set(range(25))-h=={1,11,16,21})
print('PASS')
