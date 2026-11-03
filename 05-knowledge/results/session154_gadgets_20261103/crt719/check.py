def require(x):
 if not x: raise AssertionError("CRT check failed")
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
for p,T in ((7,48),(19,171)):
 require(pw(A,T,p)==I and all(pw(A,k,p)!=I for k in range(1,T)))
 M=pw(A,T,p*p);C=tuple(tuple(((M[i][j]-(i==j))//p)%p for j in range(2)) for i in range(2));v=(25%p,34%p);reg=Counter();total=Counter()
 for k in range(T):
  h=(v[0]+6)%p;total[h]+=1
  if row(v,C,p)[0]:reg[h,k%3]+=1
  v=row(v,A,p)
 if p==7:
  require(all(sum(reg[r,j] for j in range(3))==6 for r in range(7)))
  require(all(total[r]==(6 if r==6 else 7) for r in range(7)))
 else:
  require(all(reg[r,j]==3 for r in range(19) for j in range(3)))
for a,b in ((1,1),(2,2)):
 m=7**a*19**b;period=2736*7**(a-1)*19**(b-1);v=(25%m,34%m);fib=Counter()
 for k in range(period):
  fib[(v[0]+6)%m]+=1;v=row(v,A,m)
 require(len(fib)==m and min(fib.values())>=18)
 if a==b==1:require(all(fib[r]==(18 if r%7==6 else 21) for r in range(m)))
 print('a,b=',a,b,'mod=',m,'fibre min/max=',min(fib.values()),max(fib.values()))
print('PASS')
