from itertools import combinations, permutations
from math import factorial

def ck(ok,msg):
 if not ok:raise AssertionError(msg)

def star(I,S,w,cut,u):
 n=len(I);v=n+u+w;mask=(1<<(v+1))-1; A=[]
 for i in range(v):
  if i<n:
   A.append(I[i]|(((1<<v)-1)^((1<<n)-1))|(((S>>i&1)^1)<<v))
  else:
   A.append((((1<<v)-1)^((1<<(i+1))-1))|((i<n+u or not(cut>>(i-n-u)&1))<<v))
 A.append(S|(cut<<(n+u)));return tuple(A)

def profiles(A,full):
 n=len(A); F=(1<<n)-1; H=[0]*(F+1);E=[[0]*n for _ in H];H[0]=1
 for i in range(n):E[1<<i][i]=1
 for s in range(1,F+1):
  for i,c in enumerate(E[s]):
   if not c:continue
   H[s]+=c;a=A[i]&~s
   while a:
    z=a&-a;a-=z;E[s|z][z.bit_length()-1]+=c
 if not full:return H[-1]
 Q=[0]*(F+1)
 for s in range(F+1):
  t=s
  while True:
   Q[s]+=H[t]*H[s^t]
   if not t:break
   t=(t-1)&s
 return H[-1],Q[-1],sum(H[t]*Q[F^t] for t in range(F+1))

def strong(A):
 J=(1<<len(A))-1
 for B in (A,tuple(J^a^(1<<i) for i,a in enumerate(A))):
  S=1
  while True:
   T=S
   for i,x in enumerate(B):
    if S>>i&1:T|=x
   if T==S:break
   S=T
  if S!=J:return False
 return True

def slots(A):
 for p in permutations(range(len(A))):
  for k in range(len(A)+1):
   a,b=p[:k],p[k:]
   if all(A[i]>>j&1 for x in (a,b) for i,j in zip(x,x[1:])):
    yield a,b

def L(n):
 if not n:return 0
 t,b=divmod(n,2)
 return L(t)+b*(t+(t%2))

def onto(m,j):
 from math import comb
 return sum((-1)**(j-a)*comb(j,a)*a**m for a in range(j+1))

def channels(Y):
 c=[]
 for m,y in enumerate(Y,1):
  r=y-sum(x*onto(m,j) for j,x in enumerate(c,1))
  ck(r%factorial(m)==0,(Y,m,r))
  c.append(r//factorial(m))
 return c

def to_base(c):
 from math import comb
 return [sum(x*(-1)**(j-a)*comb(j,a) for j,x in enumerate(c,1) if j>=a)
         for a in range(1,len(c)+1)]

def run():
 I=(10,0,3,6)
 l=[(p,q) for p,q in slots(I) if q and (3>>q[0]&1)]
 ck((profiles(I,False),len(l),sum(bool(p and 3>>p[-1]&1) for p,q in l))==(3,9,4),"I")
 families=[(I,3,4,5,99),((0,),1,5,7,25)]
 for I,S,w,C,A in families:
  W=tuple(((1<<w)-1)^((1<<(i+1))-1) for i in range(w))
  b=sum(1 for p,q in slots(W) if not p or not (C>>p[-1]&1))
  Wv=tuple(x|(((C>>i&1)^1)<<w) for i,x in enumerate(W))+(C,)
  d=profiles(Wv,False)
  ck((b,d)==((11,9) if A==99 else (25,23)),"W")
  Y=[]
  for m in range(1,5):
   M=star(I,S,w,C,m)
   ck(strong(M),("strong",A,m))
   p=profiles(M,True)
   ck(p[0]==A*2**m+23,("rail",A,m,p))
   Y.append(p)
  pred=[to_base(channels([p[k] for p in Y[:k+2]])) for k in (0,1,2)]
  print("rail",A,"profiles in bases 1,2,3,4",pred)
  if A==99: wanted99=pred
  else:wanted25=pred
 ck(wanted99==[[23,99],[-330,1696,3630],[0,-7755,35829,46155]],wanted99)
 ck(wanted25==[[23,25],[-58,592,334],[0,-609,5373,2073]],wanted25)
 for a in (25,99):
  print("odd slope",a,"intercepts",[
      (x,a//x,a//x-L(a//x-1)-(a//x-1)%2-L(x))
      for x in range(1,a+1) if a%x==0])
if __name__=="__main__":run()
