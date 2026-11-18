"""Run from alongside the original blue_core packet, or standalone."""
from itertools import permutations, combinations
from math import comb, factorial, isqrt
from fractions import Fraction

def dp(A):
 n=len(A); N=1<<n; V=[[0]*n for _ in range(N)]
 for i in range(n):V[1<<i][i]=1
 for s in range(1,N):
    for i,x in enumerate(V[s]):
      if not x:continue
      u=A[i]&~s
      while u:
        b=u&-u;u-=b;V[s|b][b.bit_length()-1]+=x
 return V,[1]+[sum(x) for x in V[1:]]

def transpose(A):
 return [sum((A[j]>>i&1)<<j for j in range(len(A))) for i in range(len(A))]
def tt(n):return [(1<<n)-(1<<(i+1)) for i in range(n)]
def ordinal(*parts):
 out=[];n=sum(map(len,parts))
 for A in parts:
    off=len(out)
    for row in A:out.append(row<<off | ((1<<n)-(1<<(off+len(A)))))
 return out

def cyclic(X,m):
 Y=transpose(X);n=len(X)
 return [r|(((1<<n)-1)<<n) for r in X]+[(r<<n)|(((1<<m)-1)<<(2*n)) for r in Y]+[(r<<2*n)|((1<<n)-1) for r in tt(m)]

def wrapper(X,r,m):
 n=len(X); a=ordinal(X,tt(m),transpose(X));s=n+m+r
 a[r]^=1<<s;a[s]^=1<<r
 return a

def profile(A):
 _,h=dp(A); N=len(h);p=[1]+[0]*(N-1);q=[p[-1]]
 for j in range(1,len(A)+1):
   Q=[0]*N
   for s in range(N):
     t=s
     while t:
       Q[s]+=h[t]*p[s^t]
       t=(t-1)&s
   p=Q;q.append(p[-1])
 return q

def val(q,z):return sum(comb(z,j)*q[j] for j in range(min(z,len(q)-1)+1))

def run(q,m):
 ans=0
 for L in range(3,3*m+3):
    for b in range(3):
       a=[0]*3
       for i in range(L):a[(i+b)%3]+=1
       x,y,z=a
       if z>m:continue
       c=sum((-1)**(z-j)*comb(z,j)*j**m for j in range(z+1))
       ans+=q[x]*q[y]*c
 return ans

def poly_eval(known,n,z):
 d=n//2
 if not d:return z**n
 value=Fraction(z**n)
 for i in range(1,d+1):
    t=Fraction(known[i]-i**n,i**(2-n%2))*z**(2-n%2)
    for j in range(1,d+1):
      if i!=j:t*=Fraction(z*z-j*j,i*i-j*j)
    value+=t
 return value

A=(58,4,1,6,14,30);B=(38,8,2,5,15,30)
def cycles(X,m):
 out=[]
 for S in combinations(range(len(X)),m):
    for p in permutations(S[1:]):
      w=(S[0],)+p
      if all(X[i]>>j&1 for i,j in zip(w,w[1:]+w[:1])):out.append(set(w))
 return out

if __name__=='__main__':
 for X in (A,B,tt(6),[2,4,1],tt(1),tt(2)):
    q=profile(X);n=len(X);d=n//2
    P=[val(q,z) for z in range(13)]
    if any(poly_eval(P,n,z)!=p for z,p in enumerate(P)):raise ArithmeticError('interpolation')
    for m in range(1,max(2,d)):
      padded=q+[0]*(m+2)
      N=dp(cyclic(X,m))[1][-1]
      if N!=run(padded,m):raise ArithmeticError('run')
      a=padded[:];a[m+1]=0
      dlt,res=divmod(N-run(a,m),factorial(m))
      if res or isqrt(dlt+q[m]**2)-q[m]!=padded[m+1]:raise ArithmeticError('inverse')
 for X in (A,B):
    n=len(X); E,h=dp(X); V,_=dp(transpose(X))
    q=[sum(V[S][r]*h[(1<<n)-1-S] for S in range(1<<n)) for r in range(n)]
    out=[]
    for r in range(n):
        Xr=wrapper(X,r,0); W,hW=dp(Xr); Wt,_=dp(transpose(Xr))
        if any(W[-1][i]==0 or Wt[-1][i]==0 for i in range(2*n)):raise ArithmeticError('missing Hamilton starts/ends')
        p=hW[-1];out.append(p)
        if p!=h[-1]**2+q[r]**2-E[-1][r]**2:raise ArithmeticError('marked square at zero')
        p1=dp(wrapper(X,r,1))[1][-1]
        if p1!=h[-1]**2+2*q[r]**2:raise ArithmeticError('marked square at one')
        out[-1]=(p,p1)
        s=[i+n if i<n else i-n for i in range(2*n)]
        if any((Xr[i]>>j&1)!=(Xr[s[j]]>>s[i]&1) for i in range(2*n) for j in range(2*n)):raise ArithmeticError('anti')
    c3=cycles(X,3); c5=cycles(X,5)
    print('profile',profile(X),'triangles',len(c3),'five-cycles',len(c5),
          'disjoint triangle pairs',sum(not (S&T) for S,T in combinations(c3,2)),
          'marked',q,'ends',E[-1],'probes',out)
