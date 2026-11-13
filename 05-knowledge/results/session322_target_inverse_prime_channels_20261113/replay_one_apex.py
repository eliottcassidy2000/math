from itertools import product
def need(ok):
 if not ok:raise ValueError('replay failed')
def form(w,z):
 n=len(w)
 return z**(n+1)+sum(z**(n-(j-i-1)-1)*((z+1)**(j-i-1)+(z-1)**(j-i-1)) for j in range(n) for i in range(j) if w[i] and not w[j])
def literal(w):
 n=len(w)+1; N=1<<n
 A=[sum(1<<j for j in range(i+1,n-1)) + ((not w[i])<<(n-1)) for i in range(n-1)]
 A.append(sum(v<<i for i,v in enumerate(w)))
 E=[[0]*n for _ in range(N)]
 for v in range(n):E[1<<v][v]=1
 for S in range(1,N):
  for u,v in enumerate(E[S]):
   if not v:continue
   a=A[u]&~S
   while a:
    b=a&-a;a-=b
    E[S|b][b.bit_length()-1]+=v
 H=[1]+[sum(r) for r in E[1:]]
 Q=[1]+[0]*(N-1)
 out=[]
 for z in range(1,4):
  R=Q.copy()
  for S in range(1,N):
   a=S
   while a:
    R[S]+=H[a]*Q[S^a]
    a=(a-1)&S
  Q=R
  out.append(Q[-1])
 return out
for n in range(8):
 for w in product(range(2),repeat=n):
  need([form(w,z) for z in (1,2,3)]==literal(w))
for k in range(2,20):
 w=[1]+[0]*(2*k+10)
 v=w.copy();v[k+3]=v[k+4]=1
 need((form(w,2)-form(v,2))%(2**k)==0)
 need((form(w,1)-form(v,1))%4==2)
print('one-apex formula and 2-adic response separation verified')
from fractions import Fraction as F
bounds=[]
for r in range(2,7):
 def wgt(d):return 2**(r-d-1)*(3**d+1)
 bounds.append(2**(r+1)+wgt(r-2)+sum(min(wgt(j-2),wgt(r-j-1)) for j in range(2,r)))
need(bounds==[12,32,84,216,548])
s=F(1,2)+sum(F(2**j,b) for j,b in enumerate(bounds))+F(3,2)*F(2,3)**5
need(s==F(1189691,1242864)<F(24,25))
need(5**4<3**4*8 and 20**4<9**4*27 and 7*24**49<25**49)
print('serial Dirichlet bounds verified exactly')
