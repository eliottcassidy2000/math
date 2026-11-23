from functools import cache
from math import comb,factorial

def transitive(n):
 return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def join(A,B):
 n,m=len(A),len(B)
 return tuple(x|(((1<<m)-1)<<n) for x in A)+tuple(x<<n for x in B)
def direct(A):
 n=len(A); d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,len(d)):
  for i,c in enumerate(d[S]):
   if not c:continue
   z=A[i]&~S
   while z:
    b=z&-z;z-=b
    d[S|b][b.bit_length()-1]+=c
 H=[1]+[sum(z) for z in d[1:]]
 return H[-1],sum(h*H[(len(H)-1)^s] for s,h in enumerate(H))
def wrap(A):
 n=len(A)
 return tuple(x|1<<n for x in A)+(1<<(n+1),(1<<n)-1)

def op(f):
 G=[0]*(len(f)+2)
 for j,a in enumerate(f):
  G[j+2]+=a
  for i in range((j-1)%2,j,2):G[i+1]+=2*a*comb(j,i)
 return G
def val(f,t):return sum(a*t**i for i,a in enumerate(f))

@cache
def tower(m,k):
 if m not in (1,2) or k<0:raise ValueError
 return tuple(op(tower(m,k-1))) if k else (0,)*m+(1,)

def encode(m,k,s):return val(tower(m,k),1),val(tower(m,k),2)<<s
def decode(D,Q):
 if D<3 or Q<2*D:return
 m=1 if D%3==0 else 2;k=1
 while val(tower(m,k),1)<D:k+=1
 f=tower(m,k)
 if val(f,1)!=D:return
 q,r=divmod(Q,val(f,2))
 if r or q.bit_count()!=1:return
 return m,k,q.bit_length()-1

def exterior(A):
 n=len(A)
 return tuple(x|(3<<n) for x in A)+(6<<n,8<<n,((1<<n)-1)|(10<<n),((1<<n)-1)|(1<<n))
 # b0->b1,c0; b1->c1; c0->b1,c1; c1->b0.

@cache
def S(n,k):
 if k==0:return int(n==0)
 return 0 if n<k else S(n-1,k-1)+k*S(n-1,k)
def coeff(slot):
 ans=[]
 for m in range(1,slot+3):
  y=direct(exterior(transitive(m)))[slot-1]
  y-=sum(a*factorial(j)*S(m,j) for j,a in enumerate(ans,1))
  b=y//factorial(m)
  assert y==b*factorial(m) and b>0
  ans.append(b)
 return tuple(ans)

MA=(42,124,345,432,481,68,137,295,99)
MB=(266,412,121,464,417,203,403,261,36)
if __name__=='__main__':
 for m in (1,2):
  A=transitive(m)
  for k in range(1,5):
   A=wrap(A)
   for s in range(3):
    assert direct(join(transitive(s),A))==encode(m,k,s)
    assert decode(*encode(m,k,s))==(m,k,s)
 print('tower two-reading inverse agrees with subset DP')
 print('minimum-kappa-two shell P-coefficients',coeff(1),coeff(2))
 print('two hidden-state successors',direct(exterior(MA)),direct(exterior(MB)))
