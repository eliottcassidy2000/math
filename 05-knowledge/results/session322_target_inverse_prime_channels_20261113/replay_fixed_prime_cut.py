from math import gcd
from replay_dyadic_inverse import BANK,info
A,B=4595601,1915
def bval(n,s):
 b=0
 while b<s and n%3==0:b+=1;n//=3
 return b
def compile_target(n):
 if n<=0 or not n&1:return None
 b=bval(n,n.bit_length());u=n//3**b;m=(u//A).bit_length()
 if m<1:return None
 d,r=divmod(A*2**m-u,B)
 if r or not(0<d<=2**(m-1)):return None
 out=[];cost=2*b
 for p,w in ((5,5),(19,8),(761,14)):
  e=0
  while d%p==0:d//=p;e+=1;cost+=w
  out.append(e)
 if d!=1 or cost>m:return None
 return (b,*out,m-cost)
def modular(r,s,N):
 b=bval(N,s);u=N//3**b;M=2**r;K=3**(s-b)
 d=(-u*pow(B,-1,M))%M if M>1 else 0
 e,f=next((a,c) for a in range(max(1,2**max(0,r-2)))for c in range(2)if (pow(5,a,M)*pow(19,c,M))%M==d)
 g=0
 if K>1:
  want=(A-u)*pow(B,-1,K)%K
  if (pow(5,e,3)*pow(19,f,3))%3!=want%3:
   g=1;e+=next(2*a for a in range(max(1,2**max(0,r-3)))if 761*pow(5,2*a,M)%M==1%M)
  L=2**max(r-2,1)
  e+=next(L*a for a in range(3**(s-b-1))if pow(5,e+L*a,K)*pow(19,f,K)*pow(761,g,K)%K==want)
 cost=2*b+5*e+8*f+14*g
 phi=2*3**(s-b-1) if K>1 else 1
 m=phi*((max(r,cost,1)+phi-1)//phi)
 return b,e,f,g,m-cost
def count(z):
 b,e,f,g,p=z;m=2*b+5*e+8*f+14*g+p
 return 3**b*(A*2**m-B*5**e*19**f*761**g)
for r in range(5):
 for s in range(5):
  modulus=2**r*3**s
  for N in range(1,2*modulus,2):
   z=modular(r,s,N);n=count(z)
   if n%modulus!=N%modulus or compile_target(n)!=z:raise ValueError((r,s,N,z))
print('all odd residues at 25 dyadic/ternary modulus settings and exact recovery of every recipe')
from math import isqrt
for N,z,p in ((9411754463,(0,0,1,0,3),9411754463),
              (4705885849,(0,1,0,0,5),672269407)):
 if count(z)!=N or compile_target(N)!=z or any(p%d==0 for d in range(2,isqrt(p)+1)):
  raise ValueError((N,z,p))
print('prime and seven-times-prime absolute examples verified')
