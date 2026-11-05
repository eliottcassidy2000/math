from math import gcd
from fractions import Fraction
import json

def need(q,m):
 if not q:raise RuntimeError(m)
def step(s,a):
 D,E,U,d=s
 if a=='O':return D,0,2*U,d
 if a=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D
def apply(w,s=(1,0,2,1)):
 for a in w:s=step(s,a)
 return s
def N(v,w):
 D,E,U,d=v;g,s,A,b=w
 return (D-E+U)*g-(D-E+d)*s+(D+3*U-d)*A-(U+E)*b
def coefficient(bits):
 k=len(bits);A=sum(e*2**i for i,e in enumerate(bits));R=sum(e*2**(k-i-1) for i,e in enumerate(bits))
 C=sum(bits[i]*bits[j]*2**(j-i-1) for i in range(k) for j in range(i+1,k))+sum(bits[i-1]*bits[i] for i in range(k))
 return A,R,C
def g(n):
 acc=0
 while n>1:
  if n%2==0:n//=2
  else:n=(n+1)//2;acc+=n-n%2
 return 1+acc
clockchecks=orbitchecks=0
for b in range(3,54,2):
 k=1
 while (2**k-1)%b:k+=1
 B=2**k
 for a in range(b//2+1,b):
  if gcd(a,b)!=1:continue
  bits=[(a*(B-1)//b>>i)&1 for i in range(k)];A,R,C=coefficient(bits);K=a*R-b*C
  need(K>0 and K%2==0,('decrement',a,b,K))
  w=''.join('Z' if e else 'O' for e in bits)
  for D,U,d in [(0,0,1),(1,0,0),(0,1,0),(17,29,3),(-2,5,-1)]:
   v=apply(w,(D,d,U,d));need(v==(D+A*U-C*d,d,B*U-R*d,d),('E boundary',a,b,D,U,d))
  for j in range(k):
   for shift in (bits[j:]+bits[:j],[1-e for e in bits[j:]+bits[:j]]):
    AA,RR,CC=coefficient(shift);aa=AA*b//(B-1)
    need(AA*b%(B-1)==0 and aa*RR-b*CC==K,('necklace',a,b,j))
    orbitchecks+=1
  clockchecks+=1
basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
for v in basis:
 for w in basis:
  D,E,U,d=v;A,B,V,b=w
  rho=D-E+U-d;p=D+d;sig=A-B+V-b;q=A+b
  need(N(v,step(w,'J'))-N(step(v,'J'),w)==rho*q-p*sig,'defect')
v=apply('ZJOZO');w=apply('ZOZZOZZOOOZ'+'ZZOOOZ'*2)
need(v==(17,0,42,3) and w==(6003419,1,10806159,1),'short companion')
need(g(w[2])==w[0] and 5*w[2]-9*w[0]==24,'exact independent g')
need(N(v,step(w,'J'))==N(step(v,'J'),w)==2972894005,'new exact common target')
D,E,U,d=v;Av=D-E+3*U-2*d;Bv=2*D-E+4*U-d
need((Av,Bv)==(137,199),'linear coefficients')
h=(w[2]-3)//18
need(h==600342 and 2*(Av*5+Bv*9)*h-Av+3*Bv+d-U==2972894005,'fixed target')
need(apply('ZOOZ')==(15,1,23,1) and apply('ZOZOZ')==(31,1,43,1),'small collision sources')
need(N(apply('ZOOZ'),step(apply('ZOZOZ'),'J'))==7649,'small collision target')
print(json.dumps({'status':'PASS','clocks_checked_odd_denominator_through':53,'clock_count':clockchecks,'cyclic_complement_checks':orbitchecks,'short_companion':w,'new_C2_target':2972894005,'h':h,'small_C2_target':7649},sort_keys=True,indent=2))
