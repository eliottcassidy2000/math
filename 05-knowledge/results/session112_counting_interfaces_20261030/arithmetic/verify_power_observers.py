#!/usr/bin/env python3
import math,json

def need(b,msg):
 if not b:raise RuntimeError(msg)
def S(q,j):return (q**j-1)//(q-1)
def val(x,p,k):
 x%=p**k
 if x==0:return k
 e=0
 while x%p==0:e+=1;x//=p
 return e
def decode(p,w,t,k,A,B):
 m=p**k;q=p**w;h=min(val(B-A,p,k),k-1)//w
 return (S(q,h+1)*pow(S(q,t),-1,m)*(B-q**t*A)-A)%m
def aliquot(n,w):
 return sum(d**w for d in range(1,n) if n%d==0)
def kernel(p,w,k):
 if p!=2:return math.gcd(w,(p-1)*p**(k-1))
 if k==1:return 1
 if k==2:return math.gcd(w,2)
 return math.gcd(w,2)*math.gcd(w,2**(k-2))
actual=0;rootchecks=0;formal=0;bezout=0;witness=0
for p in [2,3,5]:
 for w in [1,2,3,4]:
  for t in [1,2]:
   for n in range(1,81):
    A=aliquot(n,w);B=aliquot(n*p**t,w)
    for k in [1,2,3,4]:
     need(decode(p,w,t,k,A,B)==pow(n,w,p**k),'weighted inverse');actual+=1
for p,k in [(2,1),(2,2),(2,3),(2,5),(2,7),(3,1),(3,2),(3,3),(5,1),(5,2),(7,2)]:
 m=p**k
 for w in range(1,9):
  fibres={}
  for n in range(m):fibres.setdefault(pow(n,w,m),[]).append(n)
  for y,xs in fibres.items():
   if y==0:expected=p**(k-(k+w-1)//w)
   else:
    v=val(y,p,k);need(v%w==0,'power valuation')
    e=v//w;expected=p**((w-1)*e)*kernel(p,w,k-w*e)
   need(len(xs)==expected,'exact root fibre');rootchecks+=1
  if m<=32:
   for t in [1,2]:
    for A in range(m):
     need(len({decode(p,w,t,k,A,B) for B in range(m)})==m,'formal weighted permutation')
     for x in range(m):
      q=p**w;h=min(val(x,p,k),k-1)//w
      C=(A+x)*pow(1-q**(h+1),-1,m)%m
      B=(A+S(q,t)*pow(S(q,h+1),-1,m)*(x+q**(h+1)*A))%m
      xp=q**t*x%m;hp=min(val(xp,p,k),k-1)//w
      need((B+xp)*pow(1-q**(hp+1),-1,m)%m==C,'weighted conjugacy');formal+=1
 for n in range(m):
  if n%p:
   need(pow(n,3,m)*pow(pow(n,2,m),-1,m)%m==n,'square cube unit recovery');bezout+=1
for p,k,ns,expected in [(5,2,[2,23],(1,5)),(2,5,[97,47,17,31],(1,6))]:
 m=p**k
 pairs=[(aliquot(n,2)%m,aliquot(p*n,2)%m) for n in ns]
 need(set(pairs)=={expected},'genuine root ambiguity')
 for j in range(8):
  values=[(S(p*p,j+1)+S(p*p,j)*n*n)%m if j else 1 for n in ns]
  need(len(set(values))==1,'all future prime-source observations');witness+=1
need(pow(5,2,25)==pow(10,2,25)==pow(5,3,25)==pow(10,3,25)==0,'nonunit control')
print(json.dumps({'status':'PASS','direct_weighted_divisor_checks':actual,'exact_power_fibres':rootchecks,'formal_conjugacy_checks':formal,'square_cube_unit_recoveries':bezout,'all_future_witness_positions':witness},indent=2))
