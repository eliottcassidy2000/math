#!/usr/bin/env python3
"""Independent exact shell decoder; no producer imports; energy baseline cancelled."""
import random,json

def must(v,msg):
 if not v:raise RuntimeError(msg)
def coeff(k):
 n=k-2;a=[]
 for i in range(n):
  x,y=k-i-3,i
  num=-34104960*471**x*104**y+3945600*104**x*471**y-140160*104**(x+y)
  must(num%367==0,'integer S')
  a.append(2*13**(k-2)+2*num//367)
 b={(i,j):(-1658880*104**(k-4)-18*13**(k-2) if j-i==1 else -55296000*104**(k-j+i-3)*471**(j-i-2))for i in range(n)for j in range(i+1,n)}
 return a,b
def energy(s,a,b):return sum(x*y for x,y in zip(s,a))+sum(v*s[i]*s[j]for(i,j),v in b.items())
def decode(z,k,a,b):
 n=k-2
 if k==3:
  return [z//a[0]] if z in (0,a[0])else None
 known={}
 for t in range((n+1)//2):
  rest=z-sum(a[i]*v for i,v in known.items())-sum(v*known[i]*known[j] for(i,j),v in b.items()if i in known and j in known)
  if rest%13**t:return None
  lam=pow(8,t,13)*pow(471,k-t-5,13)%13
  digit=(rest//13**t)*pow(lam,-1,13)%13
  i,j=t,n-1-t
  if i==j:
   options=[(x,)for x in (0,1)if ((5-6*(known[0]+known[n-1]))*x)%13==digit]
  elif t==0:options=[(x,y)for x in(0,1)for y in(0,1)if (3*x+2*y-6*x*y)%13==digit]
  else:options=[(x,y)for x in(0,1)for y in(0,1)if ((3-6*known[n-1])*x+(2-6*known[0])*y)%13==digit]
  if not options:return None
  must(len(options)==1,'unique digit')
  known[i]=options[0][0];known[j]=options[0][-1]
 s=[known[i]for i in range(n)]
 return s if energy(s,a,b)==z else None

def run():
 rng=random.Random(286);total=0;off=0
 for k in range(3,101):
  a,b=coeff(k);n=k-2
  vectors=([[mask>>i&1 for i in range(n)]for mask in range(1<<n)]if k<=12 else [[rng.randrange(2)for _ in range(n)]for _ in range(12)])
  seen=set()
  for s in vectors:
   z=energy(s,a,b);must(decode(z,k,a,b)==s,'inverse');seen.add(z);total+=1
   # Any accepted offset must reconstruct exactly, not merely pass low digits.
   d=decode(z+2,k,a,b)
   must(d is None or energy(d,a,b)==z+2,'offset final equality');off+=1
  if k<=12:must(len(seen)==1<<n,'exhaustive injectivity')
 return dict(status='PASS',minimum_k=3,maximum_k=100,exact_inverse_controls=total,offset_controls=off,scope='finite controls of separately proved all-order valuation argument')
if __name__=='__main__':print(json.dumps(run(),indent=2))
