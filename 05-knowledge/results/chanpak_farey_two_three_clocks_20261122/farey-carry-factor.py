from math import isqrt,gcd
from itertools import product
from functools import lru_cache
from math import prod
from pathlib import Path
import json
import argparse
N=2**121-2**60-1
rows=[(4,2465,2293760,1344270,712944,1740,986,423385),(8,37265,144834560,84839910,57602544,24360,14906,33838015)]
def f(x):
 s=0
 while x:
  y=x>>1
  if not x&1:s+=y+(y&1)
  x=y
 return s
def carries(row,R0,ell,side=0):
 n,h,q,c,t,e,r,k=row; B=1<<ell
 high=isqrt((N-k)*(R0+1)**2//(q*R0)) if R0 else None
 words=range(1<<side)
 for a,b,u,v in product((0,1),(0,1),words,words):
  i=(u<<(ell-side)) if a==0 else B-1-(u<<(ell-side))
  j=(v<<(ell-side)) if b==0 else B-1-(v<<(ell-side))
  def rev(d):
   return int(f"{d:0{ell}b}"[::-1],2)
  ei=B-1-rev(i)-B*(1-a);ej=B-1-rev(j)-B*(1-b)
  local_high=high
  if local_high is None:
   mX=i or B;mY=j or B
   d=(N-k)//q
   if mX*mY>d:continue
   local_high=max(mX+d//mX,mY+d//mY)
  Z=(local_high-i-j)//B
  if Z<0:continue
  P=q*B; A=t+q*j+h*(1-a);C=c+q*i+h*(1-b)
  K=k+q*i*j+t*i+c*j+r*(i&1)+e*(j&1)+h*(f(i)+f(j))
  low=Z*min(0,ei,ej);upper=Z*max(0,ei+1,ej+1)
  if a==0 and b==0: assert upper==0
  R=(N-K)*pow(h,-1,B)%B
  for E in range(low+(R-low)%B,upper+1,B):
   V=(N-K-h*E)//B
   D=P*V+A*C
   if D<=0:continue
   yield i,j,A,C,P,D,E
FACTORS={}
WITNESSES={}
CERTIFIED_ONLY=False

def factor(n):
 if n not in FACTORS:
  assert not CERTIFIED_ONLY, ("missing factorization",n)
  from sympy import factorint
  FACTORS[n]={int(p):int(e) for p,e in factorint(n).items()}
 fac=FACTORS[n]
 assert all(p>=2 and e>=1 for p,e in fac.items())
 assert prod(p**e for p,e in fac.items())==n
 return fac

@lru_cache(None)
def certify(p):
 if p==2:return True
 if p%2==0:return False
 fac=factor(p-1)
 if not all(certify(q) for q in fac):return False
 ok=lambda a:pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fac)
 if p not in WITNESSES:
  assert not CERTIFIED_ONLY, ("missing witness",p)
  WITNESSES[p]=next(a for a in range(2,min(p,1000)) if ok(a))
 return ok(WITNESSES[p])

def solve(row,R0,ell,side):
 n,h,q,c,t,e,r,k=row;B=1<<ell;Nc=0;count=0
 for i,j,A,C,P,D,E in carries(row,R0,ell,side):
  Nc+=1
  fac=factor(D)
  print(n,ell,R0,i,j,E,D,dict(fac),flush=True)
  assert prod(int(q)**e for q,e in fac.items())==D
  assert all(certify(int(q)) for q in fac)
  pool=[1]
  for q,e in fac.items():pool=[a*int(q)**j for a in pool for j in range(e+1)]
  for d in pool:
   X,a=divmod(d-C,P);Y,b=divmod(D//d-A,P)
   if a or b or min(X,Y)<0:continue
   count+=1
   x=B*X+i;y=B*Y+j
   if not min(x,y) or R0 and max(x,y)>R0*min(x,y):continue
   if k+q*x*y+t*x+c*y+r*(x&1)+e*(y&1)+h*(f(x)+f(y))==N:
    print("HIT",n,x,y,flush=True)
 print("TOTAL", n,ell,R0,side,Nc,count,flush=True)
if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("ell",type=int);a.add_argument("ratio",type=int);a.add_argument("side",type=int,nargs="?",default=0);a.add_argument("--spine",type=int)
 a.add_argument("--certificate",type=Path)
 a.add_argument("--record",type=Path)
 a=a.parse_args()
 if a.certificate:
  CERTIFIED_ONLY=True
  j=json.loads(a.certificate.read_text())
  FACTORS.update({int(n):{int(p):e for p,e in fac.items()} for n,fac in j["factors"].items()})
  WITNESSES.update({int(p):a for p,a in j["witnesses"].items()})
 for row in rows:
  if a.spine is None or a.spine==row[0]:solve(row,a.ratio,a.ell,a.side)

 if a.record:
  a.record.write_text(json.dumps({"factors":FACTORS,"witnesses":WITNESSES},indent=1)+"\n")
