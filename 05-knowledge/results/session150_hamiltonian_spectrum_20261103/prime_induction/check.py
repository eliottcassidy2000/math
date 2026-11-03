import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'gadgets'))
from test import graph,profile,add
def require(condition, details):
 if not condition: raise AssertionError(details)
nchecks=0
for n in range(2,6):
 for code in range(1<<(n*(n-1)//2)):
  a=graph(n,code)
  for v in range(n):
   h,s,q,f=profile(a,v); b=a
   for k in range(1,4):
    b=add(b,v)
    hp,sp,qp,fp=profile(b,v)
    require((hp,sp,qp,fp)==(h-s+q+(2**k-1)*(f-q),q,q,q+2**k*(f-q)),(n,code,v,k))
    nchecks+=1
print('Rooted profile update verified:',nchecks,'checks')

for r in range(1,101):
 m=12*r+1
 require(m%3 and (1<<m)-1>(2*r*m)**2, (r,m))
q=31
res={(5+9*pow(2,k,q))%q for k in range(5)}
require(res=={10,14,15,23,25}, res)
require({(9*x)%q for x in res}=={2,8,11,21,28}, 'inverse residues')
require(1 not in res and 1 not in {(9*x)%q for x in res}, 'excluded residue')
p=311
require(p%31==1 and all(p%d for d in range(2,int(p**.5)+1)), 'example primality')
print('Modular example and explicit m inequality: OK')
