#!/usr/bin/env python3
"""Finite checks of the proved iterated-C3 transport family. Standard library."""
import json
from itertools import product
from check_transport import compose,inverse,conj,transport,holonomy,acts,need
import check_transport as ct

def tournament(k):
 N=3**k
 def forward(i,j):
  while i//3!=j//3:i//=3;j//=3
  return (j-i)%3==1
 return tuple(sum(1<<j for j in range(N) if i!=j and forward(i,j)) for i in range(N))
def negate(k):
 out=[]
 for i in range(3**k):
  x=i;y=0;place=1
  for _ in range(k):y+=((-x)%3)*place;x//=3;place*=3
  out.append(y)
 return tuple(out)
def extend_mark(s,k):
 # Last two digits form a tau-stable block whose prefix consists only of zeros.
 out=list(negate(k));out[:9]=s
 return tuple(out)
def run():
 s=(0,2,1,6,8,7,3,5,4);t=(0,2,1,7,6,8,4,3,5);u=(3,5,4,0,2,1,6,8,7)
 baseh=holonomy(u,t,s);reports=[]
 for k in range(5):
  N=3**k;W=tournament(k);tau=negate(k);e=tuple(range(N))
  need(acts(W,tau,True) and compose(tau,tau)==e,'digit negation anti')
  need(all(x.bit_count()==(N-1)//2 for x in W),'regular')
  if k:
   size=N//3
   for i,j in product(range(N),repeat=2):
    d=(W[i]^W[j]).bit_count()
    need((d<size)==(i//size==j//size),'intrinsic blocks')
    if i//size!=j//size:need(d==2*size,'cross-block distance')
  a=(N-1)//2;c=(N-2*k-1)//4;z=(N+2*k-1)//4
  need(a==k+2*c and z==k+c,'mass exponents')
  need((N-2*k-1)%4==0 and (N+2*k-1)%4==0,'integer exponents')
  if k:
   oldN=3**(k-1);olda=(oldN-1)//2;oldc=(oldN-2*(k-1)-1)//4
   need(a==3*olda+1 and c==oldc+olda,'wreath recurrence')
  if k>=2:
   S,T,U=(extend_mark(v,k) for v in (s,t,u))
   for v in (S,T,U):need(acts(W,v,True) and compose(v,v)==e,'embedded anti')
   h=holonomy(U,T,S)
   need(h==baseh+tuple(range(9,N)) and h!=e,'nonflat family')
  reports.append({'k':k,'vertices':N,'aut_order':3**a,'centralizer_order':3**c,'anti_involutions':3**z,'anti_involutions_per_fixed_root':3**c,'nonflat':k>=2})
 return {'status':'PASS','active_checks':ct.CHECKS,'finite_scope':'k=0..4, N=1,3,9,27,81; group order formulas proved, not exhaustive group enumeration beyond9','family':reports}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
