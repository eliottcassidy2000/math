"""Standalone exact controls for signed one-step closure of valuation prefixes."""
from itertools import product
from collections import Counter
import json

def need(c,msg):
 if not c: raise ArithmeticError(msg)
def stats(w):
 A=B=0
 for a in w:B,A=3*B+2**A,A+a
 return A,B
def step(n,s):
 z=3*n+s;need(z>0,'positive raw')
 a=(z&-z).bit_length()-1
 return z//2**a,a
def construct(w,s,y):
 A,B=stats(w);r=len(w);m=3**(r+1);L=2*3**r
 rhs=s*(2**A+3*B)%m
 hits=[b for b in range(1,L+1) if pow(2,A+b,m)*y%m==rhs]
 need(len(hits)==1,'unique exponent phase')
 b=hits[0]
 while 2**(A+b)*y-s*(2**A+3*B)<=0:b+=L
 num=2**(A+b)*y-s*(2**A+3*B)
 need(num%m==0,'integral closure')
 return num//m,b,L
counts=Counter();examples=[]
for r in range(5):
 for w in product(range(1,4),repeat=r):
  for s in (-1,1):
   for y in (1,5,7,17,25):
    n,b,L=construct(w,s,y)
    for extra in (0,1,2):
     A,B=stats(w);bb=b+extra*L
     nn=(2**(A+bb)*y-s*(2**A+3*B))//3**(r+1)
     need(nn>0 and nn%2==1,'positive odd')
     z=nn;got=[]
     for _ in range(r+1):z,a=step(z,s);got.append(a)
     need(tuple(got)==w+(bb,) and z==y,'exact replay')
     counts['closure_replays']+=1
    counts['unique_phase_checks']+=1
    if w==(1,2) and s==-1 and y in (1,5,17):
     examples.append({'prefix':w,'sign':s,'target':y,'start':n,'last_valuation':b,'exponent_period':L})
# Why the target restriction is essential: no odd-map image is divisible by 3.
for n in range(1,2000,2):
 for s in (-1,1):need(step(n,s)[0]%3!=0,'image excludes multiples of three');counts['target_obstruction_checks']+=1
print(json.dumps({'counts':dict(counts),'three_basin_examples':examples},indent=2))
