exec(open(__file__.replace('explore_modules.py','explore_gluing.py')).read().split('A=(20,')[0])
from itertools import combinations

def tour(n,bits):
 out=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if bits>>k&1:out[i]|=1<<j
  else:out[j]|=1<<i
 return tuple(out)

def covers(T):
 n=len(T);e=m.endpoints(T);h=[sum(x) for x in e];f=[0]*(1<<n);f[0]=1;v=[]
 for r in range(1,n+1):
  g=[0]*(1<<n)
  for mask in range(1,1<<n):
   A=mask
   while A:
    g[mask]+=h[A]*f[mask^A];A=(A-1)&mask
  v.append(g[-1]);f=g
 return tuple(v)

def rootprofile(T,v):
 return tuple(sum(m.endpoints(m.inflate(T,0,v,k)[0])[-1]) for k in range(0,5))

def canon(T):
 from itertools import permutations
 n=len(T);best=None
 for p in permutations(range(n)):
  z=tuple(sum(1<<j for j in range(n) if T[p[i]]>>p[j]&1) for i in range(n))
  if best is None or z<best:best=z
 return best

D={};seen=set()
for bits in range(1024):
 T=tour(5,bits);C=canon(T)
 if C in seen:continue
 seen.add(C);q=covers(C);r=sorted(rootprofile(C,v) for v in range(5));D.setdefault(q,[]).append((C,r))
for q,ls in D.items():
 if len({str(x[1]) for x in ls})>1:
  print('COLLISION',q,ls)
print('classes',len(seen),'profiles',len(D))
from math import gcd
for q,ls in D.items():
 print('PROFILE',q,'gcd',gcd(*q),'classes',len(ls))
