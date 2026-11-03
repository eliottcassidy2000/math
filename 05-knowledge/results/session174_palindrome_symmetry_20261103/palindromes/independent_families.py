from itertools import permutations
from math import isqrt
import json

def check(v,m):
 if not v:raise RuntimeError(m)
def factor(n):
 fs=[];d=2
 while d*d<=n:
  while n%d==0:fs.append(d);n//=d
  d+=1
 if n>1:fs.append(n)
 return fs
def sign_cycles(p):
 seen=set();c=0
 for i in range(len(p)):
  if i in seen:continue
  c+=1;j=i
  while j not in seen:seen.add(j);j=p[j]
 return (-1)**(len(p)-c)
factors={}
for m in range(1,5):
 n=10**(2*m+1)-10**m-1;fs=factor(n);factors[str(m)]=fs
 for q in fs:check(q!=41 and pow(q,20,41)==1,'decimal prime divisor subgroup')
for n in range(3):
 x=10**(n+1);B=100*x*x+666*x+1
 for q in factor(B):
  if q not in [2,5,7,17,19]:check(pow(2261,(q-1)//2,q)==1,'Belphegor factor discriminant')
# Orbit enumeration independent of exponentiation formula.
orbit_cases=0
for b in range(2,18):
 for q in range(3,501,2):
  if any(q%d==0 for d in range(2,isqrt(q)+1)) or b%q==0:continue
  orbit=[];x=1
  while x not in orbit:orbit.append(x);x=x*b%q
  check(x==1,'cyclic orbit return');t=len(orbit)
  roots={j for j,x in enumerate(orbit) if (b*x*x-x-1)%q==0}
  if t%2:check(len(roots)<=1,'odd phase')
  else:check(roots=={(t//2-1-r)%t for r in roots},'reflection pairing')
  orbit_cases+=1
# Additional n=6 examples, with permutation sign computed from cycle counts.
perms=list(permutations(range(6)));edges=[(i,j) for i in range(6) for j in range(i+1,6)]
for mask in [0,1,3,15,255,1023,12345,32767]:
 arcs={(i,j) if (mask>>k)&1 else (j,i) for k,(i,j) in enumerate(edges)}
 f=[0]*6;s=[0]*6
 for p in perms:
  k=sum((a,b) in arcs for a,b in zip(p,p[1:]));f[k]+=1;s[k]+=sign_cycles(p)
 check(f==f[::-1] and s==[-v for v in s[::-1]],'n6 reversal')
 check(sum(s)==0 and sum((-1)**j*v for j,v in enumerate(f))==0,'n6 factors')
print(json.dumps({'status':'PASS','factor_examples':factors,'independent_orbits':orbit_cases,'additional_n6_tournaments':8},indent=2))
