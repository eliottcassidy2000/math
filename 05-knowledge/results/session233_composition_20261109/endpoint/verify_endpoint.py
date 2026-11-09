import sys

def require(condition, message):
 if not condition:
  raise RuntimeError(message)

if "--fail" in sys.argv:
 require(False, "Deliberate verification failure sentinel")

from itertools import permutations
import json
from probe import ends

def base(n):
 a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 a[0]^=1<<(n-1); a[-1]|=1
 return a

def extend(n,N):
 a=base(n)+[0]
 for i in range(n):
  if 1<=i<=n-2 and N>>(i-1)&1: a[i]|=1<<n
  else: a[n]|=1<<i
 return a

def brute(a):
 e=[0]*len(a)
 for p in permutations(range(len(a))):
  if all(a[u]>>v&1 for u,v in zip(p,p[1:])): e[p[-1]]+=1
 return e

def reach(a):
 n=len(a)
 for i in range(n):
  s=1<<i
  while True:
   t=s
   for j in range(n):
    if s>>j&1:t|=a[j]
   if t==s:break
   s=t
  if s!=(1<<n)-1:return False
 return True

def aut(a,fix=None):
 n=len(a);d=[v.bit_count() for v in a]; out=[]
 for p in permutations(range(n)):
  if fix is not None and p[fix]!=fix:continue
  if any(d[i]!=d[p[i]] for i in range(n)):continue
  if all((a[i]>>j&1)==(a[p[i]]>>p[j]&1) for i in range(n) for j in range(n)):out.append(p)
 return len(out)
summary=[]
for n in range(4,9):
 a=base(n); expected=[1]+[1<<(i-1) for i in range(1,n-1)]+[1]
 require(ends(a)==expected, f"Base terminal vector failed: n={n}")
 if n<=7:require(brute(a)==expected and aut(a)==1, f"Independent base check failed: n={n}")
 for N in range(1,1<<(n-2)):
  b=extend(n,N);e=ends(b)
  require(e[-1]==N and reach(b), f"Endpoint/strong check failed: n={n}, N={N}")
  if n<=6:require(brute(b)==e and aut(b,n)==1, f"Independent selector check failed: n={n}, N={N}")
 summary.append({'base_n':n,'constructed':(1<<(n-2))-1,'endpoint_range':[1,(1<<(n-2))-1]})
for N in [7,21,37,101,257]:
 n=max(4,N.bit_length()+2);a=extend(n,N);e=ends(a)
 summary.append({'N':N,'order':len(a),'H':sum(e),'endpoint':e[-1]})
print(json.dumps(summary,indent=2))
