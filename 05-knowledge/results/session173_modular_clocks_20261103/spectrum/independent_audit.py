#!/usr/bin/env python3
"""Independent permutation/descent profiles and reverse subset recursion."""
import itertools,math,json
from functools import lru_cache
B=[(16,1,3,7,14),(20,1,2,7,14)]
def check(ok,label):
 if not ok:raise ValueError(label)
def profile(a):
 n=len(a);hist=[0]*n
 for p in itertools.permutations(range(n)):
  d=sum(not (a[x]>>y&1) for x,y in zip(p,p[1:]))
  hist[d]+=1
 out=[0]*(n+1)
 for k in range(1,n+1):
  q=sum(v*math.comb(n-1-d,k-1-d) for d,v in enumerate(hist) if k>d)
  check(q%math.factorial(k)==0,'unordered division')
  out[k]=q//math.factorial(k)
 return out
def extension(a):
 n=len(a);return tuple(a[i]|(1<<n if i!=0 else 0) for i in range(n))+(1,)
def substitute(D,parts):
 addresses=[(i,j) for i,p in enumerate(parts) for j in range(len(p))]
 return tuple(sum(1<<y for y,(k,l) in enumerate(addresses) if (parts[i][j]>>l&1 if i==k else D[i]>>k&1)) for i,j in addresses)
def wrapped(bits,r,K,address=None):
 parts=[extension(B[b]) if i==address else B[b] for i,b in enumerate(bits)]
 # Individual seed positions followed by a transitive feedback module.
 D=[];m=len(parts)
 for i in range(m):D.append(sum(1<<j for j in range(i+1,m)) | (1<<m if i>=K else 0))
 D.append((1<<K)-1)
 Z=tuple(sum(1<<j for j in range(i+1,r)) for i in range(r))
 return substitute(D,parts+[Z])
def hamilton(a):
 n=len(a);incoming=tuple(sum(1<<i for i in range(n) if a[i]>>j&1) for j in range(n))
 @lru_cache(None)
 def rec(mask,last):
  if mask==1<<last:return 1
  rest=mask^(1<<last);cand=incoming[last]&rest;total=0
  while cand:
   bit=cand&-cand;cand-=bit;total+=rec(rest,bit.bit_length()-1)
  return total
 ans=sum(rec((1<<n)-1,j) for j in range(n));rec.cache_clear();return ans
def restrict(a,keep):
 return tuple(sum(1<<j for j,y in enumerate(keep) if a[x]>>y&1) for x in keep)
def is_strong(a):
 n=len(a)
 for start in range(n):
  reached={start};front=[start]
  while front:
   x=front.pop()
   for y in range(n):
    if a[x]>>y&1 and y not in reached:reached.add(y);front.append(y)
  if len(reached)!=n:return False
 return True
for endpoint in (0,2):
 keep=[i for i in range(5) if i!=endpoint]
 check(restrict(B[0],keep)==restrict(B[1],keep),'bit erasure hostile')
# Equality of side size and r permits a competing minimum-sized cut.
hostile=wrapped((0,0),5,1)
check(not is_strong(restrict(hostile,list(range(5,15)))),'strict side-size boundary')
check(not is_strong(restrict(hostile,list(range(10)))),'feedback cut boundary')
p=[profile(a) for a in B];q=[profile(extension(a)) for a in B]
check(p[0]==p[1],'seed profile equality')
check(q==[[0,17,95,138,73,15,1],[0,19,109,150,75,15,1]],'query profile')
contexts=[]
for mask in range(8):
 D=[0]*3
 for k,(i,j) in enumerate(itertools.combinations(range(3),2)):
  if mask>>k&1:D[i]|=1<<j
  else:D[j]|=1<<i
 profiles=[profile(substitute(D,[seed,(0,),(0,)])) for seed in B]
 check(profiles[0]==profiles[1],'all three-position contexts')
 contexts.append({'context':D,'profile':profiles[0]})
# Independent literal read gaps at r beyond the main two-feedback formula.
reads=[]
for r in (1,2,3,4):
 for address in (0,1):
  results=[]
  for bit in (0,1):
   bits=[0,0];bits[address]=bit
   results.append(hamilton(wrapped(bits,r,1,address)))
  gap=results[1]-results[0]
  check(gap>=18,'positive arbitrary-r lower bound')
  if r==2:check(results==[585033,645915],'r2 read formula')
  reads.append({'r':r,'address':address,'responses':results,'gap':gap})
# Two separate closed-form derivations: evaluate q basis and transformed f basis.
C=((5,3,0),(3,7,2),(0,2,2));E=((15,-11,2),(-11,13,-4),(2,-4,2))
T=((1,0,0),(-2,1,0),(3,-3,1))
check(tuple(tuple(sum(T[k][i]*C[k][l]*T[l][j] for k in range(3) for l in range(3)) for j in range(3)) for i in range(3))==E,'basis transformation')
def evalp(p,u):
 return sum(v*math.prod(u-j for j in range(k)) for k,v in enumerate(p))
check(all(evalp(p[0],u)==u**5+6*u**3+2*u for u in range(8)),'seed power polynomial')
s=[evalp(p[0],u) for u in (1,2,3)];f=[[evalp(z,u) for u in (1,2,3)] for z in q]
for M in range(2,31):
 calc=sum(s[i]*E[i][j]*s[j]**(M-1) for i in range(3) for j in range(3))
 closed=33*9**(M-1)-651*84**(M-1)+504*411**(M-1)
 check(calc==closed,'storage all tested M')
 for b in (0,1):
  read0=sum(f[b][i]*E[i][j]*s[j]**(M-1) for i in range(3) for j in range(3))
  read1=sum(s[i]*E[i][j]*f[b][j]*s[j]**(M-2) for i in range(3) for j in range(3))
  check(read0==(689+2*b)*9**(M-1)-(3071+254*b)*84**(M-1)+(2036+200*b)*411**(M-1),'read0 formula')
  check(read1==(561+66*b)*9**(M-2)-(145824+20832*b)*84**(M-2)+(730296+81648*b)*411**(M-2),'read1 formula')
print(json.dumps({'status':'PASS','method':'permutation descent profiles; independent incoming subset recursion; explicit basis transform','seed_profiles':p,'query_profiles':q,'substitution_contexts':contexts,'literal_reads':reads,'closed_formula_M':[2,30],'hostiles':['either endpoint deletion erases a bit','side size equal to r permits competing r-cut']},indent=2))
