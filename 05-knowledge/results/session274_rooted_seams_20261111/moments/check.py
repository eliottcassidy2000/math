import itertools,json,math
from collections import defaultdict
from fractions import Fraction

def hp(n,mask,edges):
 a=[0]*n
 for i in range(n-1): a[i]|=1<<(i+1)
 for k,(i,j) in enumerate(edges):
  if mask>>k&1:a[j]|=1<<i
  else:a[i]|=1<<j
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i,v in enumerate(dp[s]):
   if not v:continue
   z=a[i]&~s
   while z:
    b=z&-z;j=b.bit_length()-1;z-=b;dp[s|b][j]+=v
 return sum(dp[-1])
def mean(n):
 ans=0
 for s in range(1<<(n-1)):
  runs=[];r=0
  for j in range(n):
   if s>>j&1:r+=1
   elif r:runs.append(r);r=0
  if all(x%2==0 for x in runs):ans+=2**len(runs)*math.factorial(n-s.bit_count())
 return Fraction(ans,2**(n-1))
def mean_transfer(n):
 A=[1];O=[0];E=[0]
 for _ in range(n-1):
  add=lambda a,b:[(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
  A,O,E=add(A,E),[0]+add([2*x for x in A],E),[0]+O
 F=add(A,E) if n>1 else A
 return Fraction(sum(c*math.factorial(n-k) for k,c in enumerate(F)),2**(n-1))
for n in range(1,16): assert mean(n)==mean_transfer(n)
out={}
for n in range(2,7):
 edges=[(i,j) for i in range(n) for j in range(i+2,n)]
 cuts=set(sum((((s>>i)^(s>>j))&1)<<k for k,(i,j) in enumerate(edges)) for s in range(1<<n))
 unseen=set(range(1<<len(edges)));rows=[]
 while unseen:
  m=min(unseen);orb={m^c for c in cuts};unseen-=orb
  vals=[hp(n,x,edges) for x in sorted(orb)]
  rows.append({'representative':m,'size':len(orb),'sum':sum(vals),'sumsq':sum(x*x for x in vals),'spectrum':sorted(set(vals))})
 assert all(Fraction(r['sum'],r['size'])==mean(n) for r in rows)
 out[n]={'mean':str(mean(n)),'orbit_count':len(rows),'second_moments':sorted(set(str(Fraction(r['sumsq'],r['size'])) for r in rows)),'orbits':rows}
print(json.dumps(out,indent=2))
