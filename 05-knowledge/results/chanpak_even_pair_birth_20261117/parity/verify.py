from itertools import permutations,product
from fractions import Fraction
from math import factorial
from functools import cache
import json

def require(b):
 if not b:raise ArithmeticError('parity control failed')

def path_table(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    rest=s^(1<<v)
    if rest:dp[s][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and rows[u]>>v&1)
 return [1]+[sum(r) for r in dp[1:]]
def anti(rows):
 n=len(rows);degrees=[r.bit_count() for r in rows]
 for p in permutations(range(n)):
  if any(p[p[i]]!=i or degrees[i]+degrees[p[i]]!=n-1 for i in range(n)):continue
  if all(((rows[i]>>j)&1)==((rows[p[j]]>>p[i])&1) for i in range(n) for j in range(n)):
   return p
 return None

def polynomial(rows):
 n=len(rows);H=path_table(rows)
 @cache
 def C(mask,r):
  if r==0:return int(mask==0)
  s=mask;ans=0
  while s:
   ans+=H[s]*C(mask^s,r-1);s=(s-1)&mask
  return ans
 result=[Fraction(0)]*(n+1);fall=[Fraction(1)]
 for r in range(1,n+1):
  nxt=[Fraction(0)]*(len(fall)+1)
  for i,x in enumerate(fall):nxt[i]-=(r-1)*x;nxt[i+1]+=x
  fall=nxt
  for i,x in enumerate(fall):result[i]+=C((1<<n)-1,r)*x/factorial(r)
 require(all(x.denominator==1 for x in result))
 return [int(x) for x in result]

def chain(parts):
 n=1+sum(parts);rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 l=0
 for p in parts:
  rows[l]^=1<<(l+p);rows[l+p]|=1<<l;l+=p
 return rows

def strong(rows):
 n=len(rows)
 for i in range(n):
  seen=1<<i
  while True:
   nxt=seen
   for j in range(n):
    if seen>>j&1:nxt|=rows[j]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True

symmetry=0;max5=0;max4=0
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
  rows=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if code>>k&1:rows[i]|=1<<j
    else:rows[j]|=1<<i
    k+=1
  if n==5:max5=max(max5,path_table(rows)[-1])
  if n==4:max4=max(max4,path_table(rows)[-1])
  p=anti(rows)
  if p is None:continue
  coeff=polynomial(rows);expected=[0]*(n+1);expected[n]=1
  if n%2:
   v=next(i for i in range(n) if p[i]==i)
   for length in range(3,n+1,2):
    for tail in permutations([i for i in range(n) if i!=v],length-1):
     cyc=(v,)+tail
     if all(rows[cyc[i]]>>cyc[(i+1)%length]&1 for i in range(length)) and all(p[cyc[i]]==cyc[-i%length] for i in range(length)):
      expected[n-length+1]+=2
  require(all((x-y)%4==0 for x,y in zip(coeff,expected)));symmetry+=1
chains=0
for a in range(4):
 for b in range(3):
  if a+b==0 or 2*a+3*b+1>8:continue
  total=fix=0;predicted=predfix=0
  for c in range(a//2+1):
   A=a-2*c;B=b;C=c;L=A+B+C
   predicted+=factorial(L)//(factorial(A)*factorial(B)*factorial(C))
   if A%2+B%2+C%2<=1:predfix+=factorial(L//2)//(factorial(A//2)*factorial(B//2)*factorial(C//2))
   for parts in set(permutations([2]*A+[3]*B+[4]*C)):
    rows=chain(parts);require(strong(rows) and path_table(rows)[-1]==3**a*5**b)
    sc=anti(rows) is not None;require(sc==(parts==parts[::-1]))
    total+=1;fix+=sc;chains+=1
  require((total,fix)==(predicted,predfix))
regular5=[sum(1<<((i+j)%5) for j in (1,2)) for i in range(5)]
require(max4==5 and path_table(regular5)[-1]==15 and path_table(chain((2,3)))[-1]==15)
rows33=[2,12,25,49,35,7]
require(max5==15 and path_table(rows33)[-1]==33 and strong(rows33))
require(path_table(chain((6,)))[-1]==33 and strong(chain((6,))))
print(json.dumps({'self_converse_labelled_polynomial_controls':symmetry,'literal_smooth_chains':chains,'order5_max_H':max5,'H15_odd_birth_even_chain':[5,6],'H33_even_birth_strong_odd_recovery':[6,7]},indent=2))
