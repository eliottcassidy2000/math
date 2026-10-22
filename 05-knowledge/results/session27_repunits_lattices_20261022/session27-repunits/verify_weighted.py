"""Independent labeled enumeration and SCC generating function, n<=6."""
from collections import Counter
from fractions import Fraction
from math import factorial,comb
import json

def H(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if dp[mask][v]:
    avail=out[v]&~mask
    while avail:
     b=avail&-avail;avail-=b;dp[mask|b][b.bit_length()-1]+=dp[mask][v]
 return sum(dp[-1])
def strong(out):
 n=len(out)
 for start in range(n):
  seen=1<<start;old=0
  while old!=seen:
   old=seen
   for v in range(n):
    if seen>>v&1:seen|=out[v]
  if seen!=(1<<n)-1:return False
 return n>=3
allc={}; sc={}
for n in range(1,7):
 a=Counter();s=Counter();edges=[(i,j)for i in range(n)for j in range(i+1,n)]
 for mask in range(1<<len(edges)):
  out=[0]*n
  for b,(i,j) in enumerate(edges):
   if mask>>b&1:out[i]|=1<<j
   else:out[j]|=1<<i
  h=H(out);a[h]+=1
  if strong(out):s[h]+=1
 allc[n]=a;sc[n]=s
# Weighted sequence GF: each tournament is one ordered sequence of SCCs.
f=[Counter()for _ in range(7)];f[0][1]=Fraction(1)
for n in range(1,7):
 for h,v in f[n-1].items():f[n][h]+=v
 for m in range(3,n+1):
  for hs,vs in sc[m].items():
   for h,v in f[n-m].items():f[n][h*hs]+=Fraction(vs,factorial(m))*v
 assert dict(f[n])=={h:Fraction(v,factorial(n))for h,v in allc[n].items()}
# Direct exact DP hostile/positive control for the Fibonacci fan.
fans=[]
for k in range(4,13):
 n=k+2;out=[sum(1<<j for j in range(i+1,n))for i in range(n)]
 a=k//2;b=k+1-a
 for t in (a,b,a+b):out[0]^=1<<t;out[t]|=1
 h=H(out);assert h==2**k-1;assert strong(out);fans.append([k,h])
res={'status':'ALL CHECKS PASSED','labeled_tournaments_checked':sum(2**(n*(n-1)//2)for n in range(1,7)),'weighted_SCC_identity_orders':[1,2,3,4,5,6],'fan_DP':fans,'n6_counts':dict(allc[6])}
print(json.dumps(res,indent=2))
