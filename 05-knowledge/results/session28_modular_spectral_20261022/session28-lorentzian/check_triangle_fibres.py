from fractions import Fraction as F
from math import factorial,comb
from itertools import combinations
import json

def dp_h(n,ts):
 p=[[F(1,2) if i!=j else F(0) for j in range(n)] for i in range(n)]
 for r,t in enumerate(ts):
  a,b,c=3*r,3*r+1,3*r+2
  for u,v in [(a,b),(b,c),(c,a)]:p[u][v]=(1+t)/2;p[v][u]=(1-t)/2
 dp=[[F(0)]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=F(1)
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    old=mask^(1<<v)
    if old:dp[mask][v]=sum((dp[old][u]*p[u][v] for u in range(n) if old>>u&1),F(0))
 value=sum(dp[-1]);scores=[sum(row) for row in p]
 assert scores==[F(n-1,2)]*n
 cyc=sum(p[a][b]*p[b][c]*p[c][a]+p[b][a]*p[c][b]*p[a][c] for a,b,c in combinations(range(n),3))
 assert cyc==F(comb(n,3),4)+F(3,4)*sum(t*t for t in ts)
 return value,cyc

def formula(n,ts):
 coeff=[F(1)]
 for t in ts:
  coeff=[coeff[j]+(6*t*t*coeff[j-1] if j else 0) for j in range(len(coeff))]+[6*t*t*coeff[-1]]
 return sum(c*factorial(n-2*j) for j,c in enumerate(coeff))/2**(n-1)
out=[]
for n,ts in [(3,[F(1)]),(4,[F(1,2)]),(6,[F(1),F(-1,3)]),(7,[F(2,3),F(3,5)]),(9,[F(1),F(1),F(1)]),(9,[F(1,2),F(-2,3),F(4,5)])]:
 a,c=dp_h(n,ts);b=formula(n,ts);assert a==b
 out.append({'n':n,'biases':list(map(str,ts)),'expected_H':str(a),'expected_c3':str(c)})
n=25;k=8;ratio=formula(n,[F(1)]*k)/F(factorial(n),2**(n-1))
print(json.dumps({'checks':out,'n25_eight_triangles_ratio':str(ratio),'n25_ratio_decimal':float(ratio)},indent=2))
