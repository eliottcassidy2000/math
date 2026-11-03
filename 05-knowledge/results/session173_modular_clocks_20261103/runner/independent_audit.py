from fractions import Fraction as F
from itertools import product, combinations
from math import gcd, ceil
import json

def check(x,msg):
 if not x: raise RuntimeError(msg)
def norm(x):
 x=x%1
 return min(x,1-x)
def masks(d,E,t,delta):
 return [{q for q in range(d) if norm(w*(t+F(q,d)))<delta} for w in E]
# Exhaust all membership patterns for path, star and cycle; no assert statements.
npatterns=0
for e in range(1,13):
 for s in product([False,True],repeat=e):
  k=sum(s)
  for edges in [[(i,i+1) for i in range(e-1)],[(0,i) for i in range(1,e)]]:
   selected=sum(s[i] and s[j] for i,j in edges)
   components=k-selected
   check(components>=bool(k),'tree bound')
   check(components-bool(k)>=0,'tree defect')
  if e>=3:
   selected=sum(s[i] and s[(i+1)%e] for i in range(e))
   val=k-selected+all(s)
   runs=sum(s[i] and not s[(i-1)%e] for i in range(e))
   gap=runs-1 if 0<k<e else 0
   check(val-bool(k)==gap>=0,'cycle defect')
  npatterns+=1
capacity_checks=0
for d in range(2,22):
 for w in range(-2*d,2*d+1):
  for t in [F(0),F(1,2),F(1,3),F(1,d),F(1,2*d),F(2,3*d)]:
   for delta in [F(1,15),F(1,3),F(1,2),F(1,d)]:
    b=masks(d,[w],t,delta)[0];g=gcd(d,w)
    check(len(b)<=g*ceil(2*delta*d/g),'capacity')
    capacity_checks+=1
# Independent exact fourteen-speed examples, including every good and bad shift.
cases=[(9,8,[4,2,5,3,7,13],81,[{0,2},{0},{0,7},{0,3,6},{5},{2,4}],[1,8]),(8,8,[4,2,5,3,7,13],72,[{0,2,4,6},{0,4},{3},{0},{1},{3}],None),(8,7,[1,3,5,7,13,44,527],64,[{0},{0},{3},{1},{6},set(),{2}],[4,5,7]),(9,6,[1,3,5,7,13,44,527,8],63,[{0},{0,3,6},{7},{5},{4},{6},{3},{1}],[2,8])]
for d,m,E,den,expected,good in cases:
 C=[d*j for j in range(1,m+1)];t=F(1,den);delta=F(1,15)
 check(len(C+E)==len(set(C+E))==14,'14 distinct moving speeds')
 B=masks(d,E,t,delta);check(B==expected,'example masks')
 union=set().union(*B);found=[q for q in range(d) if q not in union]
 if good is not None:check(good==found,'example good list')
 for q in range(d):check((q in found)==all(norm(v*(t+F(q,d)))>=delta for v in C+E),'exact witness iff')
# Translation identity over signed speeds, phase-aligned shifts.
translations=0
for d in range(2,20):
 for w in range(1,d):
  if gcd(w,d)!=1:continue
  for a in [-2,-1,1,2]:
   t=F(1,d*d);h=a%d
   B,Bp=masks(d,[w,w+a*d],t,F(1,5))
   check(Bp=={(q-h*pow(w,-1,d))%d for q in B},'translation')
   translations+=1
print(json.dumps(dict(membership_patterns=npatterns,capacity_checks=capacity_checks,exact_fourteen_speed_examples=len(cases),translation_checks=translations,status='pass'),sort_keys=True))
