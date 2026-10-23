"""Exact finite checks; all requirements active under python -O."""
from fractions import Fraction as F
from itertools import permutations
from math import factorial
import json

def require(x,msg):
    if not x: raise RuntimeError(msg)
def partitions(n, cap=None):
    if not n: yield (); return
    for a in range(min(n,cap or n),1-1,-1):
        for tail in partitions(n-a,a): yield (a,)+tail

def contents(lam,t):
    out=F(1)
    for i,row in enumerate(lam):
        for j in range(row): out*=t+j-i
    return out

def cycles(p):
    seen=set(); k=0
    for a in range(len(p)):
        if a not in seen:
            k+=1; b=a
            while b not in seen: seen.add(b); b=p[b]
    return k

def choose(t,k):
    a=F(1)
    for i in range(k): a*= (t-i)/F(i+1)
    return a
rows=[]; checks=0
for n in range(1,13):
    for j in range(0,4*(n+2)+1):
        t=F(j,4)
        values=[contents(lam,t) for lam in partitions(n)]
        predicted=(t.denominator==1 or t>=n-1)
        require(all(x>=0 for x in values)==predicted,('PSD classification',n,t)); checks+=1
    if n<=7:
        for t in [F(1,2),F(3,2),F(7,2),F(3)]:
            s=sum(((-1)**(n-cycles(p)))*t**cycles(p) for p in permutations(range(n)))/factorial(n)
            require(s==choose(t,n),('signed Burnside',n,t)); checks+=1
for t in [F(1,2),F(3,2),F(7,2),F(57,2),F(113,2),F(115,2),F(57)]:
    if t.denominator==1:
        rows.append({'t':str(t),'first_negative':None,'first_zero_exterior':int(t)+1})
    else:
        k=t.numerator//t.denominator+2
        require(all(choose(t,i)>0 for i in range(k)),('earlier positive',t)); checks+=1
        require(choose(t,k)<0,('first negative',t)); checks+=1
        rows.append({'t':str(t),'first_negative':k,'negative_trace':str(choose(t,k))})
print(json.dumps({'checks':checks,'finite_psd_grid':'n=1..12; t=j/4 from0 to n+2','signed_permutation_checks':'n=1..7','examples':rows},indent=2))
