"""Independent exact audit, no producer imports; checks survive python -O."""
from fractions import Fraction as F
from math import isqrt, prod
from pathlib import Path
import json

def require(test, label):
    if not test: raise RuntimeError(label)

def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    m=(3*x*x-25)/(2*y) if P==Q else (v-y)/(u-x)
    X=m*m-x-u
    return X,m*(x-X)-y

def height(P):
    x=P[0];return max(abs(x.numerator),x.denominator)

residues={}
for modulus in (3,16):
    sq={x*x%modulus for x in range(modulus)}
    good=[d for d in range(modulus) if any((b-5*d*d)%modulus in sq and (b+5*d*d)%modulus in sq for b in sq)]
    residues[str(modulus)]=good
require(residues=={'3':[0],'16':[0,4,8,12]}, 'denominator congruences')
solutions=[]; pairs=[]
for u in range(1,721):
    if 720%u:continue
    v=720//u
    if v<=u or (v-u)%2:continue
    a=(v-u)//2;b=(v+u)//2;c=isqrt(b*b+720)
    pairs.append([u,v])
    if c*c==b*b+720:solutions.append([a,b,c])
require(solutions==[[31,41,49]], 'unbounded-numerator factor enumeration')
P=(F(25,4),F(75,8));Q=(F(1681,144),F(62279,1728)); R=add(Q,Q)
require(add(P,P)==(Q[0],-Q[1]), 'sign in doubling')
require(all(y*y==x*x*x-25*x for x,y in (P,Q,R)), 'curve equation')
require(Q[0]>10 and R[0]<10, 'positive clique chamber failure')
# Halving is not unique: translation by rational 2-torsion preserves doubling.
T=(F(0),F(0));PT=add(P,T)
require(add(PT,PT)==add(P,P), 'fourfold nonunique halving control')
minusP=(P[0],-P[1])
halves=[add(minusP,T) for T in (None,(F(0),F(0)),(F(5),F(0)),(F(-5),F(0)))]
def square(z):return z>=0 and isqrt(z.numerator)**2==z.numerator and isqrt(z.denominator)**2==z.denominator
require([S[0] for S in halves]==[F(25,4),F(-4),F(45),F(-5,9)], 'independent chord-law halves')
require(all(add(S,S)==Q for S in halves), 'all four halves')
require(all(not all(square(S[0]+j) for j in (-5,0,5)) for S in halves), 'Q not in 4E')
require((Q[0]**2+25)**2 < 40*Q[0]*(Q[0]**2-25),'triangle-map positivity hostile')

# Producer graph rows are treated only as supplied finite data, never trusted counts.
data=json.loads((Path(__file__).resolve().parent.parent/'descent/CHECK.json').read_text())
from functools import lru_cache
@lru_cache(None)
def paths(rows, mask, start):
    if mask==(1<<start):return 1
    rest=mask^(1<<start)
    return sum(paths(rows, rest, v) for v in range(len(rows)) if rest>>v&1 and rows[start]>>v&1)
def pathcount(rows):
    rows=tuple(rows);mask=(1<<len(rows))-1
    return sum(paths(rows,mask,v) for v in range(len(rows)))
seeds={}
for target,record in data['seeds'].items():
    rows=record['masks'];n=len(rows)
    require(all(0<=r<(1<<n) and not(r>>i&1) for i,r in enumerate(rows)), 'row bounds')
    require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)), 'edge completeness')
    h=pathcount(rows);require(h==int(target),'independent seed count')
    seeds[target]={'order':n,'H':h}
# Independent unbounded-coin DP over a finite target universe, no producer formula.
limit=100000
possible={1}
for h in sorted(map(int,seeds)):
    for n in range(1,limit//h+1):
        if n in possible:possible.add(n*h)
smooth=[]
for n in range(1,limit+1,2):
    m=n
    for p in (3,5,7,31,41):
        while m%p==0:m//=p
    if m==1:
        smooth.append(n);require((n in possible)==(n not in (7,21)),'smooth closure finite control')
# Separate one-dimensional optimization over every 49/343 recipe up to c=200.
carry_cost=[]
for c in range(201):
    candidates=[8*v+11*t for v in range(c//2+1) for t in range(c//3+1) if 2*v+3*t==c]
    cost=min(candidates) if candidates else None
    q,r=divmod(c,3)
    formula=None if c==1 else 11*q+(0,5,8)[r]
    require(cost==formula, 'carry cost formula')
    carry_cost.append(cost)
require(all(carry_cost[c+3]==carry_cost[c]+11 for c in range(198) if c!=1),'carry periodicity')
print(json.dumps({'carry_cost_cases':len(carry_cost),'independent_seed_counts':seeds,'smooth_targets_through_100000':len(smooth),'denominator_residues':residues,'factor_pairs_checked':len(pairs),'unique_d12_positive_roots':solutions,'P':list(map(str,P)),'Q':list(map(str,Q)),'double_Q':list(map(str,R)),'x_height_P_Q_doubleQ':list(map(height,(P,Q,R))),'torsion_translate_P':list(map(str,PT))},indent=2,sort_keys=True))
