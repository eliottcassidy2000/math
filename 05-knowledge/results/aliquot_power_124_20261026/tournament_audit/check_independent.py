#!/usr/bin/env python3
"""Independent critical-boundary audit. No producer imports."""
from collections import Counter
from itertools import combinations,permutations,product
from math import prod
import json

checks=0
def require(ok,message):
    global checks
    checks+=1
    if not ok:raise RuntimeError(message)

# Reconstruct graph faces as actual squarefree divisors, starting with literal
# integer divisibility rather than the producer's exponent-box construction.
flag_count=0
for n in range(2,1001):
    divisors={d for d in range(1,n) if n%d==0}
    rest=n;pf=[];omega=0
    for p in range(2,n+1):
        if rest%p:continue
        pf.append(p)
        while rest%p==0:rest//=p;omega+=1
        if rest==1:break
    vertices=[p for p in pf if p in divisors]
    forbidden={frozenset((p,q)) for p,q in combinations(vertices,2) if p*q not in divisors}
    face_divisors=set()
    for bits in product((0,1),repeat=len(vertices)):
        S={p for p,b in zip(vertices,bits) if b}
        if not any(pair<=S for pair in forbidden):face_divisors.add(prod(S))
    flag=face_divisors==divisors
    require(flag==(omega<=2),'literal divisor graph iff boundary')
    flag_count+=flag
require(flag_count==467,'flag census agrees')

# A second feedback oracle minimizes recursively over the last vertex of
# each subset; the producer minimizes over all permutations.
rows=(106,60,89,48,97,68,10)
n=len(rows)
M=[[int(bool(rows[i]&(1<<j))) for j in range(n)] for i in range(n)]
require(all(not M[i][i] for i in range(n)),'no loops')
require(all(M[i][j]+M[j][i]==1 for i,j in combinations(range(n),2)),'valid tournament')
require(all(0<=row<1<<n for row in rows),'no extra vertices')
cost=[0]+[999]*((1<<n)-1);best=[() for _ in cost]
for mask in range(1,1<<n):
    for v in range(n):
        if not mask&(1<<v):continue
        before=mask^(1<<v)
        candidate=cost[before]+sum(M[v][u] for u in range(n) if before&(1<<u))
        if candidate<cost[mask]:cost[mask]=candidate;best[mask]=best[before]+(v,)
require(cost[-1]==5,'subset feedback DP gives5')
feedback=Counter();H=0;Aut=0
for order in permutations(range(n)):
    feedback[sum(M[order[j]][order[i]] for i in range(n) for j in range(i+1,n))]+=1
    H+=all(M[order[i]][order[i+1]] for i in range(n-1))
    Aut+=all(M[i][j]==M[order[i]][order[j]] for i in range(n) for j in range(n))
require(min(feedback)==cost[-1] and sum(feedback.values())==5040,'feedback histogram checks every order')
require(H==127 and Aut==1,'H127 rigid witness')
require(sum(M[best[-1][j]][best[-1][i]] for i in range(n) for j in range(i+1,n))==5,'attaining order')

# Count cycles with every starting point, then divide by length, unlike the
# producer's least-vertex-rooted cycles. Retain full cyclic order information.
cycle_counts={}
for k in (3,5,7):
    rooted=sum(all(M[c[i]][c[(i+1)%k]] for i in range(k)) for c in permutations(range(n),k))
    require(rooted%k==0,'cycle rotation orbit size')
    cycle_counts[k]=rooted//k
triangles=[]
for triple in combinations(range(n),3):
    if all(sum(M[v][w] for w in triple)==1 for v in triple):triangles.append(frozenset(triple))
pairs=sum(not(A&B) for A,B in combinations(triangles,2))
require(sum(cycle_counts.values())==45 and pairs==9,'independent OCF profile1,45,9')
require(1+2*sum(cycle_counts.values())+4*pairs==H,'literal cycle OCF at127')

# F8 via polynomial convolution and descending reduction x^3=x+1,
# independent of the producer's shift-and-reduce multiplication loop.
def multiply(x,y):
    a=[(x>>i)&1 for i in range(3)]
    b=[(y>>i)&1 for i in range(3)]
    c=[0]*5
    for i in range(3):
        for j in range(3):c[i+j]^=a[i]*b[j]
    for d in (4,3):
        if c[d]:c[d]=0;c[d-3]^=1;c[d-2]^=1
    return sum(c[i]<<i for i in range(3))

for x,y,z in product(range(8),repeat=3):
    require(multiply(multiply(x,y),z)==multiply(x,multiply(y,z)),'F8 associativity control')
    require(multiply(x,y^z)==multiply(x,y)^multiply(x,z),'F8 distributivity control')
powers=[1]
for k in range(7):powers.append(multiply(powers[-1],2))
require(powers==[1,2,4,3,6,7,5,1],'marked power table')
line=[powers[r] for r in (1,2,4)]
normal=[powers[(3*r)%7] for r in (1,2,4)]
def span(vectors):
    ans={0}
    for v in vectors:ans|={x^v for x in ans}
    return ans
require(line==[2,4,6] and len(span(line))==4,'Fano line rank2')
require(normal==[3,5,7] and len(span(normal))==8,'normal basis rank3')
require(line[0]^line[1]^line[2]==0 and normal[0]^normal[1]^normal[2]==1,'distinct trace markings')
require([multiply(v,v) for v in normal]==normal[1:]+normal[:1],'normal orbit is Frobenius')
require(powers[2]!=powers[1]^powers[1],'multiplicative marking is not additive')

print(json.dumps({'status':'PASS','checks':checks,'flag_count_2_to1000':flag_count,'H127':{'rows':rows,'H':H,'Aut':Aut,'minimum_feedback':cost[-1],'attaining_order':best[-1],'feedback_histogram':dict(sorted(feedback.items())),'cycle_counts':cycle_counts,'alpha':[1,45,9]},'F8':{'line':line,'line_span_size':len(span(line)),'normal':normal,'normal_span_size':len(span(normal))},'scope':'Independent finite controls for the exact proof audit; no new all-order tc or ordinary spectrum completeness claim.'},sort_keys=True,indent=2))
