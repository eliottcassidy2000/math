"""Independent metric-fold audit via common-neighbor obstruction, not row rank."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

E={tuple(sorted((v,v^(1<<i)))) for v in range(8) for i in range(3)}
data=Counter()
for a,b in combinations(range(8),2):
    collapse=lambda v: a if v==b else v
    vertices=sorted({collapse(v) for v in range(8)})
    edges={tuple(sorted((collapse(x),collapse(y)))) for x,y in E}
    if any(x==y for x,y in edges):
        data['unit_loop_impossible']+=1
        continue
    neighbors={x:{y for e in edges if x in e for y in e if y!=x} for x in vertices}
    triples=[(u,v,sorted(neighbors[u]&neighbors[v])) for u,v in combinations(vertices,2)
             if len(neighbors[u]&neighbors[v])>=3]
    if (a^b).bit_count()==2:
        need(bool(triples),'distance-two quotient contains K2,3')
        data['two_circles_three_intersections_impossible']+=1
    else:
        need(not triples,'opposite fold no false K2,3 rejection')
        need(sorted(map(len,neighbors.values()))==[3]*6+[6],'opposite quotient wheel')
        data['wheel_quotient']+=1
need(dict(data)=={'unit_loop_impossible':12,
                 'two_circles_three_intersections_impossible':12,
                 'wheel_quotient':4},'collision universe')
# Centered Walsh coordinates with Gram entries -1/2: this computes the
# entire distance table without using Eisenstein coordinates from producer.
G=[[F(1) if i==j else F(-1,2) for j in range(3)] for i in range(3)]
dist={}
for a,b in combinations(range(8),2):
    delta=[((a>>i)&1)-((b>>i)&1) for i in range(3)]
    dist[a,b]=sum(delta[i]*G[i][j]*delta[j] for i in range(3) for j in range(3))
need(Counter(dist.values())=={F(0):1,F(1):18,F(3):6,F(4):3},'Gram distance table')
need(all(dist[a,b]==1 for a,b in E),'unit edges')
need([e for e,d in dist.items() if not d]==[(0,7)],'exact collision')
# Independent backtracking coloring count: anchored first color, restore k.
def color_count(n, edges,k):
    prev=[[] for _ in range(n)]
    for a,b in edges: prev[max(a,b)].append(min(a,b))
    values=[-1]*n
    def visit(v):
        if v==n: return 1
        total=0
        for c in range(1 if v==0 else k):
            if all(values[u]!=c for u in prev[v]):
                values[v]=c
                total+=visit(v+1)
        values[v]=-1
        return total
    return k*visit(0)
wheel={(0,i) for i in range(1,7)}|{tuple(sorted((i,1+i%6))) for i in range(1,7)}
spindle={(0,1),(0,2),(1,2),(1,3),(2,3),(0,4),(0,5),(4,5),(4,6),(5,6),(3,6)}
counts={'cube2':color_count(8,E,2),'wheel2':color_count(7,wheel,2),
        'wheel3':color_count(7,wheel,3),'spindle3':color_count(7,spindle,3),
        'spindle4':color_count(7,spindle,4)}
need(counts=={'cube2':2,'wheel2':0,'wheel3':6,'spindle3':0,'spindle4':384},'color controls')
print(json.dumps({'status':'PASS','collision_census':dict(data),
                  'color_counts':counts,'method':'K2,3 obstruction and exact Gram distances'},
                 sort_keys=True,indent=2))
