"""Independent full-permutation / full-path checks of the attachment hostiles."""
from itertools import permutations, combinations
from collections import Counter
import json

def require(b,s):
    if not b:raise RuntimeError(s)

def matrix(rows):return [[bool(r>>j&1) for j in range(len(rows))] for r in rows]
def maps(a,b,roots=()):
    n=len(a)
    return [p for p in permutations(range(n)) if all(p[x]==y for x,y in roots)
            and all(a[i][j]==b[p[i]][p[j]] for i in range(n) for j in range(n))]
def reverse(a):
    n=len(a)
    return [[a[n-1-j][n-1-i] for j in range(n)] for i in range(n)]
def complement(a):
    n=len(a)
    return [[False if i==j else (a[i][j] if abs(i-j)==1 else not a[i][j]) for j in range(n)] for i in range(n)]
def reorder(a,p):return [[a[i][j] for j in p] for i in p]
def join(a,b):
    n=len(a);m=len(b)
    return [list(row)+[True]*m for row in a]+[[False]*n+list(row) for row in b]
def hp(a):return [p for p in permutations(range(len(a))) if all(a[i][j] for i,j in zip(p,p[1:]))]
def rooted_code(a):
    n=len(a)
    return min(tuple(a[p[i]][p[j]] for i,j in combinations(range(n),2))
               for mid in permutations(range(1,n-1)) for p in [(0,)+mid+(n-1,)])

x=matrix([38,52,24,19,33,12]);y=matrix([26,44,9,48,38,5])
cx,cy=complement(x),complement(y)
source=maps(x,y);target=maps(cx,cy)
first=maps(cx,cy,[(0,0)]);last=maps(cx,cy,[(5,5)])
both=maps(cx,cy,[(0,0),(5,5)])
require(source and target and first and last and not both,'joint root obstruction')
require((reverse(x)==x)==(reverse(y)==y),'same literal colour')
unit=[[False]]
zx=complement(join(join(unit,x),unit));zy=complement(join(join(unit,y),unit))
counts=[len(hp(zx)),len(hp(zy))]
require(counts==[289,253],'independent full-path sandwich counts')

a=matrix([2,4,9,19,7])
presentations={tuple(tuple(row) for row in reorder(a,p)) for p in hp(a)}
require(len(presentations)==9,'nine atom presentations')
hist=Counter(rooted_code(complement(t)) for t in presentations)
literal=sum(reverse(join(u,v))==join(u,v) for u in presentations for v in {tuple(tuple(row) for row in reverse(t)) for t in presentations})
require(literal==9 and sum(v*v for v in hist.values())==11,'reflection is a diagonal coupling')
print(json.dumps({'joint_roots':{'source_isomorphism':source[0],'target_isomorphism':target[0],
    'first_root_isomorphism':first[0],'last_root_isomorphism':last[0],
    'joint_isomorphisms':len(both),'blue':reverse(x)==x,'sandwich_hamiltonian_paths':counts},
    'reflection':{'all_independent_pairs':81,'literal_blue':literal,'wrong_square_sum':11,
                  'rooted_fibre_sizes':sorted(hist.values())}},sort_keys=True,indent=2))
