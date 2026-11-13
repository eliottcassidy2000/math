"""Independent compiler audit: no imports of the compiler or its counters."""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import json
BASE=Path(__file__).resolve().parent
SOURCE=BASE.parent/'compiler'/'six_state_audit.json'

def require(ok,why):
    if not ok: raise RuntimeError(why)

def det(mat):
    a=[[F(x) for x in row] for row in mat]
    d=F(1)
    for j in range(len(a)):
        pivot=next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None: return F(0)
        if pivot!=j:
            a[pivot],a[j]=a[j],a[pivot];d=-d
        u=a[j][j];d*=u
        for k in range(j,len(a)): a[j][k]/=u
        for i in range(j+1,len(a)):
            u=a[i][j]
            for k in range(j,len(a)): a[i][k]-=u*a[j][k]
    return d

def exact_state(rows,r):
    # Ordered partitions into 1,2,3 bins are represented uniquely by a
    # permutation and weakly increasing cut locations. Check edges only
    # inside each bin. No subset-DP/convolution is used.
    n=len(rows)
    h=e=s2=p2=e2=s3=0
    for order in permutations(range(n)):
        bad=[0]
        for u,v in zip(order,order[1:]):
            bad.append(bad[-1]+int(not (rows[u]>>v&1)))
        def path(l,u): return u-l<=1 or bad[u-1]==bad[l]
        if path(0,n):
            h+=1;e+=int(order[-1]==r)
        for cut in range(n+1):
            if path(0,cut) and path(cut,n):
                p2+=1
                if cut:
                    e2+=int(order[cut-1]==r)
                    s2+=int(order[0]==r)
        if order[0]==r:
            for left in range(1,n+1):
                for right in range(left,n+1):
                    s3+=int(path(0,left) and path(left,right) and path(right,n))
    return [h,e,s2,p2,e2,s3]

def ready(rows,r):
    n=len(rows)
    return [row|((1<<n) if i==r else 0) for i,row in enumerate(rows)]+[((1<<n)-1)^(1<<r)]

def brute_h(rows):
    return sum(all(rows[u]>>v&1 for u,v in zip(order,order[1:])) for order in permutations(range(len(rows))))

data=json.loads(SOURCE.read_text())
require(det(data['hankel'])==552960,'independent rational Hankel determinant')
require(det(data['difference_hankel'])==811008,'independent affine determinant')
controls=[]
for row in data['literal_controls']:
    got=exact_state(row['parent'],row['root'])
    require(got==row['state'],'permutation/cut six-state count '+row['word'])
    count=brute_h(ready(row['parent'],row['root']))
    require(count==row['ready_count'],'brute ready-child full count '+row['word'])
    require(got[4]>=1,'positive two-bin endpoint count')
    controls.append({'word':row['word'],'order':row['order'],'state':got,'ready_full_H':count})
record={'status':'PASS','method':'independent Fraction elimination and literal permutation/weak-cut enumeration; no imported compiler code','hankel_det':552960,'affine_difference_det':811008,'literal_controls':controls}
(BASE/'compiler_independent_receipts.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':'PASS','independent_literal_controls':len(controls),'determinants':[552960,811008]}))
