#!/usr/bin/env python3
"""Find the exhaustive seed classes in the orientation-preserving Gewirtz action."""
from collections import Counter
import argparse
from runpy import run_path
from pathlib import Path

mod = run_path(str(Path(__file__).with_name('gewirtz-model.py')))
G=mod['G']
I=tuple(range(56))
def cmp(g,h):return tuple(g[i] for i in h)
tag={I:0}; queue=[I]
for g in queue:
    for t,h in zip([0]*6+[1,2],G):
        k=cmp(g,h); v=tag[g]^t
        if k not in tag:tag[k]=v;queue.append(k)
        else:assert tag[k]==v
assert len(queue)==80640
# The accompanying rectangle certificate proves sign(g)=(-1)^popcount(tag(g)).
K=[g for g in queue if tag[g] in (0,3)]
assert len(K)==40320
D=[g for g in K if mod['derange'](g)]
adj=[{j for j,y in enumerate(mod['H']) if not x&y} for x in mod['H']]
def e(g):return sum(j in adj[i] for i,j in enumerate(g))
print('even derangements by (tag,adjacent moves):',dict(sorted(Counter((tag[g],e(g)) for g in D).items())))
candidates={g for g in D if not tag[g] and e(g)==8}
base=min(candidates)
# Conjugation here may use the full ambient group, not just K, since K is normal.
seen={base};q=[base]; pairs=[(g,mod["inverse"](g)) for g in G]
for u in q:
    for g,gi in pairs:
        v=cmp(g,cmp(u,gi))
        if v not in seen:seen.add(v);q.append(v)
assert seen==candidates and len(q)==2520
rest=[g for g in D if all(u!=v for u,v in zip(g,base))]
assert len(rest)==11212 and Counter(tag[g] for g in rest)=={0:3400,3:7812}
motion=min(sum(x!=y for x,y in enumerate(g)) for g in K if g!=I)
assert motion==48
print("one ambient seed class of size",len(q),"; remaining choices",len(rest))
print("minimum nonidentity even motion:",motion,flush=True)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--seconds",type=float)
parser.add_argument("--witness",default="moore-even-sharp.npy")
args=parser.parse_args()
if args.seconds is not None:
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    arr=np.array(rest,dtype=np.uint8)
    n=len(I)
    src=np.arange(n)*n
    holes=np.ones(n*n,dtype=bool)
    holes[src+I]=holes[src+base]=False
    rows=(arr+src).ravel()
    cols=np.repeat(np.arange(len(rest)),n)
    X=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(n*n,len(rest))).tocsr()[holes].tocsc()
    mask=np.array([tag[g]==3 for g in rest],dtype=int)
    info=milp(c=np.zeros(len(rest)),integrality=np.ones(len(rest)),bounds=Bounds(0,1),
              constraints=[LinearConstraint(X,1,1),LinearConstraint(mask[None,:],7,49)],
              options={"disp":True,"time_limit":args.seconds})
    print(info.message,flush=True)
    if info.x is not None:
        sharp=np.vstack([I,base,arr[info.x>.5]])
        assert sharp.shape==(56,56) and np.all(np.sort(sharp,axis=0)==np.arange(n)[:,None])
        np.save(args.witness,sharp)
        print("Verified a sharp set at",args.witness)
