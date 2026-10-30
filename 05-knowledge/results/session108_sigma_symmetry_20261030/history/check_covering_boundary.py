"""Replay the inherited covering counterexample, exactly over continuous time."""
from fractions import Fraction as F
from itertools import combinations
import json

def check(S):
    denominators={2*v for v in S}
    for a,b in combinations(S,2):denominators.update((a+b,abs(a-b)))
    times={F(k,d) for d in denominators if d for k in range(d+1)}
    def dist(t):
        r=t%1;return min(r,1-r)
    pairs=[(min(dist(v*t) for v in S),t) for t in times]
    m=max(x[0] for x in pairs)
    return m,{'speeds':S,'maximum':str(m),'first_witness':str(min(t for v,t in pairs if v==m)),'critical_points':len(times)}
a,ra=check([1,2,5,6,7,8]);b,rb=check([1,2,3,4,5,42])
if not(a==F(2,13)<F(7,43)==b):raise RuntimeError('inherited boundary mismatch')
print(json.dumps({'status':'PASS','scope':'Exact two finite speed sets; no general LRC result','rows':[ra,rb]},indent=2,sort_keys=True))
