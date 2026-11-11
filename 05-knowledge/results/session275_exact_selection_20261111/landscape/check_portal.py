from itertools import permutations
from experiment import hcount
from itertools import product
from collections import defaultdict
def canon(a):
 n=len(a); groups=defaultdict(list)
 for i,r in enumerate(a):groups[r.bit_count()].append(i)
 candidates=[]
 for ps in product(*(permutations(groups[d]) for d in sorted(groups))):
  p=sum(ps,());r=tuple(sum(1<<j for j in range(n) if a[p[i]]>>p[j]&1) for i in range(n));candidates.append((sum(x<<(8*i) for i,x in enumerate(r)),r))
 return min(candidates)[1]
import json,pathlib
P=pathlib.Path(__file__).parent
r=json.loads((P/'PORTAL63.json').read_text());chain=r['chain'];checks=[]
def brute(a):return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in permutations(range(len(a))))
def record(a,b,flips):
 common=list(a)
 for i,j in flips:common[i]&=~(1<<j);common[j]&=~(1<<i)
 path=next(p for p in permutations(range(8)) if all(common[p[i]]>>p[i+1]&1 for i in range(7)))
 ha,hb=brute(a),brute(b)
 if abs(ha-hb)!=2:raise ValueError('unit move')
 checks.append({'beforeH':ha,'afterH':hb,'common_path':path,'flips':flips})
for c in chain:
 a,b=c['before'],c['after_raw'];record(a,b,c['flips'])
 if canon(b)!=tuple(c['after_canonical']):raise ValueError('relabel')
upper=[144,5,17,7,10,31,191,62];record(r['seed'],upper,[(0,7),(4,6)])
for a,deleted in [(r['found'],7),(upper,6)]:
 if a[deleted].bit_count()!=7:raise ValueError('not source')
 keep=[i for i in range(8) if i!=deleted];small=[sum(1<<j for j,x in enumerate(keep) if a[v]>>x&1) for v in keep]
 if brute(small)!=brute(a):raise ValueError('not neutral')
 checks.append({'delete_source':deleted,'rows7':small,'H':brute(small)})
print(json.dumps({'status':'PASS','checks':checks},indent=2))
