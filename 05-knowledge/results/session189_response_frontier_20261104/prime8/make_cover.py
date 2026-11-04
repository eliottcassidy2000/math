from itertools import product,permutations
from collections import defaultdict
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

def canonical(T):
 groups=defaultdict(list)
 for i,r in enumerate(T):groups[r.bit_count()].append(i)
 best=None;witness=None
 for pp in product(*(permutations(groups[d]) for d in sorted(groups))):
  p=sum((tuple(v) for v in pp),());R=tuple(sum(((T[i]>>j)&1)<<k for k,j in enumerate(p)) for i in p)
  if best is None or R<best:best=R;witness=p
 return best,witness

def main():
 current=[()];levels=[dict(order=0,representatives=[[]],extension_records=[])]
 for n in range(1,9):
  keys=set();temporary=[]
  for pi,T in enumerate(current):
   for mask in range(1<<(n-1)):
    G=tuple(T[i]|((1-((mask>>i)&1))<<(n-1)) for i in range(n-1))+(mask,)
    canon,p=canonical(G);keys.add(canon);temporary.append((canon,p))
  current=sorted(keys);index={T:i for i,T in enumerate(current)}
  records=[[index[T],sum(i<<(3*j) for j,i in enumerate(p))] for T,p in temporary]
  levels.append(dict(order=n,representatives=current,extension_records=records));print(n,len(current),len(records),flush=True)
 (ROOT/'cover.json').write_text(json.dumps(dict(levels=levels),separators=(',',':'))+'\n')
if __name__=='__main__':main()
