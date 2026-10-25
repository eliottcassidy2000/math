import math,itertools,json,time,sys
from pathlib import Path
N=int(sys.argv[1]) if len(sys.argv)>1 else 2000
start=time.monotonic();A=[set() for _ in range(N+1)]
for k in range(2,math.isqrt(2*N)+1):
 s=k*k
 for a in range(max(1,s-N),min(N,(s-1)//2)+1):
  b=s-a;A[a].add(b);A[b].add(a)
hits={3:None,4:None};seen={3:set(),4:set()};eligible={3:0,4:0};checked={3:0,4:0}
for a in range(1,N):
 for b in range(a+1,N+1):
  C=A[a]&A[b]
  for k in [3,4]:
   if len(C)<k:continue
   eligible[k]+=1
   for R in itertools.combinations(sorted(C),k):
    if R in seen[k]:continue
    seen[k].add(R);checked[k]+=1
    L=set.intersection(*(A[r] for r in R))
    if len(L)<k:continue
    L=sorted(L)[:k];item=(max(L+list(R)),L,list(R))
    if hits[k] is None or item<hits[k]:hits[k]=item
result={'N':N,'smallest_k33_in_universe':hits[3],'smallest_k44_in_universe':hits[4],'eligible_seed_pairs':eligible,'unique_neighbor_tuples_checked':checked,'seconds':time.monotonic()-start}
print(json.dumps(result,indent=2));json.dump(result,open(Path(__file__).resolve().parent.parent/f'closure_{N}.json','w'),indent=2)
