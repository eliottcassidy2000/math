#!/usr/bin/env python3
from itertools import permutations,combinations
from pathlib import Path
import json
from check import require,matrix,encode,relabel,sccs,profile,profile_subset
from search_induced import chain
BASE=Path(__file__).resolve().parent
MASKS=[17928585,9314500,10846530,9320707,1712324,3244354,1718531,54952202,25475209,3342737,54542928,67668801,83954712,54804816,67684801,84478232]
n=8;ps=list(permutations(range(n)));pairs=list(combinations(range(n),2));_,expected=chain(5)
a=[matrix(n,t)for t in MASKS]
require(all(len(sccs(t))==1 for t in a),'strong seed')
r=[[encode(relabel(t,p))for p in ps]for t in a]
dist=[[0]*16 for _ in range(16)];witnesses=[]
for i,j in combinations(range(16),2):
 ds=[(t^MASKS[j]).bit_count()for t in r[i]];d=min(ds);dist[i][j]=dist[j][i]=d
 require(d>0,'isomorphic seed collision')
 require((d==1)==((i,j)in expected),'induced chain mismatch')
 if d==1:
  k=ds.index(1);delta=r[i][k]^MASKS[j]
  witnesses.append({'source':i,'target':j,'permutation':list(ps[k]),'relabelled_mask':r[i][k],'flip_pair':list(pairs[delta.bit_length()-1])})
profiles=[]
for t in a:
 p,_=profile(t);require(p==profile_subset(t),'profile algorithm mismatch');profiles.append(p)
out={'status':'PASS','tournament_order':8,'chain_diamonds':5,'masks':MASKS,'point_order':['O']+[f'{s}{i}' for i in range(5) for s in ['t','a','b']],'edges':sorted(expected),'minimum_quotient_hamming_distances':dist,'edge_witnesses':witnesses,'profiles':profiles,'pair_checks':120,'permutations_per_pair':40320}
(BASE/'chain8_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['status','tournament_order','chain_diamonds','pair_checks','permutations_per_pair']},indent=2))
