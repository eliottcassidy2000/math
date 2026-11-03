#!/usr/bin/env python3
from itertools import permutations,combinations
from pathlib import Path
import json
from check import require,matrix,encode,relabel,sccs,profile,profile_subset
from search_induced import chain
BASE=Path(__file__).resolve().parent
MASKS=[801,296,673,4712,825,809,817,1584,1589,1593]
n=7;ps=list(permutations(range(n)));pairs=list(combinations(range(n),2));_,expected=chain(3)
a=[matrix(n,t)for t in MASKS]
require(all(len(sccs(t))==1 for t in a),'strong seed')
r=[[encode(relabel(t,p))for p in ps]for t in a]
dist=[[0]*10 for _ in range(10)];witnesses=[]
for i,j in combinations(range(10),2):
 ds=[(t^MASKS[j]).bit_count()for t in r[i]];d=min(ds);dist[i][j]=dist[j][i]=d
 require(d>0,'isomorphic seed collision')
 require((d==1)==((i,j)in expected),'induced chain mismatch')
 if d==1:
  k=ds.index(1);delta=r[i][k]^MASKS[j]
  witnesses.append({'source':i,'target':j,'permutation':list(ps[k]),'relabelled_mask':r[i][k],'flip_pair':list(pairs[delta.bit_length()-1])})
profiles=[]
for t in a:
 p,_=profile(t);require(p==profile_subset(t),'profile algorithm mismatch');profiles.append(p)
out={'status':'PASS','tournament_order':7,'chain_diamonds':3,'masks':MASKS,'point_order':['O','t0','a0','b0','t1','a1','b1','t2','a2','b2'],'edges':sorted(expected),'minimum_quotient_hamming_distances':dist,'edge_witnesses':witnesses,'profiles':profiles,'pair_checks':45,'permutations_per_pair':5040}
(BASE/'chain7_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['status','tournament_order','chain_diamonds','pair_checks','permutations_per_pair']},indent=2))
