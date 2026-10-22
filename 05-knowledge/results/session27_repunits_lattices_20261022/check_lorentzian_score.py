import random,json
from itertools import combinations
rng=random.Random(27025);checks=0
cases=[]
for sample in range(64):
 out=[0]*25
 for u,v in combinations(range(25),2):
  a,b=(u,v) if rng.randrange(2) else(v,u);out[a]|=1<<b
 s=[x.bit_count() for x in out]
 c=sum(all(((out[u]>>v)&1,(out[v]>>w)&1,(out[w]>>u)&1)) or all(((out[v]>>u)&1,(out[w]>>v)&1,(out[u]>>w)&1)) for u,v,w in combinations(range(25),3))
 norm=sum(x*x for x in s)-70**2
 if norm!=-2*c or (sum(s)-70)%2:raise ValueError('identity or lattice')
 checks+=1
 if sample<3:cases.append({'scores':s,'c3':c,'norm':norm})
s=[4,4,4,3,3,3,3,4]+list(range(24,7,-1))
if len(s)!=25 or sum(s)!=300 or sum(x*x for x in s)-4900!=-40:raise ValueError('padding')
if sum(j*j for j in range(25))!=4900:raise ValueError('cannonball')
print(json.dumps({'random_controls':checks,'seed':27025,'examples':cases,'transitive_norm':0,'padded_scores':s,'padded_norm':-40,'padded_c3':20},indent=2))
