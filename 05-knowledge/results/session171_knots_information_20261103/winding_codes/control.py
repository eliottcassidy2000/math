import itertools, json, math, collections
seeds=(153,8345,8337,8369,177,185)
pairs=list(itertools.combinations(range(6),2))
def graph(m):
 g=[[False]*6 for _ in range(6)]
 for b,(i,j) in enumerate(pairs):g[i][j]=bool(m&(1<<b));g[j][i]=not g[i][j]
 return g
out=[]
for m in seeds:
 g=graph(m); counts=[0]*7
 for p in itertools.permutations(range(6)):
  for cut in range(32):
   if all(cut&(1<<i) or g[p[i]][p[i+1]] for i in range(5)):
    counts[1+cut.bit_count()]+=1
 profile=[counts[k]//math.factorial(k) for k in range(1,7)]
 if profile != [29,131,162,77,15,1]:raise Exception(profile)
 for r in range(10):
  f=sum(profile[k-1]*math.prod(range(r-k+1,r+1)) for k in range(1,min(6,r)+1))
  if f!=r**6+12*r**4+16*r*r:raise Exception(('profile',m,r,f))
 for root in range(6):
  seen={root}
  while True:
   nex=seen|{v for u in seen for v in range(6) if g[u][v]}
   if nex==seen:break
   seen=nex
  if len(seen)!=6:raise Exception(('strong',m))
 canon=min(sum((1<<b) for b,(i,j) in enumerate(pairs) if g[p[i]][p[j]]) for p in itertools.permutations(range(6)))
 out.append(dict(seed=m,profile=profile,canonical=canon))
if not all(out[i]['canonical']==out[(i+3)%6]['canonical'] for i in range(6)):raise Exception('antipodal iso')
if len({x['canonical'] for x in out})!=3:raise Exception('class census')
# Full induced seed adjacency, all starting labels and antipodes.
for i,a in enumerate(seeds):
 for j,b in enumerate(seeds):
  if ((a^b).bit_count()==1) != ((i-j)%6 in [1,5]):raise Exception(('adjacency',i,j))
# Decode all errors up to radius 3 for repetition code independently.
errors=[v for v in range(128) if v.bit_count()<=3]
syn=lambda v: tuple(((v>>i)^(v>>6))&1 for i in range(6))
if len(set(map(syn,errors)))!=64:raise Exception('repetition')
print(json.dumps(dict(seed_controls=out,repetition_radius3_syndromes=64),indent=2))
