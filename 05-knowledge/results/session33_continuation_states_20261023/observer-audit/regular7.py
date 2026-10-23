#!/usr/bin/env python3
import itertools as it, collections,json
E=list(it.combinations(range(7),2));N=(1<<21)-1
incident=[sum(1<<k for k,(a,b) in enumerate(E) if v in (a,b)) for v in range(7)]
negative=[sum(1<<k for k,(a,b) in enumerate(E) if b==v) for v in range(7)]
def arc(m,a,b):
 k=E.index((min(a,b),max(a,b)));return ((m>>k)&1)==(a<b)
def hp(m):return sum(all(arc(m,a,b) for a,b in zip(p,p[1:])) for p in it.permutations(range(7)))
# Full 2^21 scan, distinct from root score-pruned generation.
regular={m for m in range(1<<21) if all(((m^negative[v])&incident[v]).bit_count()==3 for v in range(7))}
left=set(regular); classes={};orbits=[]
while left:
 rep=min(left);h=hp(rep);orbit=set()
 for p in it.permutations(range(7)):
  z=sum(1<<k for k,(a,b) in enumerate(E) if arc(rep,p[a],p[b]));orbit.add(z)
 if not orbit<=regular:raise RuntimeError('orbit nonregular')
 for z in orbit:classes[z]=h
 left-=orbit;orbits.append({'rep':rep,'H':h,'size':len(orbit)})
rows={};checked=0
for m in regular:
 row=collections.Counter()
 for a,b,c in it.combinations(range(7),3):
  if arc(m,a,b)==arc(m,b,c)==arc(m,c,a):
   flip=sum(1<<E.index(e) for e in [(a,b),(a,c),(b,c)])
   if (m^flip) not in regular:raise RuntimeError('not regular after cyclic flip')
   row[classes[m^flip]]+=1
 key=classes[m];r=[row[x] for x in [189,171,175]]
 if key in rows and r!=rows[key]:raise RuntimeError('not lumpable')
 rows[key]=r;checked+=1
if [rows[x] for x in [189,171,175]]!=[[0,14,0],[2,9,3],[0,7,7]]:raise RuntimeError('row mismatch')
print(json.dumps({'status':'PASS','full_masks_scanned':1<<21,'regular_count':len(regular),'orbits':orbits,'all_vertices_lumpability_checked':checked,'rows_H_order_189_171_175':rows},indent=2,sort_keys=True))
