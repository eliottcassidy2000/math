import itertools as it,json
from pathlib import Path
import audit as a
result={}
for m,expected in ((40,(5,1)),(41,(13,3))):
 paths=[p for p in it.permutations(range(5)) if all(a.arc(m,5,p[i],p[i+1]) for i in range(4))]
 aut=[p for p in it.permutations(range(5)) if a.induce(m,5,p)==m]
 anti=[p for p in it.permutations(range(5)) if a.induce(m,5,p)==(1023^m)]
 orbit=lambda p:min(tuple(g[v] for v in p) for g in aut)
 reps={orbit(p) for p in paths};fix=[]
 for g in anti:
  fixed=[p for p in reps if orbit(tuple(g[v] for v in p[::-1]))==p]
  fix.append(len(fixed))
 a.require((len(reps),fix[0])==expected,'mirror seed expected')
 a.require(len(set(fix))==1,'antiisomorphism choice independence')
 result[m]={'tile_count':len(reps),'mirror_fixed_each_anti':fix}
for r in range(1,9):
 blue=0;mass=0
 for w in it.product((40,41),repeat=r):
  mass+=__import__('math').prod(result[x]['tile_count'] for x in w)
  if w==w[::-1]:
   blue+=__import__('math').prod(result[x]['tile_count'] for x in w[:r//2])*(result[w[r//2]]['mirror_fixed_each_anti'][0] if r%2 else 1)
 a.require(mass==18**r,'full mass')
 a.require(blue==18**(r//2)*(4 if r%2 else 1),'blue family formula')
 result['r'+str(r)]={'mass':mass,'blue':blue,'black':mass-blue}
Path(__file__).with_name('mirror-results.json').write_text(json.dumps(result,indent=2)+'\n')
print('MIRROR PASS core anti choices and all words r1..8')
