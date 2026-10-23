import json
N=57
unseen=set(range(1,N));orbits=[]
while unseen:
 a=min(unseen);o=tuple(sorted({a,7*a%N,49*a%N}));orbits.append(o);unseen-=set(o)
used=set();pairs=[]
for o in orbits:
 if o in used:continue
 neg=tuple(sorted(-a%N for a in o))
 if neg==o:raise RuntimeError('self negative')
 used|={o,neg};pairs.append((o,neg))
if len(pairs)!=10:raise RuntimeError('orbit count')
for mask in range(1<<len(pairs)):
 S={a for b,pair in enumerate(pairs) for a in pair[(mask>>b)&1]}
 if len(S)!=28 or S&{-a%N for a in S} or {7*a%N for a in S}!=S:raise RuntimeError('orientation')
 # Every translate has score28; symbolically M1=0.
# Independent finite mod2 elimination for the even principal off-diagonal matrix.
rows=[sum(1<<j for j in range(56) if j!=i) for i in range(56)]
rank=0
for c in range(56):
 pivot=next((r for r in range(rank,56) if rows[r]>>c&1),None)
 if pivot is None:continue
 rows[rank],rows[pivot]=rows[pivot],rows[rank]
 for r in range(56):
  if r!=rank and rows[r]>>c&1:rows[r]^=rows[rank]
 rank+=1
if rank!=56:raise RuntimeError('rank')
print(json.dumps({'N':57,'multiplier':7,'residue_orbits':orbits,'negative_orbit_pairs':pairs,'labelled_rules_checked':1024,'principal_mod2_rank':rank,'augmentation_dimension':56,'Heisenberg_dimension':57},indent=2))
