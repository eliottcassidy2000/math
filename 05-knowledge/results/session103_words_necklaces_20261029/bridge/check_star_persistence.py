from star_probe import stars
from pathlib import Path
import json
W=['']; cache={};counts=[];fails=[];birth=[]
for depth in range(8):
 count=0
 for w in W:
  if 'be' not in w and 'ee' not in w: continue
  s=set(stars(w));cache[w]=s;count+=1
  if not s:fails.append(['no_star',w])
  if w[:-1] in cache and w[-2:] not in ['be','ee'] and not s.intersection(cache[w[:-1]]):fails.append(['no_inherited_star',w])
  if w[-2:] in ['be','ee']:
   # Canonical star is anchored at oldest vertex of the top triangle: max vertex before two birth letters.
   anchored=[x for x in s if x[2]==len(w)+1]
   birth.append([w,sorted(s),anchored])
 counts.append([depth,count]);print(counts[-1], 'fails',fails,flush=True)
 W=[w+x for w in W for x in ('ab' if w.endswith(('a','c')) else 'cde')]
Path(__file__).with_name('star_receipt.json').write_text(json.dumps({'counts':counts,'fails':fails,'birth':birth},indent=2))
