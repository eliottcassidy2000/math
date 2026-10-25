#!/usr/bin/env python3
from itertools import combinations,product,permutations
from pathlib import Path
import json
from check import arr,H,V,HC,require
P=tuple(product(range(3),repeat=2))
L={frozenset(t) for t in combinations(P,3) if all(sum(v[k] for v in t)%3==0 for k in (0,1))}
require(len(L)==12,'affine line count')
for pair in combinations(P,2): require(sum(set(pair)<=l for l in L)==1,'Steiner pair')
mitres={}
for missing in combinations(P,2):
 support=set(P)-set(missing); blocks={l for l in L if l<=support}
 deg={v:sum(v in l for l in blocks) for v in support}
 require(len(blocks)==5 and sorted(deg.values())==[2]*6+[3],'mitre degrees')
 center=next(v for v,d in deg.items() if d==3)
 outer=[l for l in blocks if center not in l]
 require(len(outer)==2 and not outer[0]&outer[1],'mitre transversals')
 require(all((missing[0][k]+missing[1][k]+center[k])%3==0 for k in (0,1)),'missing pair center')
 mitres[frozenset(missing)]=(center,blocks)
patterns=set(); pairs=[]
for (m,(c,a)),(n,(d,b)) in combinations(mitres.items(),2):
 if c!=d: continue
 require(not m&n and len(a&b)==2,'same-center seam')
 union=frozenset(a|b); require(len(union)==8,'eight-line union')
 patterns.add(union);pairs.append((tuple(sorted(m)),tuple(sorted(n))))
require(len(pairs)==len(patterns)==54,'mitre pair pattern census')
# Independent magic-array enumeration via vector-valued normal form.
layouts=set(); layout_patterns={}
for E,x,z in product(P,repeat=3):
 a=tuple(tuple(arr(E[k],x[k],z[k])[i]%3 for k in (0,1)) for i in range(9))
 if len(set(a))!=9: continue
 blocks=frozenset(frozenset(a[i] for i in line) for line in set(H+V))
 require(blocks<=L,'layout embeds')
 layouts.add(a);layout_patterns[blocks]=layout_patterns.get(blocks,0)+1
require(len(layouts)==432 and set(layout_patterns)==patterns,'independent magic patterns')
require(set(layout_patterns.values())=={8},'D4 fibers')
# Literal grid omitted four lines and integer defects, independently checked.
coords=tuple((j,i) for i in range(3) for j in range(3))
gridlines={frozenset(coords[i] for i in l) for l in set(H+V)}
missing=L-gridlines
omitted=sorted(tuple(coords.index(v) for v in sorted(l)) for l in missing)
for E,x,z in product(range(-3,4),repeat=3):
 a=arr(E,x,z)
 defects=sorted(sum(a[i] for i in line)-3*E for line in omitted)
 require(defects==sorted((3*x,-3*x,3*z,-3*z)),'four missing line defects')
# Nonzero-center prime hostile searched independently and exhaustively up to37.
nonzero=[]
for p in (2,3,5,7,11,13,17,19,23,29,31,37):
 Q={r*r%p for r in range(p)}; found=None
 for E,x,z in product(sorted(Q-{0}),range(p),range(p)):
  a=tuple(v%p for v in arr(E,x,z))
  if len({a[i] for i in HC})==7 and {a[i] for i in HC}<=Q and a[3] not in Q and a[5] not in Q:
   found={'E':E,'x':x,'z':z,'values':a};break
 nonzero.append({'p':p,'first':found})
require(nonzero[-1]['first']['values']==(9,3,28,20,1,19,11,36,30),'F37 hostile')
require(all(r['first'] is None for r in nonzero[:-1]),'small-prime minimality')
result={'status':'PASS','affine_points':9,'affine_lines':12,'mitres':36,'same_center_pairs':54,'magic_layouts':432,'layouts_per_pattern':8,'omitted_grid_cell_triples':omitted,'nonzero_center_both_nonsquare_hostiles':nonzero}
Path(__file__).with_name('geometry-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
