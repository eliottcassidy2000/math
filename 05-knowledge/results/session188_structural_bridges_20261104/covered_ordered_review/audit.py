"""Independent physical-name replay of ordered-block separator receipts."""
from pathlib import Path
import importlib.util,itertools,json,hashlib
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'covered_laminar'
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
data=json.loads((SRC/'ordered.json').read_text());v.need((SRC/'ordered.json').read_bytes()==(SRC/'ordered.optimized.json').read_bytes(),'original modes agree')
def build(blocks,order):
 sizes=list(map(len,blocks));off=[sum(sizes[:i]) for i in range(len(sizes))];rows=[0]*sum(sizes)
 for i,a in enumerate(blocks):
  for u,r in enumerate(a):rows[off[i]+u]|=r<<off[i]
 for pos,i in enumerate(order):
  for j in order[pos+1:]:
   for u in range(sizes[i]):rows[off[i]+u]|=((1<<sizes[j])-1)<<off[j]
 return rows,off
def lift(rows):
 n=len(rows);return [r|(1<<n) for r in rows]+[1<<(n+1),(1<<n)-1]
ordered=0;coordinates=0
for r in data['ordered_blocks']['records']:
 sizes=r['sizes'];variant=r['variant'];k=len(sizes);a=[v.tournament(n,(17*i+variant*5+1)%(1<<(n*(n-1)//2))) for i,n in enumerate(sizes)];b=[v.opposite(x) if (i+variant)%2 else x for i,x in enumerate(a)]
 T,off=build(a,list(range(k)));U,_=build(b,r['order']);cut=r['cut'];prefix=sum(((1<<sizes[i])-1)<<off[i] for i in range(cut));X=v.extend(T,0,prefix);Y=v.extend(U,0,prefix);p=v.profile(X);q=v.profile(Y)
 v.need(all(y>=x for x,y in zip(p,q)),'independent complete profile domination');v.need(q[1]-p[1]==r['local_difference'],'independent local difference')
 outside=1<<(len(T)+1);f=v.hp(v.extend(lift(T),outside,prefix|outside));g=v.hp(v.extend(lift(U),outside,prefix|outside));v.need(g-f==r['ambient_difference'],'independent strong ambient literal count')
 ordered+=1;coordinates+=len(p)-1
laminar=0
for r in data['covered_laminar']['records']:
 T=v.tournament(r['n'],r['source_code']);U=list(T)
 for module in r['modules']:
  for i,j in itertools.combinations(range(len(T)),2):
   if module>>i&1 and module>>j&1:U[i]^=1<<j;U[j]^=1<<i
 a,c=r['query'];values=[v.hp(v.extend(rows,a,c)) for rows in (T,U)];v.need(values==r['values'] and values[1]>values[0],'independent original-module XOR physical query');laminar+=1
print(json.dumps({'status':'PASS','ordered_block_records':ordered,'profile_coordinates':coordinates,'laminar_family_records':laminar,'method':'Separately implemented predecessor Hamilton DP, independent reconstruction of physical queries and XOR targets, independent ordered-cover convolution.','source_sha256':{p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ('ORDERED_BLOCKS.md','verify_ordered.py','ordered.json','ordered.optimized.json')}},sort_keys=True,indent=2))
