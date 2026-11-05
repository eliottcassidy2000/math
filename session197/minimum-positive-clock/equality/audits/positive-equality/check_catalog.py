from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
core=[8,5,9,50,39,327,159,63,223]
cat=json.loads((P.parent.parent/'unit-clocks-through-nine/catalog.json').read_text())['catalog']
candidates=[x for x in cat if x['A']==72 and x['beta']==1]
if len(candidates)!=1:raise RuntimeError('catalog amplitude72 cardinality')
def maps(A,B,ra,rb):
 n=len(A);order=[ra]+sorted((u for u in range(n) if u!=ra),key=lambda u:sum(x.bit_count()==A[u].bit_count() for x in B));ans=[];f={};used=set()
 def go(k):
  if k==n:ans.append([f[i] for i in range(n)]);return
  u=order[k]
  for v in ([rb] if k==0 else range(n)):
   if v in used or A[u].bit_count()!=B[v].bit_count():continue
   if any(((A[u]>>w)&1)!=((B[v]>>f[w])&1) for w in f):continue
   f[u]=v;used.add(v);go(k+1);used.remove(v);del f[u]
 go(0);return ans
x=candidates[0];isos=maps(core,x['rows'],4,x['root'])
if len(isos)!=1:raise RuntimeError('unique rooted catalog isomorphism')
autos=[]
for r in range(9):autos.extend(maps(core,core,4,r))
if autos!=[list(range(9))]:raise RuntimeError('independent unmarked rigidity')
print(json.dumps({'status':'PASS','catalog_matches':len(candidates),'core_to_catalog':isos[0],'automorphisms':autos},indent=2,sort_keys=True))
