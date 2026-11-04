"""Exact anchor-battery checks for exposed laminar module reversals."""
import itertools,json
from pathlib import Path
from search import modules,edge_mask,classify
from verify import tournament,opposite,endpoints,induced,need
ROOT=Path(__file__).resolve().parent

def submasks(m):
 s=m
 while True:
  yield s
  if s==0:return
  s=(s-1)&m

def read_bank(n):
 d=json.loads((ROOT/f'full_response_n{n}.json').read_text());basis=[a|(c<<n) for a,c in d['basis']]
 polys=[[(s,v) for s,v in zip(basis,r['vector']) if v] for r in d['records']]
 cache={}
 def value(code,a,c):
  key=code,a,c
  if key not in cache:
   mask=a|(c<<n);cache[key]=sum(v for s,v in polys[code] if s&mask==s)
  return cache[key]
 return value

def main():
 records=[];families=anchors=covered=avgchecks=0;counts=[]
 for n in range(1,6):
  V=read_bank(n);local=[0,0,0];full=(1<<n)-1
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code);mods=modules(a);em={m:edge_mask(n,m) for m in mods}
   for f in range(1,1<<len(mods)):
    ms=[m for j,m in enumerate(mods) if f>>j&1]
    if classify(ms)=='overlapping':continue
    maxima=[m for m in ms if not any(m!=s and m&s==m for s in ms)];options=[]
    for B in maxima:
     proper=[m for m in ms if m!=B and m&B==m];union=0
     for m in proper:union|=m
     exposed=B&~union
     if exposed:options.append((B,exposed))
    if not options:covered+=1;local[2]+=1;continue
    target=code
    for m in ms:target^=em[m]
    families+=1;local[0]+=1
    for B,exposed in options:
     vs=[u for u in range(n) if B>>u&1];core=induced(a,vs);S=endpoints(opposite(core))[-1];E=endpoints(core)[-1];index={v:i for i,v in enumerate(vs)}
     exterior=sum(1<<u for u in range(n) if not B>>u&1 and a[u]>>vs[0]&1)
     for v in vs:
      if not exposed>>v&1:continue
      b=(S[index[v]]+E[index[v]])%2;d=(a[v]&B).bit_count();queries=[]
      for u in vs:
       if u!=v:queries.append((exterior|(1<<u)|((1<<v) if a[v]>>u&1 else 0),exterior|(1<<v)))
      if d%2:queries.append((exterior|(1<<v),exterior|(1<<v)))
      if (len(vs)-1+d)%2:queries.append((exterior,exterior|(1<<v)))
      need(len(set(queries))==len(queries),'distinct anchor battery')
      only_B=code^em[B]
      original=sum(V(code,x,c) for x,c in queries);one=sum(V(only_B,x,c) for x,c in queries);multiple=sum(V(target,x,c) for x,c in queries)
      # The n=2 boundary uses the endpoint query only; E1's constant is zero there.
      if len(vs)>=3:need((original-multiple)%4==2*(1+b)%4,'all-order anchor residue')
      need(one==multiple,'exact aggregate invisibility of other modules')
      onequery=(exterior|(1<<v),exterior)
      singles=[V(c,*onequery) for c in (code,only_B,target)]
      need(singles[1]==singles[2],'endpoint query invariant')
      if b:need((singles[0]-singles[2])%4==2,'endpoint imbalance separator')
      if len(vs)==2:need(b==1,'two-vertex boundary handled by endpoints')
      anchors+=1;local[1]+=1
      records.append({'n':n,'source_code':code,'target_code':target,'module_masks':ms,'chosen_outer_module':B,'anchor':v,'endpoint_parity':b,'outdegree_in_outer':d,'query_masks':[list(q) for q in queries],'source_query_sum':original,'single_outer_query_sum':one,'all_reversals_query_sum':multiple,'endpoint_query':list(onequery),'endpoint_counts_source_single_all':singles})
   # Independently check the exact averaging lemma over all outside port masks through n=4.
   if n<=4:
    for A in mods:
     ext=full^A;us=[u for u in range(n) if A>>u&1];target=code^em[A]
     for xout,cout,cy in itertools.product(submasks(ext),submasks(ext),(0,A)):
      s=sum(V(code,xout|(1<<u),cout|cy) for u in us);t=sum(V(target,xout|(1<<u),cout|cy) for u in us)
      need(s==t,'singleton-port averaging lemma');avgchecks+=1
  counts.append({'n':n,'exposed_laminar_families':local[0],'anchor_certificates':local[1],'fully_covered_laminar_families':local[2]})
 return {'status':'PASS','scope':'The proof covers exposed-anchor laminar families at arbitrary order; these are exact bounded controls on complete named responses. Fully covered laminar families are counted but not claimed excluded by this theorem.','exposed_laminar_families':families,'anchor_certificates':anchors,'fully_covered_laminar_families':covered,'exact_averaging_checks':avgchecks,'by_order':counts,'records':records}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
