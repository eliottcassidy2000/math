import sys,json,itertools,collections,functools,pathlib,hashlib
sys.path.insert(0,str(pathlib.Path(__file__).parent/'sources'))
from verify import Field,magic,valid,prime_powers
import sympy as sp
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def graph(F,V,zero=True):
 return tuple(sum(1<<j for j,y in enumerate(V) if i!=j and F.add(x,y) in F.squares and (zero or F.add(x,y)!=0)) for i,x in enumerate(V))
@functools.lru_cache(None)
def stats(adj):
 n=len(adj);full=(1<<n)-1
 @functools.lru_cache(None)
 def f(mask,last):
  if mask==full:return 1
  return sum(f(mask|1<<v,v) for v in range(n) if adj[last]>>v&1 and not mask>>v&1)
 counts=[f(1<<v,v) for v in range(n)]
 path=[]
 if sum(counts):
  v=next(v for v in range(n) if counts[v]);path=[v];mask=1<<v
  while mask!=full:
   v=next(w for w in range(n) if adj[v]>>w&1 and not mask>>w&1 and f(mask|1<<w,w));path.append(v);mask|=1<<v
 A=sp.Matrix([[int(adj[i]>>j&1) for j in range(n)] for i in range(n)])
 return {'edges':sum(a.bit_count() for a in adj)//2,'degrees':sorted(a.bit_count() for a in adj),'hp_unoriented':sum(counts)//2,'hp_witness_indices':path,'charpoly_coefficients':[int(x) for x in A.charpoly().all_coeffs()]}
def record(F,m):
 V=sorted(m)
 return {'matrix':m,'vertices':V,'center':m[4],'graphs':{str(z):stats(graph(F,V,z)) for z in (True,False)}}
def grids(F):
 for c in (0,1):
  for u in ([1] if c==0 else range(F.q)):
   for v in range(F.q):
    m=magic(F,c,u,v)
    if valid(F,m):yield m
out={'scope':'all normalized valid grids for odd prime powers q<=101; zero center normalized u=1, nonzero center normalized c=1','fields':[]}
for q,p,n in prime_powers(102):
 F=Field(p,n);ms=list(grids(F)); groups=collections.Counter();examples={}
 for m in ms:
  V=sorted(m)
  key=tuple((z,stats(graph(F,V,z))['edges'],stats(graph(F,V,z))['hp_unoriented'],tuple(stats(graph(F,V,z))['charpoly_coefficients'])) for z in (True,False))
  groups[key]+=1;examples[key]=record(F,m)
 row={'q':q,'p':p,'n':n,'modulus':F.modulus,'normalized_magic_count':len(ms),'spectral_hp_groups':[{'count':v,**examples[k]} for k,v in groups.items()]}
 if q in (17,19,23,25,27,31,43,47,67):
  # first nine square values: a sufficient certificate, not a census of all subsets.
  V=sorted(F.squares)[:9];row['first_nine_square_graph']={'vertices':V,'graphs':{str(z):stats(graph(F,V,z)) for z in (True,False)}}
 out['fields'].append(row);print(q,len(ms),len(groups),flush=True)
pathlib.Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
# Supplemental characteristic-three witnesses, pinned to prior census.
C=json.loads(pathlib.Path(__file__).with_name('sources').joinpath('census.json').read_text())
for q in (729,2187):
 r=next(r for r in C['rows'] if r['q']==q);F=Field(r['p'],r['n'],r['modulus_low_first']);z=record(F,r['witness']);z.update(q=q,p=r['p'],n=r['n'],modulus=F.modulus)
 pathlib.Path(__file__).with_name(f'witness-{q}.json').write_text(json.dumps(z,indent=2)+'\n')
