"""Exact fixed-name search for simultaneous reversals of original modules."""
import sys,itertools,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'response'))
from verify import tournament,polynomial,hp,extend,need,strong

def modules(a):
 n=len(a);full=(1<<n)-1;out=[]
 for m in range(1,full+1):
  if m.bit_count()<2:continue
  if all(a[u]&m in (0,m) for u in range(n) if not m>>u&1):out.append(m)
 return out

def edge_mask(n,m):
 return sum(1<<b for b,(i,j) in enumerate(itertools.combinations(range(n),2)) if (m>>i&1) and (m>>j&1))
def classify(ms):
 if len(ms)==1:return 'single'
 if all(not(a&b) for a,b in itertools.combinations(ms,2)):return 'pairwise_disjoint'
 if all(not(a&b) or a&b in (a,b) for a,b in itertools.combinations(ms,2)):return 'nested_laminar'
 return 'overlapping'

def basis(n):
 half=[m for m in range(1<<n) if m.bit_count()<=2]
 return sorted(a|(c<<n) for a in half for c in half)

def toggle(a,ms):
 n=len(a);b=list(a)
 for m in ms:
  for i,j in itertools.combinations(range(n),2):
   if (m>>i&1) and (m>>j&1):b[i]^=1<<j;b[j]^=1<<i
 return b

def certificate(a,ms,full_queries=False):
 n=len(a);need(set(ms)<=set(modules(a)),'original sets are modules');b=toggle(a,ms);B=basis(n);p=polynomial(a);q=polynomial(b);d={s:p.get(s,0)-q.get(s,0) for s in B if p.get(s,0)!=q.get(s,0)}
 queries=[];first=None;uniform=[]
 assignments=itertools.product(range(1<<n),repeat=2) if full_queries else itertools.product([m for m in range(1<<n) if m.bit_count()<=2],repeat=2)
 for x,c in assignments:
  mask=x|(c<<n);fa=hp(extend(a,x,c));fb=hp(extend(b,x,c));need(sum(z for s,z in p.items() if s&mask==s)==fa,'literal source polynomial');need(sum(z for s,z in q.items() if s&mask==s)==fb,'literal target polynomial')
  if first is None and fa!=fb:first=[x,c,fa,fb]
  queries.append([x,c,fa,fb])
 if n==6 and ms==[3,12,48,63]:
  unis=[sum(3<<(2*i) for i in range(3) if m>>i&1) for m in range(8)]
  lookup={(x,c):(fa,fb) for x,c,fa,fb in queries}
  for x,c in itertools.product(unis,repeat=2):
   fa,fb=lookup[x,c];need((fa-fb)%4==0,'uniform block modular blindness');uniform.append([x,c,fa,fb])
  need(sum(a!=b for _,_,a,b in uniform)==24,'uniform exact differences');need(lookup[3,12]==(509,545),'uniform hostile');need((lookup[1,5][0]-lookup[1,5][1])%4==2,'nonuniform modular separator')
 return {'T':a,'U':b,'original_module_masks':ms,'family_type':classify(ms),'strong_T':strong(a),'strong_U':strong(b),'H_T':hp(a),'H_U':hp(b),'complete_named_basis':[[s&((1<<n)-1),s>>n] for s in B],'F_T_vector':[p.get(s,0) for s in B],'F_U_vector':[q.get(s,0) for s in B],'difference_vector':[p.get(s,0)-q.get(s,0) for s in B],'nonzero_difference_by_degree':dict(sorted(Counter(s.bit_count() for s in d).items())),'literal_query_count':len(queries),'literal_queries':queries,'first_separator':first,'pair_uniform_queries':uniform}

def main():
 summaries=[];bank={};total=0;allclasses=Counter();identities=0;collisions=[]
 for n in range(1,6):
  B=basis(n);count=1<<(n*(n-1)//2);vectors=[]
  for code in range(count):
   p=polynomial(tournament(n,code));vectors.append(tuple(p.get(s,0) for s in B))
  need(len(set(vectors))==count,'bounded full-response injectivity control')
  counts=Counter();two=Counter();source_records=[]
  for code in range(count):
   a=tournament(n,code);mods=modules(a);flips=[edge_mask(n,m) for m in mods];targets=[0]*(1<<len(mods));local=Counter()
   for f in range(1,1<<len(mods)):
    bit=f&-f;k=bit.bit_length()-1;targets[f]=targets[f^bit]^flips[k];target=code^targets[f];ms=[m for j,m in enumerate(mods) if f>>j&1];kind=classify(ms);counts[kind]+=1;local[kind]+=1;total+=1
    if len(ms)==2:two[kind]+=1
    if target==code:identities+=1
    elif vectors[code]==vectors[target]:collisions.append({'n':n,'source_code':code,'target_code':target,'module_masks':ms})
   source_records.append({'code':code,'modules':mods,'families':sum(local.values()),'family_types':dict(sorted(local.items()))})
  allclasses.update(counts)
  data={'n':n,'basis':[[s&((1<<n)-1),s>>n] for s in B],'records':[{'code':code,'rows':tournament(n,code),'vector':list(vec)} for code,vec in enumerate(vectors)]}
  raw=json.dumps(data,sort_keys=True,separators=(',',':')).encode();(Path(__file__).resolve().parent/f'full_response_n{n}.json').write_bytes(raw)
  summaries.append({'n':n,'tournaments':count,'complete_named_coordinate_count':len(B),'families_by_type':dict(sorted(counts.items())),'two_module_families_by_type':dict(sorted(two.items())),'source_records':source_records,'full_bank_sha256':hashlib.sha256(raw).hexdigest()})
 need(not collisions,'unexpected full-response collision');need(identities==0,'unexpected nonempty identity combination')
 a=[14,12,56,48,35,3]
 examples=[certificate(a,[3,12]),certificate(a,[3,63]),certificate(a,[3,12,48,63],True),certificate(tournament(5,1023),[7,28])]
 return {'status':'PASS','operation':'Choose distinct nontrivial modules of the original tournament, then reverse an edge iff it lies in an odd number of chosen sets. This is simultaneous edge XOR, not an assumption that each overlapping set remains a module after another flip.','tournaments':sum(s['tournaments'] for s in summaries),'nonempty_family_presentations':total,'families_by_type':dict(sorted(allclasses.items())),'nonempty_identity_families':identities,'exact_named_full_response_collisions':collisions,'scope':'Exhaustive only through five old names; larger examples are individually specified controls. This is a mechanism-specific census rather than a claimed new general reconstruction bound.','census':summaries,'examples':examples}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
