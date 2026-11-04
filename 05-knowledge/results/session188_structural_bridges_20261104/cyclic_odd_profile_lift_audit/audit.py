"""Independent complete canonical and physical ordered-cover lift audit."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,json,math,subprocess,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'cyclic_odd_profile_lift'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def tt(n):return tuple(((1<<n)-1)^((1<<(u+1))-1) for u in range(n))
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def blow(Q,blocks):
 offs=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):T[offs[i]+u]=(row<<offs[i])|sum(((1<<len(D))-1)<<offs[j] for j,D in enumerate(blocks) if Q[i]>>j&1)
 return tuple(T),offs
def physical(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def build(blocks,query,outer,target=None):
 T,offs=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks if target is None else target);n=len(T);a,c=[physical(s,blocks,offs) for s in query];outside=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);T,_=blow(outer,[T]+[(0,)]*(len(outer)-1));U,_=blow(outer,[U]+[(0,)]*(len(outer)-1));return T,U,(a|outside,c|outside)
def profile(T):
 n=len(T);N=1<<n;D=[[0]*n for _ in range(N)];H=[0]*N
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    R=M^(1<<v);D[M][v]=1 if not R else sum(D[R][w] for w in range(n) if R>>w&1 and T[w]>>v&1)
  H[M]=sum(D[M])
 @lru_cache(None)
 def part(M,k):
  if not M:return int(k==0)
  if k<=0:return 0
  b=M&-M;R=M^b;S=R;total=0
  while True:
   P=S|b;total+=H[P]*part(M^P,k-1)
   if not S:break
   S=(S-1)&R
  return total
 return [math.factorial(k)*part(N-1,k) for k in range(n+1)]
def module(T,M):return all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)
def flip(T,mods):
 U=list(T)
 for M in mods:
  for v in range(len(T)):
   if M>>v&1:U[v]^=M^(1<<v)
 return tuple(U)
def main():
 paths=[SRC/x for x in('PROOF.md','verify.py','literal.cpp')]+[BASE/d/'PROOF.md' for d in('five_profile_types','mod8_transfer','monotone_dual_queries')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source receipts');types={'P':tt(2),'E':tt(4),'O':tt(3),'C':(2,4,1)};expected=set(t for t in product(types,repeat=3) if sum(len(types[x]) for x in t)%2);need(len(expected)==32,'complete theoretical type domain');seen=set();bank={};queries={};occurrences=0
 def save(T,a,c,H,Q2):
  nonlocal occurrences
  key=(tuple(T),a,c);need(key not in queries or queries[key]==[H,Q2],'consistent named actual count');queries[key]=[H,Q2];occurrences+=1
 records=[]
 for r in data['canonical_q1_q2_records']:
  typ=tuple(r['types']);need(typ in expected and typ not in seen,'exact type enumeration');seen.add(typ);blocks=[types[t] for t in typ];need(list(map(list,blocks))==r['child_rows'],'actual canonical representatives');T,offs=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks);need(T==tuple(r['source_old_rows']) and U==tuple(r['target_old_rows']),'fixed named quotient orientations');good=[]
  for j,(a,c) in enumerate(product(range(8),repeat=2)):
   am,cm=physical(a,blocks,offs),physical(c,blocks,offs);need([am,cm]==r['physical_queries'][j],'every physical uniform mask');h1,h2,q1,q2=r['values'][j];save(T,am,cm,h1,q1);save(U,am,cm,h2,q2);need((h1-h2)%8 in(0,4),'nontrivial canonical residues')
   if(h1-h2)%8==4 and(q1-q2)%8==0:good.append(j)
  need(good==r['good_query_indices'] and good,'complete good-query index set');need(r['chosen_block_query']==list(divmod(good[0],8)),'deterministic first good query');bank[typ]=r;records.append({'types':list(typ),'good_indices':good,'complete_values_sha256':hashlib.sha256(json.dumps(r['values'],separators=(',',':')).encode()).hexdigest()})
 need(seen==expected and min(len(r['good_indices']) for r in records)==4 and max(len(r['good_indices']) for r in records)==20,'all32 classes and good-query bounds')
 for r in data['larger_child_ambient_controls']:
  blocks=list(map(tuple,r['child_rows']));qs=[profile(B) for B in blocks];need(qs==r['source_profiles'],'independent full child profiles');unit=1
  for typ,B,qv in zip(r['types'],blocks,qs):
   cq=profile(types[typ]);lam=qv[1]*pow(cq[1],-1,8)%8;unit=unit*lam%8;need(all((qv[k] if k<len(qv) else 0)%8==lam*(cq[k] if k<len(cq) else 0)%8 for k in range(max(len(qv),len(cq)))),'complete residue profile through tails')
  need(unit==r['unit'],'same odd homogeneous factor');T,U,(a,c)=build(blocks,r['block_query'],r['outer_rows']);need(T==tuple(r['source_old_rows']) and U==tuple(r['target_old_rows']) and[a,c]==r['physical_query'],'root-incidence lift exact physical names');hs,ht,qs,qt=r['values'];save(T,a,c,hs,qs);save(U,a,c,ht,qt);need((hs-ht)%8==4,'all ambient/larger separators')
  if len(r['outer_rows'])==1:need((qs-qt)%8==0 and all((x-unit*y)%8==0 for x,y in zip(r['values'],r['canonical_core_values'])),'both full-cover coordinates transfer by same unit')
 blocks=[tt(2),tt(2),tt(3)];safe=bank['P','P','O']['chosen_block_query'];dual=[7^safe[1],7^safe[0]]
 for r in data['unsafe_and_safe_odd_query_controls']:
  T,U,(a,c)=build(blocks,r['block_query'],r['outer_rows']);need([a,c]==r['physical_query'],'boundary masks');hs,ht,qs,qt=r['values'];save(T,a,c,hs,qs);save(U,a,c,ht,qt);need([hs-ht,qs-qt]==r['contrast'],'retained boundary exact contrast')
  if r['safe']:need((hs-ht)%8==4 and r['block_query']==safe,'safe boundary')
  elif len(r['outer_rows'])>1:need(hs-ht==-600 and(hs-ht)%8==0,'bad local residue may vanish in ambient')
 baseprofiles=[profile(B) for B in blocks];family_sets={};strict=0;negative_by_context={}
 for r in data['root_family_records']:
  target=[];signature=[]
  for B,mods in zip(blocks,r['selected_original_child_modules']):need(all(M.bit_count()>1 and module(B,M) for M in mods),'every selected module original');target.append(op(flip(B,mods)));signature.append(tuple(mods))
  need(target==list(map(tuple,r['target_child_rows'])),'exact physical root and child XOR');tq=[profile(B) for B in target];need(tq==r['target_profiles'] and all(b>=a for qa,qb in zip(baseprofiles,tq) for a,b in zip(qa,qb)),'complete profile monotonicity');key=tuple(r['outer_rows']);family_sets.setdefault(key,set()).add(tuple(signature));d=[]
  for j,z in enumerate((safe,dual)):
   T,U,(a,c)=build(blocks,z,key);base=r['baseline_values'][j];save(T,a,c,base[0],base[2]);save(U,a,c,base[1],base[3]);T,U,(a,c)=build(blocks,z,key,target);actual=r['target_values'][j];save(T,a,c,actual[0],actual[2]);save(U,a,c,actual[1],actual[3]);d.append(base[0]-base[1])
  need(d[0]==-d[1] and d[0]%8==4,'exact dual baseline');chosen=r['negative_query_index'];need(d[chosen]<0 and r['target_values'][chosen][0]-r['target_values'][chosen][1]<=d[chosen],'one common negative query');negative_by_context.setdefault(key,chosen);need(negative_by_context[key]==chosen,'same source query for every target');need(r['profile_strict']==(tq!=baseprofiles),'strict classification');strict+=r['profile_strict']
 need(len(data['root_family_records'])==64 and all(len(s)==32 for s in family_sets.values()) and strict==16,'all finite family/context controls')
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'count.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(queries))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in queries)+'\n';out=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need([out[i:i+2] for i in range(0,len(out),2)]==list(queries.values()),'every independent predecessor Hamilton/complementary-cover count')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','canonical_type_triples':32,'complete_queries_per_type':64,'canonical_records':records,'independent_unique_physical_graph_queries':len(queries),'physical_count_occurrences_with_repeated_family_sources':occurrences,'larger_child_ambient_controls':len(data['larger_child_ambient_controls']),'family_contexts':64,'strict_contexts':strict,'same_negative_query_by_context':[[list(k),v] for k,v in negative_by_context.items()],'source_sha256':pins,'audit_sha256':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in('audit.py','count.cpp')},'scope':'All-cover homogeneous residue transfer, complete odd-type certificates and fixed-ambient dual family theorem; no exact collision or unrestricted prime reconstruction.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
