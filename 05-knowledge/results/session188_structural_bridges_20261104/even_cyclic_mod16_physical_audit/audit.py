"""Independent modulo-sixteen normalization and complete physical response audit."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,json,math,subprocess,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'even_cyclic_mod16'
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
def submasks(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M

def main():
 paths=[SRC/x for x in('PROOF.md','verify.py','literal.cpp')]+[BASE/d/'PROOF.md' for d in('five_profile_types','mod8_transfer','monotone_dual_queries','module_decomposition')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source receipts');formal=[[0,1,b,c,8,d] for b,c,d in product((6,14),(4,12),(0,8))];need(formal==data['formal_profile_upper_domain'],'complete eight-element upper domain');profiles=[]
 for r in data['actual_profile_controls']:
  T=tuple(r['rows']);need(len(T)>=4 and len(T)%2==0,'even nontrivial profile domain');q=profile(T);need(q==r['q'] and q[1]%2==1,'independent complete profile and odd unit');inv=pow(q[1],-1,16);norm=[0]+[(q[k] if k<len(q) else 0)*inv%16 for k in range(1,6)];need(norm==r['normalized_first_five'] and norm in formal and formal.index(norm)==r['type'],'exact normalization and type');need(q[4]%16==8 and all(x%16==0 for x in q[6:]),'four-slot parity and all tails');profiles.append([T,q,norm])
 need(len(profiles)==120,'all actual profile controls');S=[[0]*6 for _ in range(101)];S[0][0]=1
 for n in range(1,101):
  for k in range(1,6):S[n][k]=k*S[n-1][k]+S[n-1][k-1]
 need(all(S[n][4]%2==1 for n in range(4,101,2)),'finite controls of all-order Stirling parity proof')
 queries={};occurrences=0
 def save(T,a,c,value):
  nonlocal occurrences
  key=(tuple(T),a,c);need(key not in queries or queries[key]==value,'consistent actual physical value');queries[key]=value;occurrences+=1
 physical_records=[]
 for r in data['actual_physical_controls']:
  blocks=list(map(tuple,r['child_rows']));qs=[profile(B) for B in blocks];need(qs==r['profiles'],'actual full child profile');T,off=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks);need(T==tuple(r['source_old_rows']) and U==tuple(r['target_old_rows']),'fixed physical quotient-converse graphs');types=[]
  for q in qs:types.append(formal.index([0]+[(q[k] if k<len(q) else 0)*pow(q[1],-1,16)%16 for k in range(1,6)]))
  need(types==r['type_indices'],'normalized actual type indices')
  for z,(a,c),(v,w) in zip(r['block_queries'],r['physical_queries'],r['physical_values']):
   need([a,c]==[physical(z[0],blocks,off),physical(z[1],blocks,off)],'all actual whole-child bit expansions');save(T,a,c,v);save(U,a,c,w)
   if z in([1,3],[4,6]):need((v-w)%16==8,'actual universal query residue')
  physical_records.append({'child_orders':list(map(len,blocks)),'query_count':len(r['block_queries']),'types':types,'value_sha256':hashlib.sha256(json.dumps(r['physical_values'],separators=(',',':')).encode()).hexdigest()})
 bank=data['complete_named_TT4_bank'];T=tuple(bank['source_old_rows']);U=tuple(bank['target_old_rows']);base,off=blow((2,4,1),[tt(4)]*3);target,_=blow((4,1,2),[tt(4)]*3);need((T,U)==(base,target),'complete named physical base');half=[m for m in range(1<<12) if m.bit_count()<=2];expected=set(product(half,repeat=2));by={(r[0],r[1]):r for r in bank['records']};need(set(by)==expected and len(by)==6241,'entire bidegree22 basis including all overlaps and zeros');hist={k:0 for k in range(5)}
 for a,c in expected:
  r=by[a,c];save(T,a,c,r[2]);save(U,a,c,r[3]);coeff=[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*by[x,y][2+j] for x in submasks(a) for y in submasks(c)) for j in(0,1)];need(coeff==r[4:],'all complete named Boolean coefficients');hist[a.bit_count()+c.bit_count()]+=coeff[0]!=coeff[1]
 need({str(k):v for k,v in hist.items()}==bank['nonzero_difference_histogram'],'entire degree histogram')
 for a,c,v,w in bank['full_mask_reconstructions']:
  got=[sum(r[4+j] for(x,y),r in by.items() if x&a==x and y&c==y) for j in(0,1)];need(got==[v,w],'all full-mask reconstructions');save(T,a,c,v);save(U,a,c,w)
 blocks=[tt(4)]*3;pair=data['root_family_block_queries'];need(pair==[[1,3],[4,6]],'exact core converse-complement pair');pairphys=[(physical(a,blocks,off),physical(c,blocks,off)) for a,c in pair];baseline=data['root_family_baseline_values'];d=[a-b for a,b in baseline];need(d[0]==-d[1] and d[0]%16==8 and d[0]<0,'known common negative baseline');mods=[M for M in range(1,16) if M.bit_count()>1 and module(blocks[0],M)];need(len(mods)==6,'all original TT4 nontrivial interval modules');expected_choices={tuple(M for j,M in enumerate(mods) if bits>>j&1) for bits in range(64)};seen=set();strict=0
 for(a,c),(v,w) in zip(pairphys,baseline):save(T,a,c,v);save(U,a,c,w)
 for r in data['root_family_records']:
  selected=tuple(r['selected_original_first_child_modules']);need(selected in expected_choices and selected not in seen,'complete first-child interval-family controls');seen.add(selected);targetblocks=[op(flip(blocks[0],selected)),op(blocks[1]),op(blocks[2])];need(targetblocks==list(map(tuple,r['target_child_rows'])),'exact physical root-family target interiors');tq=[profile(B) for B in targetblocks];need(tq==r['target_profiles'],'complete target child profiles');old=profile(tt(4));need(all(y>=x for z in tq for x,y in zip(old,z)),'full hereditary profile monotonicity');isstrict=tq[0]!=old;need(isstrict==r['profile_strict'],'strict profile flag');strict+=isstrict;V,_=blow((4,1,2),targetblocks);need(r['negative_query_index']==0,'same negative query for every target')
  for j,((a,c),(v,w)) in enumerate(zip(pairphys,r['physical_values'])):save(T,a,c,v);save(V,a,c,w)
  need(r['physical_values'][0][0]-r['physical_values'][0][1]<=d[0]<0,'common strict monotone direction')
 need(seen==expected_choices and strict==40,'all64 first-child family controls and40 strict cases')
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'count.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(queries))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G))+f' {a} {c}' for G,a,c in queries)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(values==list(queries.values()),'every literal physical value independently recomputed')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','scope':'Independent domain normalization, complete actual profiles and physical named/family controls. The512x64 formal table is independently audited in a separate packet. Local-only modulo16 theorem here.','actual_complete_profiles':len(profiles),'profile_records_sha256':hashlib.sha256(json.dumps(profiles,separators=(',',':')).encode()).hexdigest(),'physical_controls':physical_records,'complete_named_basis_size':len(by),'complete_named_difference_histogram':hist,'full_mask_queries':len(bank['full_mask_reconstructions']),'root_family_controls':len(seen),'strict_families':strict,'independent_unique_actual_graph_queries':len(queries),'retained_physical_occurrences':occurrences,'source_sha256':pins,'audit_sha256':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in('audit.py','count.cpp')}},sort_keys=True,indent=2))
if __name__=='__main__':main()
