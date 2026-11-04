"""Independent physical-query audit of uniform aggregate-start tomography."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,math,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'uniform_block_starts'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def rows(n,code):
 T=[0]*n
 for j,(u,v) in enumerate(combinations(range(n),2)):
  a,b=(v,u) if code>>j&1 else (u,v);T[a]|=1<<b
 return tuple(T)
def endpoint(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)];P=op(T)
 for u in range(n):D[1<<u][u]=1
 for M in range(1,1<<n):
  rem=M
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1;S=M^bit
   if S:D[M][v]=sum(D[S][u] for u in range(n) if S>>u&1 and P[v]>>u&1)
 return D

def profile(T):
 n=len(T);out=[0]*(n+1)
 for word in permutations(range(n)):
  bad=sum(1<<j for j,(u,v) in enumerate(zip(word,word[1:])) if not T[u]>>v&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad==bad:out[cuts.bit_count()+1]+=1
 return out
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def induce(T,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
def main():
 paths=[SRC/'PROOF.md',SRC/'verify.py',BASE/'covered_quotient'/'PROOF.md',BASE/'module_decomposition'/'PROOF.md'];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());checks=[];queries={}
 def register(T,U,rec):
  I=rec['physical_mask'];s=sum(v for u,v in enumerate(endpoint(op(T))[-1]) if I>>u&1);t=sum(v for u,v in enumerate(endpoint(op(U))[-1]) if I>>u&1);need([s,t]==rec['start_counts'],'independent physical start counts')
  vals=rec['T_queries']+rec['U_queries'];need(vals[0]-vals[2]-2*(vals[1]-vals[3])==2*(s-t),'exact two-query tomography')
  for G,expected in ((T,vals[:2]),(U,vals[2:])):
   for (a,c),z in zip(((I,0),(0,I)),expected):
    key=(tuple(G),a,c);need(key not in queries or queries[key]==z,'consistent retained query counts');queries[key]=z
  checks.append([len(T),I,s-t,vals[0]-vals[2],vals[1]-vals[3]])
 for rec in data['all_labelled_converse_union_controls']:
  T=rows(rec['n'],rec['code']);register(T,op(T),rec)
 for rec in data['general_quotient_profile_preserving_controls']:
  Q=rec['quotient'];blocks=rec['child_rows'];other=rec['new_child_rows'];need([profile(A) for A in blocks]==[profile(A) for A in other],'independent child-profile preservation');T,off=blow(Q,blocks);U,_=blow(op(Q),other);mask=sum(((1<<len(A))-1)<<off[i] for i,A in enumerate(blocks) if rec['block_union']>>i&1);need(mask==rec['physical_mask'],'whole-child union physical mask');register(T,U,rec)
 for rec in data['arbitrary_third_controls']:
  a,b=rec['transitive_sizes'];blocks=[rows(a,0),rows(b,0),rec['third_rows']];T,_=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks);register(T,U,rec);need(any(rec['contrasts']),'literal separation');q=profile(blocks[-1]);need(q==rec['third_profile'],'arbitrary third explicit profile');need(abs(rec['signed_start_skew'])>=q[2]*abs(2**a-2**b)>0,'quantitative strict bound')
 blocks=[rows(2,0)]*3;T,_=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks);register(T,U,data['equal_size_boundary'])
 blocks=data['all_start_skews_zero_control']['child_rows'];T,_=blow((2,4,1),blocks);U,_=blow((4,1,2),blocks)
 for rec in data['all_start_skews_zero_control']['uniform_start_records']:register(T,U,rec)
 for G,z in zip((T,U),data['all_start_skews_zero_control']['other_values']):queries[(G,3,12)]=z
 rootcounts={'scalar_strict':0,'scalar_equal':0};blocks=[rows(2,0),rows(3,0),(2,4,1)];T,off=blow((2,4,1),blocks)
 for rec in data['root_family_controls']:
  U=list(T)
  for M in rec['selected_original_modules']:
   need(all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1),'original module verification')
   for u in range(len(T)):
    if M>>u&1:U[u]^=M^(1<<u)
  U=tuple(U);h=sum(endpoint(T)[-1]);g=sum(endpoint(U)[-1]);need([h,g]==rec['H'] and g>=h,'independent scalar branch')
  for G,z in ((T,h),(U,g)):queries[(G,0,0)]=z
  for G,vals in ((T,rec['T_queries']),(U,rec['U_queries'])):
   for (a,c),z in zip(((224,0),(0,224)),vals):queries[(G,a,c)]=z
  if h==g:
   rootcounts['scalar_equal']+=1;need(any(x!=y for x,y in zip(rec['T_queries'],rec['U_queries'])),'equal scalar separator')
   for i,A in enumerate(blocks):need(profile(A)==profile(induce(U,range(off[i],off[i]+len(A)))),'hereditary child-profile restoration')
  else:rootcounts['scalar_strict']+=1
 def q(n,k):return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1)) if 1<=k<=n else 0
 minors=0
 for a in range(1,20):
  for b in range(a+1,21):
   for k in range(1,22):need(q(b,k)*q(a,k-1)-q(b,k-1)*q(a,k)>=0,'independent inclusion-exclusion Stirling minor');minors+=1
   need(q(b,2)-q(a,2)==2**b-2**a,'strict second coordinate including order one')
 inputs=[str(len(G))+' '+' '.join(map(str,G))+f' {a} {c}' for G,a,c in queries]
 with tempfile.TemporaryDirectory(prefix='uniform-start-audit-') as tmp:
  exe=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True,capture_output=True);out=subprocess.check_output([str(exe)],input=(str(len(inputs))+'\n'+'\n'.join(inputs)+'\n').encode());need(list(map(int,out.split()))==list(queries.values()),'all independent physical query values')
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer mode replay');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','independent_start_query_controls':checks,'independent_distinct_actual_queries':len(queries),'inclusion_exclusion_minors':minors,'root_family_branches':rootcounts,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'literal.cpp')}},sort_keys=True,indent=2))
if __name__=='__main__':main()
