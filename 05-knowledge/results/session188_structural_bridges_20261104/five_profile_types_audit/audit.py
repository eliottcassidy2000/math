"""Independent full run-word and physical predecessor count audit."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import hashlib,json,math,subprocess,tempfile
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'five_profile_types'
def need(c,msg):
 if not c:raise RuntimeError(msg)
def profile(T):
 n=len(T);N=1<<n;ends=[[0]*n for _ in range(N)];H=[0]*N
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    R=M^(1<<v);ends[M][v]=1 if not R else sum(ends[R][w] for w in range(n) if R>>w&1 and T[w]>>v&1)
  H[M]=sum(ends[M])
 # Unordered set partitions: force the least remaining name into the first part.
 @lru_cache(None)
 def part(M,k):
  if not M:return int(k==0)
  if k<=0:return 0
  b=M&-M;R=M^b;S=R;out=0
  while True:
   B=S|b;out+=H[B]*part(M^B,k-1)
   if not S:break
   S=(S-1)&R
  return out
 return [math.factorial(k)*part(N-1,k) for k in range(n+1)]
def cycles(T):
 n=len(T)
 if n==1:return 0
 @lru_cache(None)
 def paths(M,v):
  if M==1<<v:return int(v==0)
  R=M^(1<<v)
  return sum(paths(R,w) for w in range(n) if R>>w&1 and T[w]>>v&1)
 return sum(paths((1<<n)-1,v) for v in range(1,n) if T[v]&1)
def skeleton(a,c,reverse):
 def arc(u,v):
  if u==v:return False
  if u<3 and v<3:return ((u-v if reverse else v-u)%3)==1
  if u==3:return v==4 or(v<3 and not(a>>v&1))
  if u==4:return v<3 and not(c>>v&1)
  return bool((a if v==3 else c)>>u&1)
 @lru_cache(None)
 def f(counts,mask,last):
  total=sum(counts)+mask.bit_count()
  if last<3:
   if not counts[last]:return 0
   cc=list(counts);cc[last]-=1;cc=tuple(cc);mm=mask
  else:
   b=1<<(last-3)
   if not mask&b:return 0
   cc=counts;mm=mask^b
  if total==1:return 1
  return sum(f(cc,mm,v) for v in range(5) if arc(v,last))
 return {counts:sum(f(counts,3,last) for last in range(5)) for counts in product(range(1,5),repeat=3)}
def blow(blocks,reverse=False):
 offsets=[];T=[]
 for B in blocks:offsets.append(len(T));T.extend([0]*len(B))
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):
   T[offsets[i]+u]|=row<<offsets[i]
   for j,C in enumerate(blocks):
    if i!=j and ((i-j if reverse else j-i)%3)==1:T[offsets[i]+u]|=((1<<len(C))-1)<<offsets[j]
 return T,offsets
def physical(mask,blocks,offs):return sum(((1<<len(B))-1)<<off for i,(B,off) in enumerate(zip(blocks,offs)) if mask>>i&1)
def subsets(a):
 s=a
 while True:
  yield s
  if not s:break
  s=(s-1)&a

def main():
 source=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source modes')
 can=source['canonical_child_rows'];cp=[profile(T) for T in can];need(cp==source['canonical_profiles'],'exact canonical profiles');norm=[[1,0,0],[1,2,0],[1,6,4],[1,6,6],[1,2,2]];typ=['S','P','E','O','C'];controls=[]
 for r in source['actual_profile_controls']:
  T=r['rows'];q=profile(T);z=cycles(T) if len(T)%2 else 0;t=0 if len(T)==1 else 1 if len(T)==2 else 2 if len(T)%2==0 else 4 if z%2 else 3;unit=q[1]*pow(cp[t][1],-1,8)%8
  need(q==r['q'] and z==r['odd_Hamilton_cycles'] and typ[t]==r['type'] and unit==r['unit'],'independent complete actual profile and cycle controls')
  need([(q[k] if k<len(q) else 0)*pow(q[1],-1,8)%8 for k in range(1,4)]==norm[t] and all(v%8==0 for v in q[4:]),'five normalized types, all tails')
  controls.append([T,q,z,t,unit])
 kernels={(a,c,rev):skeleton(a,c,rev) for a,c,rev in product(range(8),range(8),[False,True])};coarse=[];blind=[]
 for r in source['canonical_records']:
  blocks=r['child_rows'];T,offs=blow(blocks);U,_=blow(blocks,True);need(T==r['source_old_rows'] and U==r['target_old_rows'],'physical canonical construction');qs=[profile(B)+[0]*4 for B in blocks];vv=[]
  for a,c in product(range(8),repeat=2):
   vals=[sum(v*math.prod(qs[i][e[i]] for i in range(3)) for e,v in kernels[a,c,rev].items()) for rev in [False,True]];vv.append(vals)
   need([physical(a,blocks,offs),physical(c,blocks,offs)]==r['physical_queries'][8*a+c],'all physical block masks restored')
  need(vv==r['values'],'every raw canonical integer count via independent run words');res=[(u-v)%8 for u,v in vv];need(res==r['contrasts_mod8'],'complete canonical residues')
  if not any(res):blind.append(r['types'])
  if 'S' not in r['types']:need(set(res)<={0,4},'all nontrivial types residue-four')
  coarse.append({'types':r['types'],'values_sha256':hashlib.sha256(json.dumps(vv,separators=(',',':')).encode()).hexdigest(),'nonzero':sum(bool(v) for v in res)})
 need(blind==[['E']*3],'exact unique blindness')
 # Independent raw physical DP for all 2,116 named basis values plus full-mask controls.
 bank=source['complete_named_TT3_bank'];queries=[r[:2] for r in bank['records']]+[r[:2] for r in bank['full_mask_reconstructions']];cases=[(bank['source_old_rows'],queries),(bank['target_old_rows'],queries)];expected=[]
 for side in [0,1]:expected.extend([r[2+side] for r in bank['records']]+[r[2+side] for r in bank['full_mask_reconstructions']])
 for r in source['actual_profile_transfer_controls']:
  blocks=r['child_rows'];T,offs=blow(blocks);U,_=blow(blocks,True);qq=[(physical(a,blocks,offs),physical(c,blocks,offs)) for a,c in product(range(8),repeat=2)]
  for side,rows in enumerate([T,U]):cases.append((rows,qq));expected.extend(v[side] for v in r['physical_values'])
 for r in source['even_total_ambient_controls']:
  for side,key in enumerate(['source_old_rows','target_old_rows']):cases.append((r[key],[r['physical_query']]));expected.append(r['values'][side])
 lines=[]
 for rows,qq in cases:lines.extend([f'{len(rows)} {len(qq)}',' '.join(map(str,rows))]+[f'{a} {c}' for a,c in qq])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True);got=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(got==expected,'all independent raw physical counts')
 dic={(r[0],r[1]):r for r in bank['records']};hist={d:0 for d in range(5)}
 for r in bank['records']:
  a,c=r[:2];cc=[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*dic[x,y][2+s] for x in subsets(a) for y in subsets(c)) for s in [0,1]];need(cc==r[4:],'every named cubic/quartic coefficient including shared names');hist[a.bit_count()+c.bit_count()]+=cc[0]!=cc[1]
 for a,c,u,v in bank['full_mask_reconstructions']:
  z=[sum(r[4+s] for r in bank['records'] if r[0]&a==r[0] and r[1]&c==r[1]) for s in [0,1]];need(z==[u,v],'full response reconstruction')
 need({str(k):v for k,v in hist.items()}==bank['nonzero_difference_histogram'],'complete degree histogram')
 pins={p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ['PROOF.md','verify.py','literal.cpp','normal.json','optimized.json']}
 print(json.dumps({'status':'PASS','source_pins':pins,'independent_profile_cycle_controls':len(controls),'profile_records_sha256':hashlib.sha256(json.dumps(controls,separators=(',',':')).encode()).hexdigest(),'coarse_integer_counts':125*64*2,'coarse_records':coarse,'only_blind':blind,'independent_physical_counts':len(got),'full_named_basis_size':len(bank['records']),'nonzero_coefficient_difference_histogram':hist,'literal_sha256':hashlib.sha256((HERE/'literal.cpp').read_bytes()).hexdigest(),'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Five-type block-uniform reduction and complete canonical C3 certificate; ambient theorem only even total order in frozen source.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
