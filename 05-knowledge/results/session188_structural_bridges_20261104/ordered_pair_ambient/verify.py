"""Complete physical pair-child ambient and root-family controls."""
from pathlib import Path
from itertools import product
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def tt(n):return tuple(((1<<n)-1)^((1<<(u+1))-1) for u in range(n))
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):T[off[i]+u]=(row<<off[i])|sum(((1<<len(D))-1)<<off[j] for j,D in enumerate(blocks) if Q[i]>>j&1)
 return tuple(T),off
def endpoint(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for M in range(1,1<<n):
  for v in range(n):
   if M>>v&1:
    S=M^(1<<v);D[M][v]=1 if not S else sum(D[S][u] for u in range(n) if S>>u&1 and T[u]>>v&1)
 return D
def hp(T):return sum(endpoint(T)[-1])
def induced(T,M):
 vs=[u for u in range(len(T)) if M>>u&1];return tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def profile(T):
 n=len(T);H=list(map(sum,endpoint(T)));N=1<<n;D=[[0]*N for _ in range(n+1)];D[0][0]=1
 for k in range(1,n+1):
  for M in range(1,N):
   S=M
   while S:D[k][M]+=H[S]*D[k-1][M^S];S=(S-1)&M
 return [D[k][-1] for k in range(n+1)]
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  R[v if a>>v&1 else n]|=1<<(n if a>>v&1 else v);R[v if c>>v&1 else n+1]|=1<<(n+1 if c>>v&1 else v)
 return tuple(R)
def quartet(X,Y):return [(2*bool(s&1)|(X if s&2 else 0),int(bool(s&4))|(Y if s&8 else 0),(-1)**s.bit_count()) for s in range(16)]
def ambient(T,K):
 R,_=blow(K,[T]+[(0,)]*(len(K)-1));bits=sum(1<<(len(T)+j-1) for j in range(1,len(K)) if K[j]&1);return R,bits
def skew(T,X,Y):
 R=tuple(r>>2 for r in T[2:]);X>>=2;Y>>=2;M=X|Y;n=len(R);N=(1<<n)-1;need(module(R,M) and all(R[u]&Y==Y for u in range(n) if X>>u&1),'ordered module survives arbitrary ambient');E=endpoint(R);S=endpoint(op(R));gap=0
 for A in range(1,N):
  gap+=sum(E[A][u] for u in range(n) if Y>>u&1)*sum(S[N^A][v] for v in range(n) if X>>v&1)-sum(E[A][u] for u in range(n) if X>>u&1)*sum(S[N^A][v] for v in range(n) if Y>>v&1)
 x=profile(induced(R,X));y=profile(induced(R,Y));representatives=[(M&-M).bit_length()-1]+[v for v in range(n) if not M>>v&1];Q=tuple(sum(1<<j for j,v in enumerate(representatives) if R[u]>>v&1) for u in representatives);lower=hp(Q)*(x[2]*y[1]+x[1]*y[2]+x[2]*y[2]);need(gap>=lower>0,'global ordered-cover skew lower bound');return gap,lower

def main():
 allcover=[]
 for children in[(tt(2),tt(1),tt(2)),(tt(2),tt(2),tt(1))]:
  Q=(2,4,1);T,offs=blow(Q,children);U,_=blow(op(Q),[children[0],op(children[1]),op(children[2])]);W=flip(T,3);X=((1<<len(children[1]))-1)<<offs[1];Y=((1<<len(children[2]))-1)<<offs[2];corners=[]
  for a,c,s in quartet(X,Y):corners.append({'a':a,'c':c,'sign':s,'source':profile(extend(T,a,c)),'quotient':profile(extend(U,a,c)),'local':profile(extend(W,a,c))})
  lhs=[sum(r['sign']*(r['source'][h]-r['quotient'][h]) for r in corners) for h in range(len(T)+3)];rhs=[sum(r['sign']*(r['source'][h]-r['local'][h]) for r in corners) for h in range(len(T)+3)];need(lhs==rhs and lhs[1]>0,'every ordered-cover quartic restoration');allcover.append({'children':children,'source_rows':T,'quotient_rows':U,'local_rows':W,'all_cover_coefficients':lhs,'corners':corners})
 C=(2,4,1);C5=(6,12,24,17,3);W5=(10,28,25,16,1);systems=[('C3_crossing',C,[tt(2),tt(3),tt(2)],1,2),('cyclic5',C5,[tt(2)]*5,2,3),('W5',W5,[tt(2)]*5,1,2)];requests=[];ids={}
 def add(T,a,c):
  key=(tuple(T),a,c)
  if key not in ids:ids[key]=len(requests);requests.append(key)
  return ids[key]
 controls=[];families=[]
 for label,Q,children,xi,yi in systems:
  core,offs=blow(Q,children);target,_=blow(op(Q),children);X=((1<<len(children[xi]))-1)<<offs[xi];Y=((1<<len(children[yi]))-1)<<offs[yi];quart=quartet(X,Y)
  for K in((0,),C,C5):
   T,bits=ambient(core,K);U,_=ambient(target,K);L=flip(T,3);gap,lower=skew(T,X,Y);rec={'name':label,'outer':K,'source_rows':T,'quotient_rows':U,'local_rows':L,'X':X,'Y':Y,'outside_bits':bits,'skew':gap,'lower_bound':lower,'corners':[]}
   for a,c,s in quart:rec['corners'].append({'a':a|bits,'c':c|bits,'sign':s,'source_id':add(T,a|bits,c|bits),'quotient_id':add(U,a|bits,c|bits),'local_id':add(L,a|bits,c|bits)})
   controls.append(rec)
  if label=='W5':continue
  blocks=[((1<<len(B))-1)<<o for B,o in zip(children,offs)];mods=[M for M in range(1,1<<len(core)) if M.bit_count()>1 and M!=(1<<len(core))-1 and module(core,M)];compressed=[2]+blocks[1:];d=sum(bool(core[0]&B) for B in compressed);anchor=[(B|(1 if core[0]&B else 0),1) for B in compressed]
  if d%2:anchor.append((1,1))
  if(len(compressed)+d)%2:anchor.append((0,1))
  anchor.append((1,0));battery=list(dict.fromkeys([(a,c) for a,c,s in quart]+anchor));need(len(battery)<=len(Q)+18,'odd-prime battery bound');baseH=hp(core)
  for K in(C,C5):
   T,outside=ambient(core,K);baseline=[add(T,a|outside,c|outside) for a,c in battery];gap,_=skew(T,X,Y);bank={'name':label,'outer':K,'core_rows':core,'original_proper_modules':mods,'source_rows':T,'outside_bits':outside,'battery':battery,'anchor':anchor,'quartet':quart,'baseline_ids':baseline,'records':[],'expected_selected_pair_coefficient':gap}
   for choice in range(1<<len(mods)):
    U=op(core);selected=[]
    for j,M in enumerate(mods):
     if choice>>j&1:U=flip(U,M);selected.append(M)
    H=hp(U);need(H>=baseH,'hereditary scalar monotonicity');V,_=ambient(U,K);rec={'selected_modules':selected,'target_core_rows':U,'target_rows':V,'core_H':H,'source_core_H':baseH,'pair_selected':3 in selected,'value_ids':[add(V,a|outside,c|outside) for a,c in battery]};bank['records'].append(rec)
   families.append(bank)
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True,capture_output=True);graphs=[extend(T,a,c) for T,a,c in requests];wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';vals=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(vals)==len(graphs),'all actual physical counts')
 for rec in controls:
  q=l=0
  for r in rec['corners']:
   for field in('source','quotient','local'):r[field]=vals[r.pop(field+'_id')]
   q+=r['sign']*(r['source']-r['quotient']);l+=r['sign']*(r['source']-r['local'])
  need(q==l==rec['skew'],'ambient quotient/local equality and global skew');rec['actual_coefficient']=q
 branches={'strict':0,'pair_selected':0,'pair_unselected':0}
 for bank in families:
  base=bank['baseline_values']=[vals[t] for t in bank.pop('baseline_ids')];at={tuple(x):i for i,x in enumerate(bank['battery'])}
  for r in bank['records']:
   v=r['values']=[vals[t] for t in r.pop('value_ids')];need(v!=base,'every root family detected')
   if r['core_H']>r['source_core_H']:need(v[at[0,0]]>base[at[0,0]],'lifted internal-zero scalar strictness');branches['strict']+=1;r['branch']='strict'
   elif r['pair_selected']:
    z=sum(s*(base[at[a,c]]-v[at[a,c]]) for a,c,s in bank['quartet']);need(z==bank['expected_selected_pair_coefficient']>0,'ambient pair selected exact quartic');branches['pair_selected']+=1;r['branch']='pair_selected';r['coefficient']=z
   else:res=[(base[at[tuple(z)]]-v[at[tuple(z)]])%4 for z in bank['anchor']];need(2 in res,'ambient unselected anchor');branches['pair_unselected']+=1;r['branch']='pair_unselected';r['anchor_residues']=res
 pins=[ROOT/x for x in('PROOF.md','verify.py','count.cpp')]+[BASE/d/'PROOF.md' for d in('quartic_cover_transfer','ordered_module_skew','marked_pair_quotient','prime_pair_root_family','module_decomposition','critical_quotient_synthesis')]+[BASE/'returns_next'/'UNIFORM.md']
 print(json.dumps({'status':'PASS','all_cover_restoration_controls':allcover,'physical_ambient_quartics':controls,'root_family_contexts':families,'family_count':sum(len(b['records']) for b in families),'branches':branches,'distinct_actual_augmented_graphs':len(requests),'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins},'scope':'Specific ordered-module pair-child exact quartic and root-family theorem in every fixed ambient; not arbitrary quartic sign transfer.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
