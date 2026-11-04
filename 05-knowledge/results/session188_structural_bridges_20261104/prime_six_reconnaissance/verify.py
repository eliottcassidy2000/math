from pathlib import Path
from itertools import combinations,permutations,product
from collections import Counter
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;PAIRS=list(combinations(range(6),2))
def need(c,m):
 if not c:raise RuntimeError(m)
def rows(code):
 T=[0]*6
 for k,(i,j) in enumerate(PAIRS):
  u,v=(i,j) if code>>k&1 else(j,i);T[u]|=1<<v
 return T
def op(Q):return [((1<<len(Q))-1)^(1<<i)^v for i,v in enumerate(Q)]
def module(Q,M):return all(Q[i]&M in (0,M) for i in range(len(Q)) if not M>>i&1)
def literal_endpoints(Q):
 S=[0]*len(Q);E=S[:]
 for w in permutations(range(len(Q))):
  if all(Q[i]>>j&1 for i,j in zip(w,w[1:])):S[w[0]]+=1;E[w[-1]]+=1
 return S,E
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 T=[]
 for i,b in enumerate(B):T.extend((v<<off[i])|sum(((1<<len(z))-1)<<off[j] for j,z in enumerate(B) if Q[i]>>j&1) for v in b)
 return T,off
def expand(mask,B,off):return sum(((1<<len(b))-1)<<off[i] for i,b in enumerate(B) if mask>>i&1)
def extend(T,a,c):
 n=len(T);U=T[:]+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:U[i]|=1<<n
  else:U[n]|=1<<i
  if c>>i&1:U[i]|=1<<(n+1)
  else:U[n+1]|=1<<i
 return U
def ambient(T,a,c):
 n=len(T);T=[v|(1<<n) for v in T]+[1<<(n+1),(1<<n)-1];return extend(T,a|(1<<(n+1)),c|(1<<(n+1)))
def main():
 with tempfile.TemporaryDirectory(prefix='prime6-complete-') as td:
  ep=Path(td)/'end';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'endpoints.cpp'),'-o',str(ep)],check=True);raw=subprocess.check_output([str(ep)])
 need(raw==(ROOT/'endpoint_table.tsv').read_bytes(),'complete endpoint table byte equality');table=[list(map(int,s.split())) for s in raw.decode().splitlines()];need(len(table)==32768,'all codes')
 hist=Counter();phist=Counter();prime=[]
 for code,r in enumerate(table):
  need(len(r)==15 and r[0]==code,'exact code and table schema');Q=rows(code);wit=next((M for M in range(1,63) if M.bit_count()>=2 and module(Q,M)),0);need(r[1]==wit,'all proper subsets module completeness');S,E=r[3:9],r[9:15];need(all(0<=v<=120 for v in S+E) and sum(S)==sum(E) and sum(S)%2==1,'Hamilton accounting, endpoint bounds and oddness');p=sum(((S[i]+E[i])%2)<<i for i in range(6));need(p==r[2] and p>0,'nonzero endpoint parity for EVERY six-name tournament');hist[p.bit_count()]+=1
  if not wit:prime.append(code);phist[p.bit_count()]+=1
 need(dict(hist)=={2:13920,4:12480,6:6368} and len(prime)==10320 and dict(phist)=={6:3840,2:720,4:5760},'complete finite histograms')
 ps=list(permutations(range(6)));left=set(prime);orbits=[]
 while left:
  code=min(left);Q=rows(code);orb={sum(bool(Q[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(PAIRS)) for p in ps};need(orb<=left,'disjoint complete prime relabeling orbit');left-=orb;S,E=literal_endpoints(Q);need(S==table[code][3:9] and E==table[code][9:15],'representative direct permutation endpoints');orbits.append(dict(code=code,rows=Q,codes=sorted(orb),S=S,E=E,p=table[code][2]))
 need(len(orbits)==15 and sorted(len(r['codes']) for r in orbits)==[240]+[720]*14,'fifteen complete prime orbits')
 B8=[238,4,240,118,195,146,34,74];S8,E8=literal_endpoints(B8);need(sum(S8)==sum(E8) and sum(S8)%2==1 and all((x+y)%2==0 for x,y in zip(S8,E8)),'eight-name endpoint boundary');need(not any(module(B8,M) for M in range(1,255) if M.bit_count()>=2),'boundary is prime')
 jobs=[];controls=[]
 def enqueue(T,h):jobs.append((T,h));return len(jobs)-1
 for r in orbits:
  Q=r['rows'];i=(r['p']&-r['p']).bit_length()-1
  for name,B,outer in [('single',[tt(1)]*6,False),('pairs',[tt(2)]*6,False),('mixed',[tt(1),tt(2),[2,4,1],tt(1),tt(2),tt(1)],False),('mixed_ambient',[tt(1),tt(2),[2,4,1],tt(1),tt(2),tt(1)],True)]:
   ids=[];graphs=[];masks=[]
   for P in (Q,op(Q)):
    T,off=blow(P,B);a=expand(1<<i,B,off);G=ambient(T,a,0) if outer else extend(T,a,0);ids.append(enqueue(G,3));graphs.append(G);masks.append([a,0])
   controls.append(dict(code=r['code'],name=name,child_rows=B,anchor=i,outer_cyclic=outer,core_masks=masks,graphs=graphs,ids=ids))
 # Complete native physical response data, all cubic/quartic terms included.
 Q=orbits[0]['rows'];full=[]
 for a,c in product(range(64),repeat=2):full.append([a,c,enqueue(extend(Q,a,c),1),enqueue(extend(op(Q),a,c),1)])
 # Every root-plus-child-pair selection for this prime core, with fixed physical names.
 B=[tt(2)]*6;T,off=blow(Q,B);anchor=(orbits[0]['p']&-orbits[0]['p']).bit_length()-1;a=expand(1<<anchor,B,off);srcid=enqueue(extend(T,a,0),1);family=[]
 for bits in range(64):
  C=[tt(2) if bits>>i&1 else op(tt(2)) for i in range(6)];U,_=blow(op(Q),C);family.append([bits,enqueue(extend(U,a,0),1)])
 # Crossing original intervals inside a TT3 child; source, restored baseline and actual target.
 B=[tt(3)]+[tt(1)]*5;T,off=blow(Q,B);U,_=blow(op(Q),B);C=[[2,4,1]]+[tt(1)]*5;V,_=blow(op(Q),C);M=(1<<len(T))-1;a=expand(1<<anchor,B,off);cross=[]
 # Verify the actual target equals XOR reversal of root, {0,1}, {1,2} in the ORIGINAL graph.
 W=T[:]
 for i,j in combinations(range(len(T)),2):
  if (1+int(i in (0,1) and j in (0,1))+int(i in (1,2) and j in (1,2)))%2:W[i]^=1<<j;W[j]^=1<<i
 need(W==V and module(T,3) and module(T,6),'actual crossing original-module target')
 for aa,cc in [(a,0),(M,M^a)]:cross.append([aa,cc,[enqueue(ambient(X,aa,cc),3) for X in (T,U,V)]])
 with tempfile.TemporaryDirectory(prefix='prime6-physical-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);wire=str(len(jobs))+'\n'+''.join(f'{len(T)} {h}\n'+' '.join(map(str,T))+'\n' for T,h in jobs);values=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=wire.encode()).decode().splitlines()]
 need(len(values)==len(jobs),'all physical jobs consumed')
 for r in controls:
  r['profiles']=[values[j] for j in r.pop('ids')];need([(x-y)%4 for x,y in zip(*r['profiles'])]==[2,0,0],'complete actual physical profile separation')
 full=[[a,c,values[i][0],values[j][0]] for a,c,i,j in full];vals=[r[2:] for r in full];co=[v[:] for v in vals]
 for bit in range(12):
  for mask in range(4096):
   if mask>>bit&1:
    for side in (0,1):co[mask][side]-=co[mask^(1<<bit)][side]
 need(all(v==[0,0] for mask,v in enumerate(co) if (mask>>6).bit_count()>2 or (mask&63).bit_count()>2),'full named bidegree bound');coeff=[[a,c,*co[64*a+c]] for a in range(64) if a.bit_count()<=2 for c in range(64) if c.bit_count()<=2];need(len(coeff)==484,'all named coefficients')
 for r in family:r[1]=values[r[1]][0];need((values[srcid][0]-r[1])%4==2,'all64 root-plus-pair selections')
 for r in cross:r[2]=[values[j] for j in r[2]]
 basegaps=[r[2][0][0]-r[2][1][0] for r in cross];need(basegaps[0]==-basegaps[1] and basegaps[0]!=0,'fixed-ambient baseline duality');neg=next(i for i,z in enumerate(basegaps) if z<0);need(cross[neg][2][0][0]-cross[neg][2][2][0]<0,'negative query persists under crossing profile increase')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','endpoints.cpp','physical.cpp','endpoint_table.tsv']};deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in ['parity_profile_transfer/PROOF.md','module_decomposition/PROOF.md','monotone_dual_queries/PROOF.md']}
 return dict(status='PASS',source_pins=pins,dependency_pins=deps,complete_codes=32768,all_endpoint_support_histogram=dict(hist),prime_endpoint_support_histogram=dict(phist),prime_orbits=orbits,prime8_boundary=dict(rows=B8,S=S8,E=E8,p=0),actual_profile_controls=controls,complete4096_named_response_bank=full,complete484_named_coefficient_bank=coeff,root_pair_family_source_count=values[srcid][0],all64_original_root_pair_families=family,crossing_ambient_queries=cross,baseline_gaps=basegaps,sign_selected_query_index=neg,physical_graphs=len(jobs),scope='Finite exact six-name endpoint theorem, arbitrary-child and fixed-ambient consequences via inherited transfer, source-relative prime-root family consequence. No general even-order theorem or unrestricted reconstruction.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
