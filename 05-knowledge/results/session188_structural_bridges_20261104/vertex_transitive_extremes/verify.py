"""Vertex-transitive module/profile preservation and exact W root-family controls."""
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  a,b=(u,v) if code>>k&1 else(v,u);T[a]|=1<<b
 return tuple(T)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def induce(T,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def modules(T,whole=True):return [M for M in range(1,(1<<len(T)) if whole else (1<<len(T))-1) if M.bit_count()>1 and module(T,M)]
def flip(T,M):return tuple(row^(M^(1<<u) if M>>u&1 else 0) for u,row in enumerate(T))
def endpoint(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for u in range(n):D[1<<u][u]=1
 for M in range(1,1<<n):
  for u in range(n):
   if not M>>u&1:continue
   R=M^(1<<u)
   if R:D[M][u]=sum(D[R][v] for v in range(n) if R>>v&1 and T[v]>>u&1)
 return D
@lru_cache(None)
def marked(T):
 n=len(T);N=1<<n;E=endpoint(T);S=endpoint(op(T));H=[sum(x) for x in E];H[0]=1;Q=[[0]*(n+1) for _ in range(N)];Q[0][0]=1
 for mask in range(1,N):
  sub=mask
  while sub:
   h=H[sub]
   if h:
    for k in range(1,mask.bit_count()+1):Q[mask][k]+=h*Q[mask^sub][k-1]
   sub=(sub-1)&mask
 ends=[[0]*(n+1) for _ in T];starts=[[0]*(n+1) for _ in T]
 for mask in range(1,N):
  for u in range(n):
   if not mask>>u&1:continue
   for k in range(1,n+1):ends[u][k]+=E[mask][u]*Q[(N-1)^mask][k-1];starts[u][k]+=S[mask][u]*Q[(N-1)^mask][k-1]
 need(all(sum(x[k] for x in ends)==sum(x[k] for x in starts)==Q[-1][k] for k in range(n+1)),'marked normalization');return Q[-1],ends,starts
@lru_cache(None)
def automorphisms_to(T):
 n=len(T);degree=[r.bit_count() for r in T];witness=[]
 for target in range(n):
  f=[target]+[-1]*(n-1)
  def rec(u,used):
   if u==n:return tuple(f)
   for v in range(n):
    if used>>v&1 or degree[u]!=degree[v]:continue
    if any(bool(T[u]>>j&1)!=bool(T[v]>>f[j]&1) for j in range(u)):continue
    f[u]=v;out=rec(u+1,used|1<<v)
    if out is not None:return out
   return None
  if degree[0]!=degree[target]:return None
  fnd=rec(1,1<<target)
  if fnd is None:return None
  witness.append(fnd)
 return witness
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def extend(T,a,c):
 n=len(T);out=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:out[u]|=1<<n
  else:out[n]|=1<<u
  if c>>u&1:out[u]|=1<<(n+1)
  else:out[n+1]|=1<<u
 return tuple(out)
def W(k):return tuple(sum(1<<v for v in range(k+1) if(u and v and u<v) or(u==0 and v and (v-1)%2==0) or(v==0 and u and (u-1)%2==1)) for u in range(k+1))

def main():
 C3=(2,4,1);D9,_=blow(C3,[C3]*3);regular=[];vt=[];all_small=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);all_small+=1
   if len({r.bit_count() for r in T})!=1:continue
   ms=modules(T);need(all(not(M&N) or M&N in(M,N) for M,N in combinations(ms,2)),'regular module laminarity')
   for M in ms:
    A=induce(T,[u for u in range(n) if M>>u&1]);need(len({r.bit_count() for r in A})==1,'every regular ambient module regular')
   regular.append({'rows':T,'modules':ms});auts=automorphisms_to(T)
   if auts is not None:
    for M in ms:need(automorphisms_to(induce(T,[u for u in range(n) if M>>u&1])) is not None,'every vertex-transitive module vertex-transitive')
    vt.append({'rows':T,'automorphisms_from_zero':auts})
 need(automorphisms_to(D9) is not None,'explicit nested vertex-transitive child');dmods=modules(D9);need(dmods==[7,56,448,511],'complete nested original module list');baseline=marked(D9);need(baseline[1]==baseline[2],'nested child endpoint balance');profiles=[]
 for bits in range(1<<len(dmods)):
  U=D9;selected=[]
  for k,M in enumerate(dmods):
   if bits>>k&1:U=flip(U,M);selected.append(M)
  p=marked(U);need(p==baseline,'complete marked-profile preservation under all nested XOR families');profiles.append({'selected_original_modules':selected,'target_rows':U,'q':p[0],'E':p[1],'S':p[2]})
 ambient=[]
 for code in range(8):
  Q=rows(3,code);T,_=blow(Q,[C3,(2,0),(0,)]);U=flip(T,7);need(marked(T)==marked(U),'balanced module converse preserves all ambient marked profiles');ambient.append({'quotient':Q,'T':T,'U':U})
 graphlist=[];graphids={}
 def add(G):
  if G not in graphids:graphids[G]=len(graphlist);graphlist.append(G)
  return graphids[G]
 rootcases=[]
 for label,k,children in [('nested_extreme',2,[C3,D9,C3]),('crossing_middle',4,[C3,C3,(6,4,0),(0,),C3])]:
  Q=W(k);T,off=blow(Q,children);n=len(T);m=len(Q);blocks=[((1<<len(A))-1)<<off[t] for t,A in enumerate(children)];ms=modules(T,whole=False);need(all(any(M&B==M for B in blocks) for M in ms),'prime root proper-module containment');need(len(ms)==6,'six proper modules in retained family case');i,j=off[1],off[-1];quartet=[]
  for bits in range(16):
   a=(1 if bits&1 else 0)|((1<<i) if bits&2 else 0);c=(2 if bits&4 else 0)|((1<<j) if bits&8 else 0);quartet.append((a,c,(-1)**bits.bit_count()))
  compressed=[2,4]+blocks[1:];d=sum(bool(T[0]&B) for B in compressed);anchor=[(B|(1 if T[0]&B else 0),1) for B in compressed]
  if d%2:anchor.append((1,1))
  if(len(compressed)+d)%2:anchor.append((0,1))
  anchor.append((1,0));queries=list(dict.fromkeys([(a,c) for a,c,s in quartet]+anchor));need(len(queries)<=m+20,'fixed actual battery bound');records=[]
  for bits in range(1<<len(ms)):
   U=op(T);selected=[]
   for t,M in enumerate(ms):
    if bits>>t&1:U=flip(U,M);selected.append(M)
   for z in(1,len(children)-1):
    new=induce(U,range(off[z],off[z]+len(children[z])));need(marked(new)==marked(children[z]),'physical extreme marked-profile restoration')
   records.append({'selected_proper_modules':selected,'hub_selected':blocks[0] in selected,'target_rows':U,'value_ids':[add(extend(U,a,c)) for a,c in queries]})
  rootcases.append({'label':label,'quotient':Q,'children':children,'T':T,'proper_modules':ms,'queries':queries,'quartic':quartet,'anchor':anchor,'baseline_ids':[add(extend(T,a,c)) for a,c in queries],'families':records})
 with tempfile.TemporaryDirectory(prefix='vertex-transitive-extremes-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(graphlist))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphlist)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphlist),'complete literal count response')
 totals={'families':0,'scalar_strict':0,'equal_hub_selected':0,'equal_hub_unselected':0}
 for rec in rootcases:
  base=[values[t] for t in rec.pop('baseline_ids')];rec['baseline_values']=base;at={tuple(q):i for i,q in enumerate(rec['queries'])};h=base[at[(0,0)]];off=[sum(map(len,rec['children'][:i])) for i in range(len(rec['children']))]
  for f in rec['families']:
   z=[values[t] for t in f.pop('value_ids')];f['values']=z;g=z[at[(0,0)]];need(g>=h and z!=base,'scalar monotonicity and root-family separation');totals['families']+=1
   if g>h:f['branch']='scalar_strict';totals['scalar_strict']+=1
   else:
    for i,A in enumerate(rec['children']):need(marked(tuple(A))[0]==marked(induce(f['target_rows'],range(off[i],off[i]+len(A))))[0],'equality child q profiles')
    if f['hub_selected']:
     gap=sum(s*(base[at[(a,c)]]-z[at[(a,c)]]) for a,c,s in rec['quartic']);need(gap>0,'selected hub positive disjoint quartic');f['quartic_gap']=gap;f['branch']='equal_hub_selected';totals['equal_hub_selected']+=1
    else:
     residues=[(base[at[q]]-z[at[q]])%4 for q in rec['anchor']];need(2 in residues,'unselected hub actual anchor');f['anchor_residues']=residues;f['branch']='equal_hub_unselected';totals['equal_hub_unselected']+=1
 paths=[ROOT/n for n in('PROOF.md','verify.py','count.cpp')]+[ROOT.parent/d/'PROOF.md' for d in('cyclic_child_disjoint','module_decomposition','alternating_prime_quotient')]+[ROOT.parent/'returns_next'/'UNIFORM.md']
 print(json.dumps({'status':'PASS','small_tournaments_checked':all_small,'regular_module_controls':regular,'vertex_transitive_controls':vt,'nested_child':{'rows':D9,'automorphisms_from_zero':automorphisms_to(D9),'modules':dmods,'all_family_marked_profiles':profiles},'balanced_module_ambient_controls':ambient,'root_families':rootcases,'totals':totals,'distinct_literal_augmented_tournaments':len(graphlist),'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
