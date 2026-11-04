"""Ordered exterior module: exact identity, physical profiles, and critical root families."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,tempfile
from algebra import symbolic_identity,positive
from kernel import kernel
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,c):
 T=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  a,b=(u,v) if c>>k&1 else(v,u);T[a]|=1<<b
 return tuple(T)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def extend(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:U[u]|=1<<n
  else:U[n]|=1<<u
  if c>>u&1:U[u]|=1<<(n+1)
  else:U[n+1]|=1<<u
 return tuple(U)
def marked(T,u=0):
 n=len(T);q=[0]*(n+1);E=[0]*(n+1);S=[0]*(n+1)
 for word in permutations(range(n)):
  bad=sum(1<<j for j,(a,b) in enumerate(zip(word,word[1:])) if not T[a]>>b&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   h=1+cuts.bit_count();last=(cuts&-cuts).bit_length()-1 if cuts else n-1;q[h]+=1;E[h]+=word[last]==u;S[h]+=word[0]==u
 need(E==S and E[1]>0,'physical endpoint balance');return q,E
def setup(E,Aout,Mout):
 # names: triangle, ordered side X, ordered side Y, then E.
 n=len(E);Q=[2,4,1]+[r<<3 for r in E]
 for j in range(n):
  if Aout>>j&1:Q[0]|=1<<(j+3)
  else:Q[j+3]|=1
  if Mout>>j&1:Q[1]|=1<<(j+3);Q[2]|=1<<(j+3)
  else:Q[j+3]|=6
 return tuple(Q)
def Tfamily(r):return tuple(sum(1<<((u+d)%(2*r+1)) for d in range(1,r+1)) for u in range(2*r+1))
def Ufamily(r):
 Q=[0]*(2*r+1)
 for i in range(r):
  for j in range(i):Q[i]|=1<<j
  for j in range(r+1):
   if j<=i:Q[i]|=1<<(r+j)
   else:Q[r+j]|=1<<i
 for j in range(r+1):
  for k in range(j+1,r+1):Q[r+j]|=1<<(r+k)
 return tuple(Q)
def Wfamily(r):
 k=2*r;return tuple(sum(1<<v for v in range(k+1) if(u and v and u<v) or(u==0 and v and (v-1)%2==0) or(v==0 and u and (u-1)%2==1)) for u in range(k+1))
def witness(Q,z,u,v):
 need(Q[z]>>u&1 and Q[u]>>v&1 and Q[v]>>z&1,'cyclic directions');need(all(bool(Q[t]>>u&1)==bool(Q[t]>>v&1) for t in range(len(Q)) if t not in(z,u,v)),'homogeneous deleted-root pair')
def main():
 identity=symbolic_identity();C3=(2,4,1);C5=Tfamily(2);one=(0,);profiles={A:marked(A)[1] for A in(one,C3,C5)};small=[];graphs=[];gids={};literal=[];families=[]
 def add(G):
  if G not in gids:gids[G]=len(graphs);graphs.append(G)
  return gids[G]
 def quart(T,U,p,q,i,j):
  rec=[]
  for bits in range(16):
   a=((1<<p) if bits&1 else 0)|((1<<i) if bits&2 else 0);c=((1<<q) if bits&4 else 0)|((1<<j) if bits&8 else 0);rec.append({'a':a,'c':c,'sign':(-1)**bits.bit_count(),'source_id':add(extend(T,a,c)),'target_id':add(extend(U,a,c))})
  return rec
 for n in range(4):
  for code in range(1<<(n*(n-1)//2)):
   E=rows(n,code)
   for ao in range(1<<n):
    for mo in range(1<<n):
     Q=setup(E,ao,mo);R=tuple(r>>1 for r in Q[1:]);out=[u for u in range(n+2) if Q[0]>>(u+1)&1];inn=[u for u in range(n+2) if not Q[0]>>(u+1)&1];f=positive(E,ao,mo,profiles[one],profiles[one]);terms=kernel(R,out,inn,0,1);need(sum(terms)==f['coefficient'],'all small arbitrary exterior physical kernel');small.append({'E':E,'Aout':ao,'Mout':mo,'coefficient':sum(terms),'lower_bound':f['lower_bound']})
 specs=[((),0,0,X,Y) for X in(one,C3,C5) for Y in(one,C3,C5)]
 for n,code,ao,mo in[(1,0,0,0),(1,0,0,1),(1,0,1,0),(1,0,1,1),(2,0,0,1),(2,1,2,3),(3,2,3,5),(3,5,6,2),(4,21,9,6),(4,53,3,12)]:
  specs.append((rows(n,code),ao,mo,C3,C3))
 specs +=[(rows(2,1),1,2,C5,C3),(rows(2,0),2,1,C3,C5)]
 for E,ao,mo,X,Y in specs:
  Q=setup(E,ao,mo);blocks=[C3,X,Y]+[one]*len(E);T,off=blow(Q,blocks);Ulocal=flip(T,7);Uroot,_=blow(op(Q),blocks);R=tuple(r>>3 for r in T[3:]);out=[u for u in range(len(R)) if T[0]>>(u+3)&1];inn=[u for u in range(len(R)) if not T[0]>>(u+3)&1];f=positive(E,ao,mo,profiles[X],profiles[Y]);terms=kernel(R,out,inn,0,len(X));need(sum(terms)==f['coefficient'],'complete profile functional equals physical kernel');literal.append({'E':E,'Aout':ao,'Mout':mo,'X':X,'Y':Y,'Q':Q,'T':T,'local_target':Ulocal,'quotient_target':Uroot,'physical_marks':[0,off[1],1,off[2]],'kernel_terms':terms,'positive_formula':f,'local_queries':quart(T,Ulocal,0,1,off[1],off[2]),'quotient_queries':quart(T,Uroot,0,1,off[1],off[2])})
 structures=[]
 for r in range(1,16):
  for family,fun in [('T',Tfamily),('U',Ufamily),('W',Wfamily)]:
   Q=fun(r);ws=[(0,r,r+1)] if family=='T' else([(i,r+i,r+i+1) for i in range(r)] if family=='U' else[(0,2*s+1,2*s+2) for s in range(r)])
   for z,u,v in ws:witness(Q,z,u,v)
   structures.append({'family':family,'r':r,'order':len(Q),'witnesses':ws})
 for label,Q,z,u,v in [('C3',C3,0,1,2),('T5',Tfamily(2),0,2,3),('U5',Ufamily(2),0,2,3),('W5',Wfamily(2),0,1,2)]:
  witness(Q,z,u,v);m=len(Q);T,off=blow(Q,[C3]*m);n=len(T);blocks=[7<<(3*t) for t in range(m)];mods=[M for M in range(1,(1<<n)-1) if M.bit_count()>1 and module(T,M)];need(mods==blocks,'exact all-C3 prime original module list');p=off[z];q=p+1;i=off[u];j=off[v];quartet=[(a['a'],a['c'],a['sign']) for a in quart(T,T,p,q,i,j)];compressed=[1<<(p+1),1<<(p+2)]+[blocks[t] for t in range(m) if t!=z];d=sum(bool(T[p]&B) for B in compressed);anchor=[(B|((1<<p) if T[p]&B else 0),1<<p) for B in compressed]
  if d%2:anchor.append((1<<p,1<<p))
  if(len(compressed)+d)%2:anchor.append((0,1<<p))
  anchor.append((1<<p,0));battery=list(dict.fromkeys([(a,c) for a,c,s in quartet]+anchor));need(len(battery)<=m+20,'safe fixed root-family query bound');rec={'label':label,'Q':Q,'witness':[z,u,v],'T':T,'quartic':quartet,'anchor':anchor,'battery':battery,'baseline_ids':[add(extend(T,a,c)) for a,c in battery],'families':[]}
  for bits in range(1<<m):
   U=op(T)
   for t,M in enumerate(blocks):
    if bits>>t&1:U=flip(U,M)
   rec['families'].append({'selected_children_bits':bits,'target_rows':U,'value_ids':[add(extend(U,a,c)) for a,c in battery]})
  families.append(rec)
 with tempfile.TemporaryDirectory(prefix='ordered-cyclic-module-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'all actual values')
 for rec in literal:
  for key in('local_queries','quotient_queries'):
   delta=0
   for q in rec[key]:q['source']=values[q.pop('source_id')];q['target']=values[q.pop('target_id')];delta+=q['sign']*(q['source']-q['target'])
   need(delta==rec['positive_formula']['coefficient']>0,'actual local/quotient positive coefficient')
 for rec in families:
  base=rec['baseline_values']=[values[t] for t in rec.pop('baseline_ids')];at={q:i for i,q in enumerate(rec['battery'])};z=rec['witness'][0]
  for f in rec['families']:
   vals=f['values']=[values[t] for t in f.pop('value_ids')];need(vals!=base and vals[at[(0,0)]]==base[at[(0,0)]],'all-C3 scalar-equal root-family detection')
   if f['selected_children_bits']>>z&1:
    gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in rec['quartic']);need(gap>0,'hub-restored quartic');f['quartic_gap']=gap
   else:
    f['anchor_residues']=[(base[at[q]]-vals[at[q]])%4 for q in rec['anchor']];need(2 in f['anchor_residues'],'hub-conversed anchor')
 paths=[ROOT/n for n in('PROOF.md','verify.py','algebra.py','kernel.py','count.cpp')]+[BASE/d/'PROOF.md' for d in('cyclic_child_disjoint','vertex_transitive_extremes','module_decomposition','critical_quotient_synthesis')]+[BASE/'returns_next'/'UNIFORM.md']
 print(json.dumps({'status':'PASS','symbolic_identity':identity,'small_arbitrary_exterior_controls':small,'balanced_profiles':[[list(T),b] for T,b in profiles.items()],'literal_positive_controls':literal,'critical_witness_structures':structures,'root_families':families,'distinct_literal_augmented_tournaments':len(graphs),'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
