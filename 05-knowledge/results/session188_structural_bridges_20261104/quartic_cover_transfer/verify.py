#!/usr/bin/env python3
"""All-cover quartic transfer and complete fixed-ambient physical controls."""
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,tempfile
from functools import lru_cache
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j) if code>>k&1 else(j,i);T[a]|=1<<b
 return T
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def flip(T,M):return [r^(M^(1<<i) if M>>i&1 else 0) for i,r in enumerate(T)]
def module(T,M):return all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)
def covers(T):
 n=len(T);N=1<<n;d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for M in range(1,N):
  for v in range(n):
   if not M>>v&1:continue
   rest=M^(1<<v)
   if rest:d[M][v]=sum(d[rest][u] for u in range(n) if T[u]>>v&1)
 H=list(map(sum,d));q=[[0]*N for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,N):
   S=M
   while S:q[h][M]+=H[S]*q[h-1][M^S];S=(S-1)&M
 return [q[h][-1] for h in range(n+1)]
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for v,r in enumerate(B):
   T[off[i]+v]=r<<off[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+v]|=((1<<len(C))-1)<<off[j]
 return T,off
def augment(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for v in range(n):
  if a>>v&1:U[v]|=1<<n
  else:U[n]|=1<<v
  if c>>v&1:U[v]|=1<<(n+1)
  else:U[n+1]|=1<<v
 return tuple(U)
def setup(E,ao,mo):
 Q=[2,4,1]+[r<<3 for r in E]
 for j in range(len(E)):
  if ao>>j&1:Q[0]|=1<<(j+3)
  else:Q[j+3]|=1
  if mo>>j&1:Q[1]|=1<<(j+3);Q[2]|=1<<(j+3)
  else:Q[j+3]|=6
 return Q
def family(kind,r=2):
 n=2*r+1
 if kind=='T':return [sum(1<<((i+d)%n) for d in range(1,r+1)) for i in range(n)]
 if kind=='W':return [sum(1<<v for v in range(n) if(i and v and i<v)or(i==0 and v and v%2)or(v==0 and i and i%2==0)) for i in range(n)]
 Q=[0]*n
 for i,j in combinations(range(n),2):
  a,b=(j,i) if j<r else((i,j) if i>=r or j-r<=i else(j,i));Q[a]|=1<<b
 return Q
def battery(T,block_masks,z,u,v):
 Z=block_masks[z];p=(Z&-Z).bit_length()-1;q=next(t for t in range(len(T)) if Z>>t&1 and T[p]>>t&1);X=block_masks[u];Y=block_masks[v]
 quartet=[((1<<p if bits&1 else 0)|(X if bits&2 else 0),(1<<q if bits&4 else 0)|(Y if bits&8 else 0),(-1)**bits.bit_count()) for bits in range(16)]
 parts=[1<<t for t in range(len(T)) if Z>>t&1 and t!=p]+[M for i,M in enumerate(block_masks) if i!=z];d=sum(bool(T[p]&M) for M in parts);anchor=[(M|((1<<p) if T[p]&M else 0),1<<p) for M in parts]
 if d%2:anchor.append((1<<p,1<<p))
 if (len(parts)+d)%2:anchor.append((0,1<<p))
 anchor.append((1<<p,0));qs=list(dict.fromkeys([(a,c) for a,c,s in quartet]+anchor));need(len(qs)<=len(block_masks)+20,'m+20 root battery')
 return quartet,anchor,qs

@lru_cache(None)
def profile(T):return covers(T)
def induced(T,names):return [sum(1<<j for j,w in enumerate(names) if T[v]>>w&1) for v in names]
def endpoints(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for M in range(1,1<<n):
  for i in range(n):
   if M>>i&1 and M!=(1<<i):d[M][i]=sum(d[M^(1<<i)][j] for j in range(n) if T[j]>>i&1)
 return d
def marked(T,i):
 n=len(T);N=1<<n;end=endpoints(T);start=endpoints(op(T));qs=[profile(tuple(induced(T,[v for v in range(n) if M>>v&1]))) if M else [1] for M in range(N)];vals=[]
 for table in(end,start):
  p=[0]*(n+1)
  for h in range(1,n+1):
   for M in range(1,N):
    q=qs[(N-1)^M]
    if h-1<len(q):p[h]+=table[M][i]*q[h-1]
  vals.append(p)
 need(vals[0]==vals[1] and vals[0][1]>0,'complete physical marked balance');return vals[0]

def main():
 allcover=[]
 for n in(2,3,4):
  pairs=[(3,3)]+([(3,5)] if n>=3 else [])+([(3,12)] if n>=4 else [])
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);U=op(T)
   for am,cm in pairs:
    av=[i for i in range(n) if am>>i&1];cv=[i for i in range(n) if cm>>i&1];coeff=[[0]*(n+3) for _ in range(2)];controls=[]
    for bits in range(16):
     a=sum(1<<av[j] for j in range(2) if bits>>j&1);c=sum(1<<cv[j] for j in range(2) if bits>>(j+2)&1);ps=[profile(augment(V,a,c)) for V in(T,U)];controls.append(dict(a=a,c=c,source=ps[0],converse=ps[1]));sign=(-1)**bits.bit_count()
     for j in range(2):
      for h,z in enumerate(ps[j]):coeff[j][h]+=sign*z
    need(coeff[0]==coeff[1],'every ordered-cover quartic equals whole converse');allcover.append(dict(T=T,a_support=am,c_support=cm,quartic_profile=coeff[0],all16_profiles=controls))
 graphs=[];gid={}
 def value(T,a,c):
  G=augment(T,a,c)
  if G not in gid:gid[G]=len(graphs);graphs.append(G)
  return gid[G]
 C3=[2,4,1];TT2=[2,0];TT3=[6,4,0];W5=[14,28,8,16,5];prime5=[4,13,24,17,3];local=[]
 specs=[]
 for E,ao,mo in [([],0,0),([0],0,1),([2,0],1,2)]:
  for outer in([0],C3,prime5):specs.append((E,ao,mo,TT2,TT3,0,0,outer,True))
 for X,Y,mx,my in [(C3,C3,0,1),([0],W5,0,4),(W5,C3,4,2)]:
  for outer in([0],C3):specs.append(([0],1,0,X,Y,mx,my,outer,False))
 for E,ao,mo,X,Y,mx,my,outer,grouped in specs:
  Q=setup(E,ao,mo);blocks=[C3,X,Y]+[[0]]*len(E);T,off=blow(Q,blocks);U,_=blow(op(Q),blocks);W=flip(T,7);n=len(T);old=[blow(outer,[V]+[[0]]*(len(outer)-1))[0] for V in(T,U,W)];ext=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1)
  ax=((1<<len(X))-1)<<off[1] if grouped else 1<<(off[1]+mx);cy=((1<<len(Y))-1)<<off[2] if grouped else 1<<(off[2]+my);queries=[]
  for bits in range(16):
   a=ext|(1 if bits&1 else 0)|(ax if bits&2 else 0);c=ext|(2 if bits&4 else 0)|(cy if bits&8 else 0);queries.append(dict(a=a,c=c,sign=(-1)**bits.bit_count(),ids=[value(V,a,c) for V in old]))
  bx,by=(profile(tuple(X)),profile(tuple(Y))) if grouped else(marked(X,mx),marked(Y,my));exterior=[0]+list(range(3+len(X)+len(Y),len(old[0])));augE=induced(old[0],exterior);lower=2*bx[1]*by[1]*sum(endpoints(augE)[-1]);local.append(dict(grouped=grouped,E=E,Aout=ao,Mout=mo,X=X,Y=Y,marked_X=mx,marked_Y=my,marked_profiles=[bx,by],outer=outer,core_source=T,core_quotient=U,core_local=W,ambient_rows=old,strict_lower_bound=lower,queries=queries))
 # The master identity itself needs no C3, cyclic or positivity assumption.
 generic=[]
 for code in range(8):
  Q=rows(3,code)
  for A in([0],TT2,TT3,C3):
   T,off=blow(Q,[A,TT2,[0]]);U,_=blow(op(Q),[A,TT2,[0]]);W=flip(T,(1<<len(A))-1);p=0;q=0 if len(A)==1 else 1;ax=3<<off[1];cy=1<<off[2]
   for outer in(C3,prime5):
    n=len(T);ext=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);old=[blow(outer,[V]+[[0]]*(len(outer)-1))[0] for V in(T,U,W)];queries=[]
    for bits in range(16):
     a=ext|((1<<p) if bits&1 else 0)|(ax if bits&2 else 0);c=ext|((1<<q) if bits&4 else 0)|(cy if bits&8 else 0);queries.append(dict(a=a,c=c,sign=(-1)**bits.bit_count(),ids=[value(V,a,c) for V in old]))
    generic.append(dict(Q=Q,distinguished_child=A,outer=outer,ambient_rows=old,queries=queries))
 source=BASE/'aggregated_cyclic_module'/'normal.json';data=json.loads(source.read_text());roots=[]
 for base in data['root_families']:
  T=base['T'];B=base['blocks'];Q=base['Q'];n=len(T);masks=base['block_masks'];z,u,v=base['witness'];proper=[M for M in range(1,(1<<n)-1) if M.bit_count()>1 and module(T,M)];need(proper==base['proper_original_modules'],'all inherited original modules');outer=C3;ext=1<<(n+1);S,_=blow(outer,[T,[0],[0]]);qs=[(a|ext,c|ext) for a,c in base['battery']];rec=dict(label=base['label'],Q=Q,blocks=B,witness=[z,u,v],proper_modules=proper,core_source=T,ambient_source=S,outer=outer,internal_battery=base['battery'],lifted_battery=qs,quartic=base['quartic'],anchor=base['anchor'],baseline_ids=[value(S,a,c) for a,c in qs],families=[])
  for f in base['families']:
   U=op(T)
   for k,M in enumerate(proper):
    if f['selection']>>k&1:U=flip(U,M)
   need(U==f['target_rows'],'complete original-module XOR target');V,_=blow(outer,[U,[0],[0]]);rec['families'].append(dict(selection=f['selection'],core_target=U,ambient_target=V,core_scalar_increase=f['scalar_increase'],ids=[value(V,a,c) for a,c in qs]))
  roots.append(rec)
 # Entire orientation cube in a strongly connected C3 outer context.
 src_cube=data['whole_child_cube'];B=src_cube['blocks'];baseT,off=blow(C3,B);n=len(baseT);ext=1<<(n+1);cubeqs=[(a|ext,c|ext) for a,c in src_cube['queries']];cube=[]
 for s in src_cube['states']:
  bits=s['bits'];Q=op(C3) if bits&1 else C3;blocks=[op(b) if bits>>(i+1)&1 else b for i,b in enumerate(B)];T,_=blow(Q,blocks);need(T==s['T'],'independent cube state rows');S,_=blow(C3,[T,[0],[0]]);cube.append(dict(bits=bits,core_rows=T,ambient_rows=S,ids=[value(S,a,c) for a,c in cubeqs]))
 with tempfile.TemporaryDirectory(prefix='quartic-cover-transfer-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'complete literal output')
 for r in local:
  coeff=[0,0,0]
  for q in r['queries']:
   q['values']=[values[i] for i in q.pop('ids')]
   for j,v in enumerate(q['values']):coeff[j]+=q['sign']*v
  need(coeff[0]-coeff[1]==coeff[0]-coeff[2]>=r['strict_lower_bound']>0,'ambient quotient restoration equals positive local quartic');r['quartics']=coeff;r['positive_gap']=coeff[0]-coeff[1]
 for r in generic:
  coeff=[0,0,0]
  for q in r['queries']:
   q['values']=[values[i] for i in q.pop('ids')]
   for j,v in enumerate(q['values']):coeff[j]+=q['sign']*v
  need(coeff[1]==coeff[2],'general master transfer with arbitrary distinguished child');r['quartics']=coeff;r['contrast']=coeff[0]-coeff[1]
  if len(r['distinguished_child'])==1:need(coeff[0]==coeff[1],'singleton distinguished-child zero boundary')
 strict=equal=0
 for r in roots:
  base=r['baseline_values']=[values[i] for i in r.pop('baseline_ids')];at={tuple(q):i for i,q in enumerate(r['internal_battery'])};z=r['witness'][0];Z=((1<<len(r['blocks'][z]))-1)<<sum(map(len,r['blocks'][:z]));zi=r['proper_modules'].index(Z)
  for f in r['families']:
   vals=f['values']=[values[i] for i in f.pop('ids')];need(vals!=base,'every ambient root family full vector differs')
   if f['core_scalar_increase']:
    strict+=1;gap=vals[at[(0,0)]]-base[at[(0,0)]];need(gap>0,'lifted internal-zero strict branch');f['lifted_zero_gap']=gap
   else:
    equal+=1
    if f['selection']>>zi&1:
     gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in r['quartic']);need(gap>0,'restored hub ambient positive grouped quartic');f['quartic_gap']=gap
    else:
     residues=[(base[at[tuple(q)]]-vals[at[tuple(q)]])%4 for q in r['anchor']];need(2 in residues,'ambient exposed anchor');f['anchor_residues']=residues
 for s in cube:s['values']=[values[i] for i in s.pop('ids')]
 need(all(a['values']!=b['values'] for a,b in combinations(cube,2)),'complete16-state strong-ambient cube exact vectors')
 paths=[ROOT/n for n in('PROOF.md','verify.py','count.cpp')]+[BASE/'ordered_cyclic_module'/'PROOF.md',BASE/'aggregated_cyclic_module'/'PROOF.md',source,BASE/'module_decomposition'/'PROOF.md',BASE/'returns_next'/'UNIFORM.md',BASE/'returns_next'/'PROOF.md']
 return dict(status='PASS',cover_quartic_controls=allcover,general_master_controls=generic,ambient_local_quotient_controls=local,ambient_root_families=roots,root_family_count=strict+equal,strict_families=strict,equal_families=equal,strong_ambient_cube=dict(queries=cubeqs,states=cube),distinct_actual_augmented_graphs=len(graphs),source_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
