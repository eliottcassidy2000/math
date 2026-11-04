#!/usr/bin/env python3
"""Exact marked circulation and actual fixed-name response controls."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,tempfile
from positive import positive,symbolic_identity
from kernel import kernel
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 R=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else(j,i);R[u]|=1<<v
 return R
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def flip(T,M):return [r^(M^(1<<i) if M>>i&1 else 0) for i,r in enumerate(T)]
def blow(Q,B):
 off=[sum(map(len,B[:i])) for i in range(len(B))];T=[0]*sum(map(len,B))
 for i,A in enumerate(B):
  for v,r in enumerate(A):
   T[off[i]+v]=r<<off[i]
   for j,C in enumerate(B):
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
def covers(T):
 n=len(T);H=[0]*(1<<n)
 for p in permutations(range(n)):
  for k in range(1,n+1):
   if k>1 and not T[p[k-2]]>>p[k-1]&1:break
   # Every prefix occurs (n-k)! times; use only increasing unused suffix.
   if tuple(sorted(p[k:]))==p[k:]:H[sum(1<<v for v in p[:k])]+=1
 q=[[0]*(1<<n) for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,1<<n):
   S=M
   while S:q[h][M]+=H[S]*q[h-1][M^S];S=(S-1)&M
 return [q[h][-1] for h in range(n+1)]
def formal_table():
 def add(*P):
  D={}
  for p in P:
   for k,v in p.items():D[k]=D.get(k,0)+v
  return {k:v for k,v in D.items() if v}
 def scale(c,P):return {k:c*v for k,v in P.items() if c*v}
 def mul(A,B):
  D={}
  for (m,t),a in A.items():
   for (s,u),b in B.items():
    if m&s:continue
    k=(m|s,t+u);D[k]=D.get(k,0)+a*b
  return {k:v for k,v in D.items() if v}
 def M(T):
  G=[[{} for _ in T] for _ in T]
  for k in range(1,4):
   for p in permutations(range(3),k):
    z={(sum(1<<v for v in p),0):1}
    for u,v in zip(p,p[1:]):z=mul(z,{(0,0):int(bool(T[u]>>v&1)),(0,1):1})
    G[p[0]][p[-1]]=add(G[p[0]][p[-1]],z)
  return G
 out=[]
 for T,pairs,expect in [([6,4,0],[(0,1),(1,2),(0,2)],[[-2,0,0,1],[2,0,0,1],[0,-2,1,1]]),([4,1,2],[(0,2)],[[0,2,-1,1]])]:
  G=M(T);C=[add(*(G[u][v] for u in range(3))) for v in range(3)];R=[add(*G[u]) for u in range(3)];H=add(*C)
  for (p,q),ex in zip(pairs,expect):
   S=add(mul(C[p],C[q]),scale(-1,mul(R[p],R[q])));D=add(mul(C[p],R[q]),scale(-1,mul(C[q],R[p])));E=add(G[p][q],scale(-1,G[q][p]));HE=mul(H,E);tab=[{t:v for (m,t),v in P.items() if m==7} for P in(S,D,E,HE)]
   target=[({0:v} if v else{}) for v in ex]
   if p==0 and q==2:target[2][1]=2*ex[2]
   need(tab==target,'complete formal tau table');need(not mul(H,S) and not mul(H,D) and not mul(H,HE),'all omitted terms vanish')
   out.append(dict(rows=T,p=p,q=q,coefficients=tab))
 return out

def main():
 formal=formal_table();positive_identity=symbolic_identity();TT=[6,4,0];C=[4,1,2];one=[0];two=[2,0];graphs=[];ids={}
 def val(T,a,c):
  G=augment(T,a,c)
  if G not in ids:ids[G]=len(graphs);graphs.append(G)
  return ids[G]
 def quart(T,U,X,Y,p,q,ext=0):
  return [dict(a=ext|(1<<p if k&1 else 0)|(X if k&2 else 0),c=ext|(1<<q if k&4 else 0)|(Y if k&8 else 0),sign=(-1)**k.bit_count(),ids=[val(S,ext|(1<<p if k&1 else 0)|(X if k&2 else 0),ext|(1<<q if k&4 else 0)|(Y if k&8 else 0)) for S in(T,U)]) for k in range(16)]
 def circ(T,U,X,Y,ext=0):
  pats=[(1,2,1),(2,4,1),(1,4,-1),(2,0,-1),(0,2,-1),(0,0,1)];R=[]
  for a,c,s in pats:
   for u in range(2):
    for v in range(2):aa=ext|a|(X if u else 0);cc=ext|c|(Y if v else 0);R.append(dict(a=aa,c=cc,sign=s*(-1)**(u+v),ids=[val(S,aa,cc) for S in(T,U)]))
  need(len({(r['a'],r['c']) for r in R})==24,'distinct battery');return R
 universal=[]
 specs=[(rows(n,c),a,i,j) for n in(2,3) for c in range(1<<(n*(n-1)//2)) for a in range(1<<n) for i in range(n) for j in range(n) if i!=j]
 specs +=[(rows(4,c),a,i,j) for c in(0,1,13,27,42,63) for a in(0,5,10,15) for i,j in((0,1),(2,3))]
 for R,ao,i,j in specs:
  Q=[ao<<1]+[r<<1 | (0 if ao>>v&1 else 1) for v,r in enumerate(R)];T,_=blow(Q,[TT]+[one]*len(R));V,_=blow(Q,[C]+[one]*len(R));U=flip(T,7);W=flip(V,7);X=1<<(i+3);Y=1<<(j+3);terms=kernel(R,[v for v in range(len(R)) if ao>>v&1],[v for v in range(len(R)) if not ao>>v&1],i,j)
  universal.append(dict(R=R,Aout=ao,i=i,j=j,T=T,V=V,kernel_terms=terms,tt_quartics=[quart(T,U,X,Y,p,q) for p,q in((0,1),(1,2),(0,2))],cycle_quartic=quart(V,W,X,Y,0,2),circulation=circ(T,U,X,Y)))
 small=[rows(n,c) for n in range(1,4) for c in range(1<<(n*(n-1)//2))]
 ws=[([],0,0,X,Y) for X in small for Y in small]
 ws +=[(rows(n,c),a,m,two,one) for n in range(3) for c in range(1<<(n*(n-1)//2)) for a in range(1<<n) for m in range(1<<n)]
 ws +=[(rows(3,0),0,4,two,two)]
 witness=[]
 for E,ao,mo,X,Y in ws:
  Q=setup(E,ao,mo);B=[TT,X,Y]+[one]*len(E);T,off=blow(Q,B);U,_=blow(op(Q),B);V,_=blow(Q,[C]+B[1:]);XM=((1<<len(X))-1)<<off[1];YM=((1<<len(Y))-1)<<off[2];f=positive(E,ao,mo,covers(X),covers(Y));witness.append(dict(E=E,Aout=ao,Mout=mo,blocks=B,Q=Q,T=T,U=U,V=V,formula=f,local=circ(T,flip(T,7),XM,YM),quotient=circ(T,U,XM,YM),cycle=quart(V,flip(V,7),XM,YM,0,2)))
 # Literal strong ambient comparisons, with fixed nonzero exterior port bits.
 ambient=[];prime5=[sum(1<<((i+d)%5) for d in(1,2)) for i in range(5)]
 for E,ao,mo in[([],0,0),([2,0],1,1)]:
  Q=setup(E,ao,mo);B=[TT,two,one]+[one]*len(E);T,_=blow(Q,B);U,_=blow(op(Q),B);V,_=blow(Q,[C]+B[1:]);XM=24;YM=32
  for K in([2,4,1],prime5):
   t,_=blow(K,[T]+[one]*(len(K)-1));u,_=blow(K,[U]+[one]*(len(K)-1));v,_=blow(K,[V]+[one]*(len(K)-1));n=len(T);ext=sum(1<<(n+j-1) for j in range(1,len(K)) if K[j]&1);ambient.append(dict(outer=K,core_Q=Q,blocks=B,T=t,U=u,V=v,exterior_port_mask=ext,local=circ(t,flip(t,7),XM,YM,ext),quotient=circ(t,u,XM,YM,ext),cycle=quart(v,flip(v,7),XM,YM,0,2,ext)))
 # Entire named bidegree bank, retaining shared supports, all cubic/quartic terms and zeroes.
 Q=[2,4,1];B=[TT,two,two];T,_=blow(Q,B);U,_=blow(op(Q),B);n=len(T);masks=[s for s in range(1<<n) if s.bit_count()<=2];bank=[dict(a=a,c=c,ids=[val(S,a,c) for S in(T,U)]) for a in masks for c in masks];extra=[dict(a=(37*k+5)&127,c=(59*k+13)&127,ids=[val(S,(37*k+5)&127,(59*k+13)&127) for S in(T,U)]) for k in range(32)]
 with tempfile.TemporaryDirectory(prefix='tt3-circulation-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(g))+' '+' '.join(map(str,g)) for g in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'all literal graphs')
 def fill(qs):
  for q in qs:q['values']=[values[i] for i in q.pop('ids')]
  return sum(q['sign']*(q['values'][0]-q['values'][1]) for q in qs)
 for r in universal:
  ds=[fill(q) for q in r['tt_quartics']];c=fill(r['cycle_quartic']);s=fill(r['circulation']);need(ds[0]+ds[1]-ds[2]==c==s==sum(r['kernel_terms']),'universal exact circulation');r['tt_coefficients']=ds;r['cycle_coefficient']=c
 for r in witness:
  a,b,c=[fill(r[k]) for k in('local','quotient','cycle')];need(a==b==c==r['formula']['coefficient']>0,'strict exact witness');r['coefficient']=a
 for r in ambient:
  a,b,c=[fill(r[k]) for k in('local','quotient','cycle')];need(a==b==c>0,'strong ambient exact identity');r['coefficient']=a
 table={}
 for r in bank:r['values']=[values[i] for i in r.pop('ids')];table[r['a'],r['c']]=r['values']
 def subsets(s):
  t=s
  while True:
   yield t
   if not t:break
   t=(t-1)&s
 for r in bank:r['coefficients']=[sum((-1)**(r['a'].bit_count()+r['c'].bit_count()-a.bit_count()-c.bit_count())*table[a,c][j] for a in subsets(r['a']) for c in subsets(r['c'])) for j in range(2)]
 for r in extra:
  r['values']=[values[i] for i in r.pop('ids')];need(r['values']==[sum(z['coefficients'][j] for z in bank if z['a']&~r['a']==0 and z['c']&~r['c']==0) for j in range(2)],'complete named polynomial')
 counts={d:sum(z['coefficients'][0]!=z['coefficients'][1] for z in bank if z['a'].bit_count()+z['c'].bit_count()==d) for d in range(5)}
 paths=[ROOT/f for f in('PROOF.md','verify.py','positive.py','kernel.py','count.cpp')]+[BASE/'cyclic_child_disjoint'/'PROOF.md',BASE/'ordered_cyclic_module'/'PROOF.md',BASE/'aggregated_cyclic_module'/'PROOF.md',BASE/'quartic_cover_transfer'/'PROOF.md']
 return dict(status='PASS',formal_effective_path_table=formal,positive_identity=positive_identity,universal_controls=universal,positive_witness_controls=witness,strong_ambient_controls=ambient,complete_named_bank=dict(T=T,U=U,coordinates=bank,full_mask_controls=extra,nonzero_contrast_by_degree=counts),distinct_actual_graphs=len(graphs),source_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
