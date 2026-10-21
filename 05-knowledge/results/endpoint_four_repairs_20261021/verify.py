#!/usr/bin/env python3
"""Portable exact controls for the endpoint-four movable-cut packet.
No inherited imports, external packages, asserts, or tournament census above n9.
"""
from collections import Counter, deque
from itertools import combinations
from pathlib import Path
import hashlib,json,random

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def edge(u,v):return (min(u,v),max(u,v))
def edges(q):return {edge(u,v) for u,v in zip(q,q[1:])}
def blue(T,q):
 r=dict(zip(q,q[::-1]));return all((r[v],r[u]) in T for u,v in T)
def directed(T,q):return len(set(q))==len(q) and all((u,v) in T for u,v in zip(q,q[1:]))
def move(T,q):
 E=edges(q);return {(u,v) if edge(u,v) in E else (v,u) for u,v in T}
def source(n,I,rng,mask=None):
 T=set();seen=set();index=0
 for u,v in combinations(range(1,n-1),2):
  if (u,v) in seen:continue
  orb={(u,v),(n-1-v,n-1-u)};seen|=orb
  bit=True if v==u+1 else bool((mask>>index)&1) if mask is not None else bool(rng.getrandbits(1))
  if v!=u+1:index+=1
  T.update((a,b) if bit else (b,a) for a,b in orb)
 J={n-1-i for i in I}
 for v in range(1,n-1):T.add((v,0) if v in I else (0,v));T.add((n-1,v) if v in J else (v,n-1))
 T.add((n-1,0));return T,index

def ledger(T,n,I):
 a,b,x,y,c,z=0,1,2,n-3,n-2,n-1; W=tuple(range(2,n-2));D=tuple(range(3,n-2))
 if (b,c) in T:
  if y not in I:return 'open_release',((b,)+W+(c,z,a),(z,)+W+(c,b,a),(a,z)+W+(c,b))
  if c not in I:return 'open_penultimate',((a,b,c,z,x)+D,(b,c,a,z,x)+D,(b,x)+D+(c,a,z))
  h=next(iter(I-{y,c}));cut=y-1 if h==y-1 else y
  L,R=tuple(range(x,cut)),tuple(range(cut,y+1))
  return 'open_endrun',((a,b,c,z)+L+R,(b,c,z)+L+(a,)+R,(b,)+L+(a,)+R+(c,z))
 if c in I:
  if y not in I:
   cut=next(i for i in range(2,c) if (z,i) in T);L,R=tuple(range(1,cut)),tuple(range(cut,c))
   return 'closed_last_movable',((a,)+L+(z,)+R+(c,),L+(z,)+R+(a,c),L+R+(a,c,z))
  h=next(iter(I-{y,c}))
  if (x,c) in T and h not in {x,n-4}:
   Q=(a,x)+D+(c,z,b);S=(a,b,x,c)+D[::-1]+(z,);V=(b,z)+D+(c,x,a)
   return 'closed_endpair_guarded',(Q,S[::-1],V)
  return None
 if (c,x) in T:return 'closed_zero',((c,)+W+(b,z,a),(z,c)+W+(b,a),(a,z,b)+W+(c,))
 if x not in I:return 'closed_one',((c,b,z,a)+W,(z,c,b,a)+W,(a,z,b)+W+(c,))
 if y not in I:
  cut=next(i for i in range(3,c) if (z,i) in T);name='closed_second_movable'
 elif 3 not in I:cut=y;name='closed_second_fixed'
 else:cut=n-4;name='closed_second_tail4'
 L,R=tuple(range(2,cut)),tuple(range(cut,c))
 return name,((a,b)+L+(z,)+R+(c,),(a,)+L+(z,)+R+(b,c),(a,)+L+R+(b,c,z))

def graph(T,q):
 return {(i,j-1) for i in range(len(q)) for j in range(i+2,len(q)) if (q[j],q[i]) in T}
def replay(n,I,T,strict,counts,digest):
 P=tuple(range(n));need(len(T)==n*(n-1)//2 and directed(T,P) and blue(T,P),'source')
 need(sum(v==0 for u,v in T)==4 and sum(u==n-1 for u,v in T)==4,'endpoint4')
 G=graph(T,P);reflect=lambda H:{(n-2-v,n-2-u) for u,v in H};host=set(combinations(range(n-1),2))
 need(reflect(G)==G and graph(move(T,P),P)==host-G,'trijection complement/reflection')
 got=ledger(T,n,I)
 if got is None:counts['unclaimed']+=1;return
 name,words=got
 if not strict and name not in {'closed_second_movable','open_endrun','closed_second_tail4'}:
  counts['below_threshold_not_replayed']+=1;return
 if not strict and n<8:
  counts['below_threshold_not_replayed']+=1;return
 now=T;xor=set();flags=[]
 for q in words:
  need(set(q)==set(P) and directed(now,q),('direction',n,I,name,q))
  flags.append(blue(now,q));xor^=edges(q);now=move(now,q)
 need(xor==edges(P) and now==move(T,P),('literal target',n,I,name))
 if strict:need(not any(flags),('blackness',n,I,name,flags))
 counts[name]+=1;counts['small_blue_steps']+=sum(flags) if not strict else 0
 digest.update(json.dumps([n,sorted(I),name,words,flags,sorted(T)],separators=(',',':')).encode())

# Two independent solvers of the same affine *edge* equations.
# A fixed node 0 carries value zero; node i+1 is the forward arc bit on pair i.
def equations(n,kind):
 pairs=list(combinations(range(n),2));idx={p:i+1 for i,p in enumerate(pairs)};eq=[]
 def node(u,v):return idx[edge(u,v)]
 def fixed(u,v,bit,label):eq.append((node(u,v),0,int(bit)^int(u>v),label))
 for i in range(n-1):fixed(i,i+1,1,'path')
 for u,v in pairs:eq.append((node(u,v),node(n-1-v,n-1-u),0,'blueP'))
 fixed(n-1,0,1,'closed')
 a,b,x,y,c,z=0,1,2,n-3,n-2,n-1;D=tuple(range(3,n-2));W=tuple(range(2,n-2))
 if kind=='open':
  q1=(a,b,c,z,x)+D;q2=(b,c,a,z,x)+D
  fixed(b,c,1,'open middle');fixed(a,y,0,'penultimate port');fixed(a,c,1,'last nonport')
 else:
  q1=(c,b,z,a)+W;q2=(z,c,b,a)+W
  fixed(b,c,0,'closed middle');fixed(x,c,1,'epsilon one');fixed(a,x,1,'second nonport');fixed(a,c,1,'last nonport')
 E=edges(q1);r=dict(zip(q2,q2[::-1]))
 for u,v in pairs:
  a1,b1=r[v],r[u];im=edge(a1,b1)
  eq.append((node(u,v),idx[im],int(a1>b1)^int((u,v) in E)^int(im in E),'blueQ2'))
 return pairs,idx,eq

def project_dsu(n,kind):
 pairs,idx,eq=equations(n,kind);p=list(range(len(pairs)+1));d=[0]*len(p)
 def find(i):
  if p[i]!=i:
   root,sgn=find(p[i]);d[i]^=sgn;p[i]=root
  return p[i],d[i]
 for u,v,s,_ in eq:
  ru,du=find(u);rv,dv=find(v)
  if ru==rv:
   if du^dv!=s:return {'inconsistent':True}
  else:p[ru]=rv;d[ru]=du^dv^s
 root0,sign0=find(0);groups={}
 for v in range(2,n-1):
  root,sgn=find(idx[0,v]);groups.setdefault(root,[]).append((v,sgn))
 cost=0
 for root,vs in groups.items():
  if root==root0:cost+=sum(s==sign0 for v,s in vs)
  else:
   ones=sum(s for v,s in vs);cost+=min(ones,len(vs)-ones)
 return {'inconsistent':False,'minimum_middle_inports':cost}

def project_bfs(n,kind):
 pairs,idx,eq=equations(n,kind);G=[[] for _ in range(len(pairs)+1)]
 for u,v,s,_ in eq:G[u].append((v,s));G[v].append((u,s))
 value={};component={};groups={};root=0
 for start in range(len(G)):
  if start in value:continue
  value[start]=0;component[start]=root;queue=deque([start])
  while queue:
   u=queue.popleft()
   for v,s in G[u]:
    t=value[u]^s
    if v in value:
     if value[v]!=t:return {'inconsistent':True}
    else:value[v]=t;component[v]=root;queue.append(v)
  root+=1
 for v in range(2,n-1):
  i=idx[0,v];groups.setdefault(component[i],[]).append(value[i])
 cost=0
 for r,vs in groups.items():cost+=vs.count(0) if r==component[0] else min(vs.count(0),vs.count(1))
 return {'inconsistent':False,'minimum_middle_inports':cost}

def check_port_chains(n,kind):
 pairs,idx,eq=equations(n,kind)
 # Only the two blue involutions, no degree, path, or branch equations.
 graph_eq={}
 for u,v,s,label in eq:
  if label in {'blueP','blueQ2'}:graph_eq[u,v,label]=s
 def value(i,j):return idx[edge(i,j)]
 if kind=='open':
  chains=[[(0,4),(n-5,n-1),(0,n-6)],[(0,8),(n-9,n-1),(4,n-6),(5,n-5),(0,n-10)],[(0,12),(n-13,n-1),(8,n-6),(5,n-9),(4,n-10),(9,n-5),(0,n-14)]]
 else:
  chains=[[(0,5),(n-6,n-1),(0,n-3)],[(0,9),(n-10,n-1),(5,n-3),(2,n-6),(0,n-7)],[(0,13),(n-14,n-1),(9,n-3),(2,n-10),(5,n-7),(6,n-6),(0,n-11)]]
 for chain in chains:
  parity=0
  for t,(a,b) in enumerate(zip(chain,chain[1:])):
   label='blueP' if t%2==0 else 'blueQ2';u,v=value(*a),value(*b)
   need((u,v,label) in graph_eq,('symbolic chain',n,kind,a,b,label));parity^=graph_eq[u,v,label]
  need(parity==1,('opposite endpoint bits',n,kind))
 # Disjoint endpoint pairs; forced fourth port is y for open, 4 for closed.
 ends=[v for c in chains for u,v in (c[0],c[-1])];forced=n-3 if kind=='open' else 4
 need(len(set(ends+[forced]))==7,('port budget disjointness',n,kind))


def cycles(p):
 seen=set();out=[]
 for i in range(len(p)):
  if i in seen:continue
  c=[];v=i
  while v not in seen:seen.add(v);c.append(v);v=p[v]
  out.append(c)
 return out

def permutation_controls():
 for n in range(7,201):
  P=tuple(range(n));q=(0,)+tuple(range(2,n-2))+(1,n-2,n-1);r=dict(zip(q,q[::-1]));g=[r[n-1-i] for i in range(n)];cs=cycles(g)
  target=sorted([1,1,n-2] if n%2 else [1,1,(n-2)//2,(n-2)//2])
  need(sorted(map(len,cs))==target,'relocated-b orbit identity')
  lengths=[len(c) for c in cs if len(c)>1]
  # An invariant set of size four inside the middle would be a union of
  # odd cycles. The possible cycle lengths themselves refute this.
  need(any(l%2==0 for l in lengths) or not any(sum(lengths[j] for j in range(len(lengths)) if mask>>j&1)==4 for mask in range(1<<len(lengths))),'four-port invariant obstruction')
 for n in range(8,201):
  q=tuple(range(n-4))+(n-1,n-4,n-3,n-2);r=dict(zip(q,q[::-1]));g=[r[n-1-i] for i in range(n)]
  need(sorted(len(c) for c in cycles(g) if len(c)>1)==[4,4],'tail-four commutator')

def hostile_controls():
 rng=random.Random(71);out={}
 n=9;I={n-4,n-3,n-2};T,d=source(n,I,rng,0);T,_=source(n,I,rng,(1<<d)-1)
 q1=(0,1,n-2,n-1)+tuple(range(2,n-2))
 oldq2=(1,n-2,n-1)+tuple(range(2,n-3))+(0,n-3)
 need(directed(T,q1) and not directed(move(T,q1),oldq2),'old endpair insertion hostile')
 out['old_endpair_insertion']={'n':n,'ports':sorted(I),'rows':[sum(1<<v for u,v in T if u==i) for i in range(n)],'illegal_join':[n-4,0]}
 n=11;I={2,3,4}
 for attempt in range(1000):
  T,_=source(n,I,rng)
  if (n-2,1) in T and (2,n-2) in T:break
 else:raise RuntimeError('hostile construction')
 oldq1=tuple(range(n-4))+(n-1,n-4,n-3,n-2)
 need(not directed(T,oldq1),'fixed adjacent-pair cut hostile')
 got=ledger(T,n,I);need(got[0]=='closed_second_movable','movable repair selected')
 now=T
 for q in got[1]:need(directed(now,q) and not blue(now,q),'hostile repaired');now=move(now,q)
 need(now==move(T,tuple(range(n))),'hostile literal target')
 out['old_adjacent_pair_cut']={'n':n,'ports':sorted(I),'rows':[sum(1<<v for u,v in T if u==i) for i in range(n)],'illegal_join':[n-5,n-1]}
 return out

def main():
 digest=hashlib.sha256();counts=Counter();rng=random.Random(20261021)
 for n in range(7,10):
  _,dim=source(n,{2,3,4},rng,0)
  for I in combinations(range(2,n-1),3):
   for mask in range(1<<dim):
    T,_=source(n,set(I),rng,mask)
    # Larger inherited templates have stated thresholds. Check target and
    # legality, but do not demand blackness below n15.
    replay(n,set(I),T,False,counts,digest)
 small=dict(counts);counts.clear()
 for n in range(15,61):
  y,c=n-3,n-2
  ports=[set(rng.sample(range(2,n-1),3)) for _ in range(160)]
  ports += [{y,c,h} for h in (2,3,n-5,n-4)]
  ports += [{2,3,y},{2,3,4},{2,y,n-4}]
  for I in ports:
   for _ in range(2):
    T,_=source(n,I,rng);replay(n,I,T,True,counts,digest)
 finite=[]
 for n in range(15,27):
  for kind in ('open','closed'):
   d,b=project_dsu(n,kind),project_bfs(n,kind);need(d==b,'independent projection equality')
   need(d.get('inconsistent') or d['minimum_middle_inports']>3,('finite endpoint budget',n,kind,d))
   finite.append({'n':n,'branch':kind,**d})
 for n in range(27,201):
  for kind in ('open','closed'):check_port_chains(n,kind)
 permutation_controls()
 receipt={'status':'PASS; structural proofs separate from finite controls','hostile_controls':hostile_controls(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'small_exhaustive_n7_to9':small,'large_deterministic_n15_to60':dict(counts),'finite_boundary_affine_systems':finite,'port_chain_control_orders':[27,200],'permutation_control_orders':[7,200],'replay_digest_sha256':digest.hexdigest()}
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
