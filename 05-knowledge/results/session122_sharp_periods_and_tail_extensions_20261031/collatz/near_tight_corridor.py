#!/usr/bin/env python3
"""Complete slack<=2 predecessor graph, retaining exact original clock/length pairs."""
import mmap,struct,json,sys
from collections import deque
E=int(sys.argv[1]);R=int(sys.argv[3]) if len(sys.argv)>3 else 2
B=int(sys.argv[4]) if len(sys.argv)>4 else 2
LIMIT=int(sys.argv[5]) if len(sys.argv)>5 else 1000000
END=int(sys.argv[6]) if len(sys.argv)>6 else E
if not (0<=B<=2):raise ValueError('only exact entry bank for slack<=2 implemented')
q=5**E;pw=[5**i for i in range(E+1)]
def val(x):
 if not x:return E
 k=0
 while x%5==0:x//=5;k+=1
 return min(k,E)
base=[];width=[];size=0
for h in range(E+1):
 base.append(size);width.append(4*5**(max(R,h)-1))
 size+=(1 if h==E else 4*q//pw[h+1])*width[-1]
def index(z,m):
 h=val(z);j=0 if not z else z//pw[h]
 return base[h]+(0 if not z else (j-1-j//5)*width[h])+m-1-m//5
f=open(sys.argv[2],'rb');data=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)
if len(data)!=8*size:raise ValueError('length')
def V(z,m):return struct.unpack_from('<q',data,8*index(z,m))[0]
def ones(k):return sum(k//(2*5**(e-1)) for e in range(1,E+1))
entries={};entrytests=0;longtests=0
for a in range(2,65):
 inv=pow(pow(2,a,q),-1,q);z=(1-2*inv)%q;h=val(z);ell=max(R,h);mod=5**ell;T=2*5**(ell-1)
 for k in range(T):
  m=pow(3*pow(2,-1,mod)%mod,k,mod)*3*inv%mod
  reward=min(h,val(m-1));score=8*(ones(k)+reward)-5*(k+a);gap=V(z,m)-score
  if gap<0:raise ValueError('invalid entry certificate')
  if gap<=B:entries.setdefault((z,m,gap),[]).append((k+a,k+1,k,a,None))
  if gap+2<=B:
   for e in range(ell,E+1):
    K=k+2*5**(e-1)
    score2=8*(ones(K)+reward)-5*(K+a)
    if score2!=score-2:raise ValueError('entry phase score')
    m2=pow(3*pow(2,-1,mod)%mod,K,mod)*3*inv%mod
    if m2!=m:raise ValueError('entry phase state')
    entries.setdefault((z,m,gap+2),[]).append((K+a,K+1,K,a,e));longtests+=1
  entrytests+=1
p2=[pow(2,a,q) for a in range(65)];inv3=pow(3,-1,q)
roots=[(0,1,b) for b in range(B+1)]
nodes={r:[] for r in roots};todo=deque(roots);scanned=0;edgecount=0
while todo:
 z,m,b=todo.popleft();h=val(z);mod=5**max(R,h);target=V(z,m);reward=min(h,val(m-1));pred=[]
 for a in range(1,65):
  want=target-8*reward+5*a
  if want-b>0 or want<-30000:continue
  z0=((p2[a]*(z-1)+2)*inv3)%q;h0=val(z0);mod0=5**max(R,h0);common=min(mod0,mod)
  residue=m*p2[a]*inv3%common
  for m0 in range(residue,mod0,common):
   scanned+=1;gap=want-V(z0,m0)
   if 0<=gap<=b:
    key=(z0,m0,b-gap);pred.append((key,a,gap))
    if key not in nodes:nodes[key]=[];todo.append(key)
 nodes[(z,m,b)]=pred;edgecount+=len(pred)
 if len(nodes)>LIMIT:raise ValueError(f'incomplete resource guard: nodes={len(nodes)} edges={edgecount} scanned={scanned}')
 if len(nodes)%10000==0:print(f'nodes={len(nodes)} pending={len(todo)} scanned={scanned}',file=sys.stderr)
out={k:[] for k in nodes};degree={k:len(v) for k,v in nodes.items()}
pairs={k:{(x[0],x[1]) for x in entries.get(k,[])} for k in nodes}
paths={k:len(entries.get(k,[])) for k in nodes}
for dest,preds in nodes.items():
 for src,a,gap in preds:out[src].append((dest,a))
ready=deque(k for k,v in degree.items() if not v);processed=0
while ready:
 src=ready.popleft();processed+=1
 for dest,a in out[src]:
  pairs[dest].update((A+a,L+1) for A,L in pairs[src]);paths[dest]+=paths[src]
  if len(pairs[dest])>1000000:raise ValueError('pair-set resource guard')
  degree[dest]-=1
  if not degree[dest]:ready.append(dest)
if processed!=len(nodes):raise ValueError('near-tight augmented cycle')
records=[]
for root in roots:
 result=sorted(pairs[root]);survivors=[]
 for A,L in result:
  if not 1<=L<A:raise ValueError('positive non-one word clock')
  if pow(3,L,5**END)==pow(2,A,5**END):survivors.append([A,L])
 records.append(dict(slack=root[2],score=-V(0,1)+root[2],abstract_paths=paths[root],
  pairs=result,multiplier_survivors=survivors))
print(json.dumps(dict(E=E,R=R,slack_budget=B,endpoint_test_exponent=END,
 status='complete near-tight graph; requires separately verified potential and strict omitted-letter bounds',
 terminal_score=V(0,1),nodes=len(nodes),edges=edgecount,scanned_candidates=scanned,
 entries_checked=entrytests,long_entry_controls=longtests,
 retained_entries=sum(len(entries.get(k,[])) for k in nodes),acyclic=True,records=records),indent=2))
