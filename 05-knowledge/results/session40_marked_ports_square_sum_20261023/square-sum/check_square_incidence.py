#!/usr/bin/env python3
"""Exact square-sum path census, edge usage, and endpoint-conditioned sensor banks.
Run from this directory: python3 check_square_incidence.py > receipt.json
No dependency other than Python standard library. Checks survive python -O.
"""
from math import isqrt
from itertools import combinations
from collections import Counter, defaultdict
import json

def check(p,msg):
 if not p:raise RuntimeError(msg)
def graph(n):return {v:tuple(w for w in range(1,n+1) if v!=w and isqrt(v+w)**2==v+w) for v in range(1,n+1)}
def edges_of(p):return frozenset(tuple(sorted(e)) for e in zip(p,p[1:]))
def paths(n):
 G=graph(n);leaves=[v for v in G if len(G[v])==1]; full=(1<<n)-1; P=[]
 if n>1 and (any(not G[v] for v in G) or len(leaves)>2):return P
 def dfs(p,mask):
  if mask==full:
   if leaves or p[0]<p[-1] or n==1:P.append(tuple(p))
   return
  for w in G[p[-1]]:
   bit=1<<(w-1)
   if not mask&bit:dfs(p+[w],mask|bit)
 for s in (leaves[:1] if leaves else range(1,n+1)):dfs([s],1<<(s-1))
 return P

def minimum_bank(E,X,P):
 """Exact hitting-set search on symmetric differences within endpoint fibres."""
 pairs=[(a,b) for a,b in combinations(range(len(P)),2) if set((P[a][0],P[a][-1]))==set((P[b][0],P[b][-1]))]
 masks=[sum(1<<j for j,(a,b) in enumerate(pairs) if (e in X[a])!=(e in X[b])) for e in E]
 target=(1<<len(pairs))-1;seen=set();cand=[];nodes=[0]
 for i,m in enumerate(masks):
  if m and m not in seen:cand.append((i,m));seen.add(m)
 def solve(rem,k,allowed):
  nodes[0]+=1
  if not rem:return []
  if not k:return None
  if rem.bit_count()>sum(sorted(((rem&m).bit_count() for i,m in allowed),reverse=True)[:k]):return None
  opts=None;z=rem
  while z:
   bit=z&-z;z-=bit;o=[(i,m) for i,m in allowed if m&bit]
   if not o:return None
   if opts is None or len(o)<len(opts):opts=o
  excluded=set()
  for i,m in opts:
   tail=solve(rem&~m,k-1,[(j,q) for j,q in allowed if j!=i and j not in excluded])
   if tail is not None:return [i]+tail
   excluded.add(i)
  return None
 trace=[]
 for k in range(len(E)+1):
  before=nodes[0];sol=solve(target,k,cand);trace.append(dict(k=k,search_nodes=nodes[0]-before,exists=sol is not None))
  if sol is not None:
   sig=[(tuple(sorted((p[0],p[-1]))),tuple(E[i] in x for i in sol)) for p,x in zip(P,X)]
   check(len(set(sig))==len(P),'sensor failed injectivity')
   return dict(minimum=k,bank=[E[i] for i in sol],same_endpoint_pairs=len(pairs),proof_search=trace)

def forced_edge_exclusion(n,edge):
 """Check all endpoint choices by degree propagation and proper-cycle rejection.
 Returns a fully replayable trace, not an assumption about degree-two endpoints.
 """
 G=graph(n);E={tuple(sorted((v,w))) for v in G for w in G[v]};leaves={v for v in G if len(G[v])==1};out=[]
 for endpoints in combinations(G,2):
  if not leaves<=set(endpoints):continue
  A=set(E);F={edge};target={v:1 if v in endpoints else 2 for v in G};trace=[];reason=None
  while reason is None:
   changed=False
   for v in G:
    av={e for e in A if v in e};fv=av&F
    if len(fv)>target[v]:reason=['excess',v,len(fv),target[v]];break
    if len(av)<target[v]:reason=['deficit',v,len(av),target[v]];break
    if len(fv)==target[v] and av-F:
     removed=av-F;A-=removed;trace.append(['delete',v,sorted(removed)]);changed=True
    elif len(av)==target[v] and av-F:
     added=av-F;F|=av;trace.append(['force',v,sorted(added)]);changed=True
   if reason:break
   unseen=set(G)
   while unseen:
    S={min(unseen)};unseen-=S;todo=list(S)
    while todo:
     v=todo.pop();adj={b if a==v else a for a,b in F if v in (a,b)};fresh=adj&unseen;unseen-=fresh;S|=fresh;todo.extend(sorted(fresh))
    if sum(a in S and b in S for a,b in F)>=len(S):reason=['forced_cycle',sorted(S)];break
   if reason:break
   if not changed:reason=['OPEN'];break
  out.append(dict(endpoints=endpoints,trace=trace,contradiction=reason))
 return out

def tree_decode(G, E, endpoints, cotree_bits):
 """Canonical spanning tree, GF(2) leaf elimination, and exact edge recovery."""
 if not G:return frozenset(),0
 root=min(G);parent={root:None};order=[root];tree=set()
 for v in order:
  for w in G[v]:
   if w not in parent:parent[w]=v;order.append(w);tree.add(tuple(sorted((v,w))))
 check(len(order)==len(G),'tree decode requires connected graph')
 chosen=set(cotree_bits);rhs={v:int(v in endpoints) for v in G}
 for a,b in chosen:rhs[a]^=1;rhs[b]^=1
 for v in order[:0:-1]:
  if rhs[v]:chosen.add(tuple(sorted((v,parent[v]))));rhs[parent[v]]^=1
 check(rhs[root]==0,'boundary odd')
 return frozenset(chosen),len(E)-len(tree)

rows=[];allP=[];allX=[]
for n in range(1,33):
 P=paths(n);G=graph(n);E=[(v,w) for v in G for w in G[v] if v<w];X=[edges_of(p) for p in P];C=Counter(e for x in X for e in x);ends=Counter(tuple(sorted((p[0],p[-1]))) for p in P)
 check(all(len(p)==n and len(set(p))==n and all(b in G[a] for a,b in zip(p,p[1:])) for p in P),f'invalid path {n}')
 check(len(set(min(p,p[::-1]) for p in P))==len(P),f'duplicate {n}')
 check(sum(C.values())==len(P)*(n-1),f'edge mass {n}')
 paircounts=[dict(edges=[e,f],count=sum(e in x and f in x for x in X)) for e,f in combinations(E,2)]
 check(sum(z['count'] for z in paircounts)==len(P)*(n-1)*(n-2)//2,f'pair mass {n}')
 birth=Counter()
 for p in P:
  i=p.index(n)
  if n==1:birth[('base',)]+=1
  elif i in (0,n-1):birth[('endpoint',p[1] if i==0 else p[-2])]+=1
  else:birth[('interior',*sorted((p[i-1],p[i+1])))]+=1
 # Independent endpoint-term test from prior full path list; forest terms are certified by deletion bijection.
 if n>=3:
  for a in G[n]:check(birth[('endpoint',a)]==sum(a in (p[0],p[-1]) for p in allP[-1]),f'endpoint term {n},{a}')
 if n>=2 and P:
  # Build the same canonical tree once to select cotree coordinates.
  seen={1};order=[1];tree=set()
  for v in order:
   for w in G[v]:
    if w not in seen:seen.add(w);order.append(w);tree.add(tuple(sorted((v,w))))
  for p,x in zip(P,X):
   recovered,beta=tree_decode(G,E,{p[0],p[-1]},x-tree)
   check(recovered==x,f'parity tree decoder {n}')
 row=dict(n=n,count=len(P),endpoints=[dict(pair=e,count=c) for e,c in sorted(ends.items())],edge_counts=[dict(edge=e,count=C[e],status=('NO_PATH' if not P else 'UNUSED' if not C[e] else 'FORCED' if C[e]==len(P) else 'OPTIONAL')) for e in E],pair_counts=paircounts,birth_terms=[dict(type=k[0],ports=k[1:],count=c) for k,c in sorted(birth.items())],paths=P)
 if n>=25:row['endpoint_sensor']=minimum_bank(E,X,P)
 rows.append(row);allP.append(P);allX.append(X)
expected={1:1,15:1,16:1,17:1,23:3,25:10,26:12,27:35,28:52,29:19,30:20,31:349,32:392}
check({r['n']:r['count'] for r in rows if r['count']}==expected,'count benchmark')
check([r['endpoint_sensor']['minimum'] for r in rows[24:]]==[0,1,3,4,1,1,5,6],'bank benchmarks')
check(all(z['status']=='OPTIONAL' for z in rows[31]['edge_counts']),'n32 full support')
# Every edge born by25: exact active horizons and missing intervals within observed successful orders.
history=[]
for z in rows[24]['edge_counts']:
 e=tuple(z['edge']);events=[]
 for r in rows:
  if r['count'] and max(e)<=r['n']:
   q=next(w for w in r['edge_counts'] if tuple(w['edge'])==e);events.append(dict(n=r['n'],count=q['count'],total=r['count'],status=q['status']))
 active=[q['n'] for q in events if q['count']];changes=[q for i,q in enumerate(events) if i==0 or bool(q['count'])!=bool(events[i-1]['count'])]
 history.append(dict(edge=e,birth=max(e),first_observed_activation=min(active),usage_changes=changes,events=events))
# User's explicit colors; '13:,12' interpreted as13--12; literal3--12 separately rejected.
c1={(13,23),(2,23),(11,25),(24,25),(1,15),(10,15),(8,17),(9,16)}
c2={(3,13),(12,13),(2,14),(11,14),(6,19),(17,19),(5,20),(16,20),(4,21),(15,21),(3,22),(14,22),(1,24),(12,24)}
E25={tuple(z['edge']) for z in rows[24]['edge_counts']};check(c1<=E25 and c2<=E25 and not c1&c2,'colors invalid')
check(c1=={e for e in E25 if rows[max(e)-1]['count']>0},'positive birth colors')
check(c2=={e for e in E25 if rows[max(e)-1]['count']==0 and len(graph(max(e))[max(e)])==2},'negative double birth colors')
check((3,12) not in E25,'3--12 must be rejected')
cert=forced_edge_exclusion(25,(1,15));check(len(cert)==24 and all(x['contradiction'][0]!='OPEN' for x in cert),'unused edge forcing certificate')
print(json.dumps(dict(status='FINITE-EXACT, all Hamiltonian paths up to reversal for1<=n<=32',scope='No claim of permanent use or nonuse beyond32. Sensor minimum is only within global edge-membership-bit banks with actual endpoints supplied.',rows=rows,old_edge_history=history,colors=dict(first=sorted(c1),second=sorted(c2),third=sorted(E25-c1-c2),invalid_literal=[3,12]),edge_1_15_exclusion_at25=cert),indent=2))
