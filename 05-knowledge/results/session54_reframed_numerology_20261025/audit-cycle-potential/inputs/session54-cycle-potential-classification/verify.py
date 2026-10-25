#!/usr/bin/env python3
import pathlib,json,itertools as it,collections as co
D=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def determinant(A):
 A=[r[:] for r in A];n=len(A);previous=1;sgn=1
 for k in range(n-1):
  if A[k][k]==0:
   j=next((j for j in range(k+1,n) if A[j][k]),None)
   if j is None:return 0
   A[k],A[j]=A[j],A[k];sgn=-sgn
  pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    z=A[i][j]*pivot-A[i][k]*A[k][j];ck(z%previous==0,'Bareiss exact division');A[i][j]=z//previous
   A[i][k]=0
  previous=pivot
 return sgn*A[-1][-1]
def hp(n,E,B):
 inc=[0]*n
 for j,(u,v) in enumerate(E):
  if not(B>>j&1):inc[v]|=1<<u
 end=[0]*(1<<n)
 for v in range(n):end[1<<v]=1<<v
 for S in range(1,1<<n):
  if S&(S-1)==0:continue
  for v in range(n):
   if S>>v&1 and end[S^(1<<v)]&inc[v]:end[S]|=1<<v
 if not end[-1]:return None
 S=(1<<n)-1;v=(end[-1]&-end[-1]).bit_length()-1;P=[v]
 while S&(S-1):
  S^=1<<v;q=end[S]&inc[v];v=(q&-q).bit_length()-1;P.append(v)
 return P[::-1]
out=[]
for n in (7,11):
 cert=json.loads((D/f'rank-{n}.json').read_text());V=set(range(n));Q={x*x%n for x in V if x};E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)};r=(n-1)//2;lam=(n-3)//4
 ck(cert['edges']==[list(e) for e in E],'edge convention')
 def mask(edges):return sum(1<<EI[e] for e in edges)
 def cyclevec(C):
  ck(len(C)==n and set(C)==V,'cycle vertices');edges=list(zip(C,C[1:]+C[:1]));ck(set(edges)<=set(E),'cycle arcs')
  return [int(e in edges) for e in E]
 base=cyclevec(cert['base_cycle']);cycles=[cyclevec(C) for C in cert['independent_cycles']]
 A=[[x-y for x,y in zip(C,base)] for C in cycles];piv=cert['pivot_columns'];target=len(E)-2*n+1
 ck(len(A)==target and len(piv)==target and len(set(piv))==target,'minor size')
 det=determinant([[row[p] for p in piv] for row in A]);ck(det!=0,'full cycle-difference rank')
 # Connected bipartite support, with explicit BFS tree.
 graph={i:[] for i in range(2*n)}
 for u,v in E:graph[u].append(n+v);graph[n+v].append(u)
 reached={0};todo=[0];tree=[]
 while todo:
  u=todo.pop(0)
  for v in graph[u]:
   if v not in reached:reached.add(v);todo.append(v);tree.append([u,v])
 ck(len(reached)==2*n,'bipartite support connected')
 IN={v:{u for u,w in E if w==v} for v in V}
 ck(all(len(IN[v])==r for v in V),'regular')
 ck(all(len(IN[u]&IN[v])==lam for u,v in it.combinations(V,2)),'double regular')
 # Verify arc transitivity directly; two endpoint orbit witnesses.
 orbit={(b,(a+b)%n) for a in Q for b in V};ck(orbit==set(E),'arc transitive')
 for P in cert['endpoint_paths']:
  ck(len(P)==n and set(P)==V and all(e in EI for e in zip(P,P[1:])),'endpoint path')
 ck(cert['endpoint_paths'][0][0]==0 and cert['endpoint_paths'][0][-1]==1,'forward endpoints')
 ck(cert['endpoint_paths'][1][0]==1 and cert['endpoint_paths'][1][-1]==0,'reverse endpoints')
 # Explicit difference-cycle decomposition supplies beta=2r.
 for d in Q:cyclevec([k*d%n for k in range(n)])
 # All deficiency-one/two compatible cuts, reconstructed independently by neighbor unions.
 cuts={1:{},2:{}}
 for delta in (1,2):
  for t in range(n-delta+1):
   for T0 in it.combinations(sorted(V),t):
    T=set(T0);U=set().union(*(IN[v] for v in T))
    if len(U)>t+delta:continue
    for X in it.combinations(sorted(V-U),t+delta-len(U)):
     S=U|set(X);B=mask([e for e in E if e[0] in S and e[1] not in T]);ck(B.bit_count()==r*delta,'cut cost')
     cuts[delta].setdefault(B,[]).append({'S':sorted(S),'T':sorted(T)})
 ck(len(cuts[1])==2*n,'deficiency one count')
 rows={v:mask([e for e in E if e[0]==v]) for v in V};cols={v:mask([e for e in E if e[1]==v]) for v in V}
 ck(set(cuts[1])==set(rows.values())|set(cols.values()),'deficiency one classification')
 stars={rows[v]|cols[v] for v in V};mixed={rows[v]|cols[w] for v in V for w in V if v!=w and not(rows[v]&cols[w])}
 candidates=set(cuts[2])|stars|mixed
 ck(len(candidates)==(84 if n==7 else 176),'complete potential candidate count')
 blocked=[];surviving=[]
 for B in sorted(candidates):
  ck(B.bit_count()==2*r and all(sum(C[j] for j in range(len(E)) if B>>j&1)==2 for C in [base]+cycles),'constant cycle intersection')
  P=hp(n,E,B)
  if P:surviving.append({'mask':B,'path':P})
  else:blocked.append(B)
 ck(set(blocked)==set(cuts[2])|stars,'full minimum classification')
 ck(len(blocked)==(63 if n==7 else 121),'full count')
 # Validate nested two-level possibilities directly, rejecting duplicates/overlap.
 nested=set()
 for A,ca in cuts[1].items():
  for B,cb in cuts[1].items():
   if A&B:continue
   if any(set(y['S'])<=set(x['S']) and set(y['T'])<=set(x['T']) for x in ca for y in cb):nested.add(A|B)
 ck(nested==stars|mixed,'nested level classification')
 out.append({'n':n,'beta':2*r,'minimum_blocker_count':len(blocked),'Hall_minimum_count':len(cuts[2]),'stars':len(stars),'mixed_survivors':len(surviving),'cycle_difference_rank':target,'exact_minor_determinant':det,'support_spanning_tree':tree,'potential_candidate_count':len(candidates),'Hall_certificates':[{'mask':B,'certificates':cs} for B,cs in sorted(cuts[2].items())],'all_minimum_blockers':blocked,'mixed_surviving_paths':surviving})
print(json.dumps({'status':'PASS','classification':out},indent=2,sort_keys=True))
