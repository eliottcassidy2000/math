import itertools,json,math
from functools import lru_cache

def need(b,m):
 if not b:raise RuntimeError(m)
def words(k):return list(itertools.product(range(3),repeat=k))
def arc(a,b):
 for x,y in zip(a,b):
  if x!=y:return (y-x)%3==1
 return False
def matrix(k):
 w=words(k);return w,[[arc(a,b) for b in w] for a in w]
def rotation(k,prefix):
 w=words(k);idx={a:i for i,a in enumerate(w)};d=len(prefix)
 return tuple(idx[a[:d]+((a[d]+1)%3,)+a[d+1:]] if a[:d]==prefix else i for i,a in enumerate(w))
def auto(A,p):return all(A[i][j]==A[p[i]][p[j]] for i in range(len(A)) for j in range(len(A)))
@lru_cache(None)
def wreath(k):
 if k==0:return ((0,),)
 prev=wreath(k-1);r=3**(k-1)
 return tuple(tuple(((a+c)%3)*r+gs[a][x] for a in range(3) for x in range(r)) for c in range(3) for gs in itertools.product(prev,repeat=3))
def components(m,depths,threshold):
 adj=[set() for _ in range(m)]
 for (i,j),l in depths.items():
  if l>=threshold:adj[i].add(j);adj[j].add(i)
 unseen=set(range(m));out=[]
 while unseen:
  todo=[min(unseen)];comp=set(todo);unseen-=comp
  while todo:
   for v in adj[todo.pop()] & unseen:unseen.remove(v);comp.add(v);todo.append(v)
  out.append(sorted(comp))
 return out

def host(k,m,depths):
 w,A=matrix(k);n=len(w)
 def outer(i,j):return 0<(j-i)%m<=(m-1)//2
 N=m*n;H=[[False]*N for _ in range(N)]
 for u in range(N):
  i,a=divmod(u,n)
  for v in range(N):
   j,b=divmod(v,n)
   if i==j:H[u][v]=A[a][b]
   else:
    l=depths[tuple(sorted((i,j)))];eq=w[a][:l]==w[b][:l]
    H[u][v]=eq if outer(i,j) else not eq
 need(all(H[u][v]!=H[v][u] for u in range(N) for v in range(u+1,N)),'not tournament')
 # Weighted divergence formula, independently compared with actual row sums.
 for i in range(m):
  div=0
  for j in range(m):
   if i==j:continue
   l=depths[tuple(sorted((i,j)))];weight=2*3**(k-l)-n
   div+=weight if outer(i,j) else -weight
  need(all(2*sum(H[i*n+a])==m*n-1+div for a in range(n)),'weighted degree formula')
 generators=0
 comps=[]
 for d in range(k):
  cd=components(m,depths,d+1);comps.append(len(cd))
  for prefix in itertools.product(range(3),repeat=d):
   r=rotation(k,prefix)
   for comp in cd:
    cs=set(comp);p=tuple(i*n+(r[a] if i in cs else a) for i in range(m) for a in range(n))
    need(auto(H,p),'mixed prefix generator not actual automorphism');generators+=1
 exponent=sum(3**d*comps[d] for d in range(k))
 need(generators==exponent,'generator coordinate count')
 P=[i*n+a for i in reversed(range(m)) for a in range(n)]
 path_ok=all(H[a][b] for a,b in zip(P,P[1:]))
 need(path_ok,'reverse lex path not directed')
 fullcomps=components(m,depths,k);which={i:c for c,C in enumerate(fullcomps) for i in C}
 slots=[(which[i],a//3) for i in range(m) for a in range(n)]
 distinct=sorted(set(slots));sid={v:i for i,v in enumerate(distinct)};t=len(distinct);adj=[set() for _ in range(t)];edge_positions={}
 for pos,(a,b) in enumerate(zip(P,P[1:])):
  x,y=sid[slots[a]],sid[slots[b]];edge_positions.setdefault(tuple(sorted((x,y))),pos)
  if x!=y:adj[x].add(y);adj[y].add(x)
 @lru_cache(None)
 def matching(mask):
  if not mask:return ()
  low=mask&-mask;i=low.bit_length()-1;best=matching(mask^low)
  for j in adj[i]:
   if mask>>j&1:
    cand=((i,j),)+matching(mask^low^(1<<j))
    if len(cand)>len(best):best=cand
  return best
 M=matching((1<<t)-1);covered={v for e in M for v in e};flagpos=[edge_positions[tuple(sorted(e))] for e in M]
 for v in range(t):
  if v not in covered:
   edge=next(e for e in edge_positions if v in e);flagpos.append(edge_positions[edge])
 seen={sid[slots[P[p+q]]] for p in flagpos for q in (0,1)}
 need(len(seen)==t and len(set(flagpos))==t-len(M),'edge cover not exact construction')
 return {'k':k,'cells':m,'order':N,'regular':all(sum(row)==(N-1)//2 for row in H),'threshold_component_counts':comps,'chart_order_power_of_3':exponent,'owner_base':t,'arc_matching':len(M),'arc_flag_base':len(flagpos),'arc_flag_positions':flagpos,'slot_labels':distinct,'matching_edges':M,'slot_edges':[list(e) for e in sorted(edge_positions)],'checked_generators':generators}

out={'status':'PASS','internal':[],'uniform_hosts':[]}
for k in range(1,4):
 w,A=matrix(k);n=len(w);r=n//3
 need(all(sum(row)==(n-1)//2 for row in A),'internal not regular')
 need(all(A[i][i+1] for i in range(n-1)),'lex not path')
 for i in range(n):
  for j in range(i+1,n):
   common=sum(A[i][z] and A[j][z] for z in range(n))
   need((common>(r-1)//2)==(w[i][0]==w[j][0]),'intrinsic top classes')
 out['internal'].append({'k':k,'n':n,'full_automorphism_order_formula':3**((n-1)//2)})
 for l in range(1,k+1):
  R=host(k,3,{e:l for e in itertools.combinations(range(3),2)})
  need(R['regular'],'uniform host not regular')
  expected=3*(n//3) if l<k else n//3
  need(R['owner_base']==expected and R['arc_flag_base']==(expected+1)//2,'uniform bases')
  out['uniform_hosts'].append(R)
small=[]
for k in (1,2):
 w,A=matrix(k);n=len(w);G=wreath(k)
 need(len(G)==3**((n-1)//2) and len(set(G))==len(G),'wreath order')
 need(all(auto(A,g) for g in G),'wreath not automorphism')
 for selections in itertools.product(range(3),repeat=n//3):
  marked=[3*b+selections[b] for b in range(n//3)]
  need(sum(all(g[v]==v for v in marked) for g in G)==1,'bottom base failure')
 for l in range(k+1):
  restriction=lambda g:tuple(w[g[i]][:l] for i in range(n))
  identity=restriction(tuple(range(n)))
  ker=[g for g in G if restriction(g)==identity]
  need(len(ker)==3**((n-3**l)//2),'kernel size')
  compatible=0
  B=[[a[:l]==b[:l] for b in w] for a in w]
  for g in G:
   for h in G:
    actual=all(B[i][j]==B[g[i]][h[j]] for i in range(n) for j in range(n))
    need(actual==(restriction(g)==restriction(h)),'prefix alignment criterion')
    compatible+=actual
  need(compatible==3**((3**l-1)//2)*len(ker)**2,'pair order')
  small.append({'k':k,'l':l,'kernel_order':len(ker),'compatible_chart_pairs':compatible})
out['small_exhaustive']=small
# Regular heterogeneous order63 star-slot hostile.
m=7;k=2;C={1,3,5};D={e:(2 if set(e)<=C else 1) for e in itertools.combinations(range(m),2)}
r=host(k,m,D)
need(r['regular'] and r['owner_base']==15 and r['chart_order_power_of_3']==16 and r['arc_matching']==6 and r['arc_flag_base']==9,'order63 hostile')
# Explicit odd-component deficiency witness after deleting central endpoint slots.
center=next(c for c,comp in enumerate(components(m,D,2)) if set(comp)==C)
labels=[tuple(x) for x in r['slot_labels']];deleted={labels.index((center,0)),labels.index((center,2))}
remaining=set(range(len(labels)))-deleted;edges=r['slot_edges'];odd=[]
while remaining:
 start=min(remaining);remaining.remove(start);comp={start};todo=[start]
 while todo:
  u=todo.pop()
  for e in edges:
   if u in e:
    for v in e:
     if v in remaining:remaining.remove(v);comp.add(v);todo.append(v)
 if len(comp)%2:odd.append(sorted(comp))
need(len(odd)==5,'matching deficiency witness')
r['deleted_central_endpoint_slots']=sorted(deleted);r['remaining_odd_components']=odd
out['mixed_order63_hostile']=r
# Positive step-class mixed regular example.
m=9;k=2;D={e:(2 if min((e[1]-e[0])%m,(e[0]-e[1])%m)==3 else 1) for e in itertools.combinations(range(m),2)}
r=host(k,m,D)
need(r['regular'] and r['owner_base']==9 and r['chart_order_power_of_3']==10,'circulant positive')
out['mixed_order81_positive']=r
print(json.dumps(out,indent=2,sort_keys=True))
