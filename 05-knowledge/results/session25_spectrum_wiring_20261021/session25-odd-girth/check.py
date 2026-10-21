from collections import Counter,defaultdict

def edge(u,v):return tuple(sorted((u,v)))
def construct(r):
 k=2*r+1; mod=2*k; inf=mod; removed={0,k-1}; paths={}; triples=[]; pairs=[]
 cycles={}
 for i in range(k):
  word=[inf,i]
  for h in range(1,k):word += [(i-h)%mod,(i+h)%mod]
  word += [(i-k)%mod]
  assert len(set(word))==mod+1
  cycles[i]={edge(word[t],word[(t+1)%len(word)]) for t in range(len(word))}
 assert Counter(e for es in cycles.values() for e in es)==Counter({(u,v):1 for u in range(mod+1) for v in range(u+1,mod+1)})
 repairs=defaultdict(list)
 for j in range(r+1,2*r+1):
  a,b=j,(j+r)%k;c=2*j%k;e=(2*j-1,2*j)
  assert len({a,b,c})==3 and e in cycles[c]
  repairs[a].append(e);repairs[b].append(e);triples.append((a,b,c,e));pairs.append((a,b))
 assert set(sum((list(p) for p in pairs),[]))==set(range(k))-{r}
 for i in range(k):
  es={e for e in cycles[i] if not(set(e)&removed)}
  for e in repairs[i]: assert e not in es;es.add(e)
  adj=defaultdict(list)
  for u,v in es:adj[u].append(v);adj[v].append(u)
  ends=[u for u in adj if len(adj[u])==1]
  assert len(ends)==2 and all(len(v)<=2 for v in adj.values())
  word=[min(ends)];prev=None
  while True:
   ns=[v for v in adj[word[-1]] if v!=prev]
   if not ns:break
   prev=word[-1];word.append(ns[0])
   assert len(word)<=mod-1
  assert len(word)==mod-1 and len(set(word))==mod-1
  paths[i]=word
 count=Counter(edge(u,v) for p in paths.values() for u,v in zip(p,p[1:]))
 vertices=set(range(mod+1))-removed
 assert count==Counter({(u,v):(3 if (u,v) in {t[3] for t in triples} else 1) for u in vertices for v in vertices if u<v})
 # Functional graph on repair pairs; -1 is base.
 owner={v:h for h,p in enumerate(pairs) for v in p}
 nxt=[owner.get(t[2],-1) for t in triples]
 cyc=[];seen=set()
 for h in range(r):
  chain=[];where={};v=h
  while v!=-1 and v not in seen and v not in where:
   where[v]=len(chain);chain.append(v);v=nxt[v]
  if v in where:cyc.append(chain[where[v]:])
  seen.update(chain)
 active={}
 for cc in cyc:
  for h in cc:active[nxt[h]]=triples[h][2]
 spares={h:next(v for v in pairs[h] if v!=active[h]) for h in active}
 order=[v for v in range(k) if v not in set(active.values())]+list(active.values())
 directions={i:{edge(u,v):int(u>v) for u,v in zip(p,p[1:])} for i,p in paths.items()}
 def relation(h,a,b):
  t=triples[h];e=t[3];rank={v:p for p,v in enumerate(sorted(t[:3],key=order.index))}
  return directions[a][e]^directions[b][e]^(rank[a]%2)^(rank[b]%2)
 def obstruction(cc):
  ans=0
  for h in cc:ans ^= relation(h,active[h],triples[h][2])
  return ans
 for cc in cyc:
  if obstruction(cc):
   h=cc[0];s=spares[h];a=active[h];b=triples[h][2]
   order.remove(s);lo=min(order.index(a),order.index(b));order.insert(lo+1,s)
  assert obstruction(cc)==0
 # Independent full constraint propagation, using actual order.
 adj=defaultdict(list)
 for h,(a,b,c,e) in enumerate(triples):
  for v in (b,c):
   z=relation(h,a,v);adj[a].append((v,z));adj[v].append((a,z))
 sol={}
 for root in range(k):
  if root in sol:continue
  sol[root]=0;stack=[root]
  while stack:
   u=stack.pop()
   for v,z in adj[u]:
    x=sol[u]^z
    if v in sol:assert sol[v]==x
    else:sol[v]=x;stack.append(v)
 directed=[]
 for pos,i in enumerate(order):directed.append(paths[i][::-1] if sol[i]^(pos%2) else paths[i])
 # Derive initial tournament from all occurrences; verify every directed step.
 initial={};earlier=Counter()
 for pos,p in enumerate(directed):
  for u,v in zip(p,p[1:]):
   e=edge(u,v);z=int(u>v)^((pos-earlier[e])%2)
   if e in initial:assert initial[e]==z
   else:initial[e]=z
   earlier[e]+=1
 state=initial.copy();states=[tuple(sorted(state.items()))]
 for p in directed:
  keep={edge(u,v) for u,v in zip(p,p[1:])}
  assert all(state[edge(u,v)]==int(u>v) for u,v in zip(p,p[1:]))
  state={e:z^(e not in keep) for e,z in state.items()}
  states.append(tuple(sorted(state.items())))
 assert state==initial and len(set(states[:-1]))==k
 return k,len(cyc),[len(c) for c in cyc]
if __name__=='__main__':
 for r in range(1,51):
  k,c,lens=construct(r);print('PASS n=%d k=%d functional_cycles=%s'%(4*r+1,k,lens),flush=True)
