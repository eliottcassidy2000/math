from collections import defaultdict,Counter

def run(r):
 k=2*r+1; m=2*k; V=set(range(m+1))-{0,k-1}; paths=[]; triples=[]
 # Independent edge-sum construction including infinity incidence.
 source={}
 for u in range(m+1):
  for v in range(u+1,m+1):
   i=(u%k) if v==m else ((u+v+1)//2)%k
   source[u,v]=i
 cycles=[{e for e,i in source.items() if i==j} for j in range(k)]
 for i,E in enumerate(cycles):
  adj=defaultdict(set)
  for u,v in E:adj[u].add(v);adj[v].add(u)
  assert len(E)==m+1 and all(len(adj[u])==2 for u in range(m+1))
  seen={m}; prev=m;cur=next(iter(adj[m]))
  while cur!=m:
   assert cur not in seen;seen.add(cur);prev,cur=cur,next(iter(adj[cur]-{prev}))
  assert len(seen)==m+1
  word=[m,i]+[x%m for h in range(1,k) for x in (i-h,i+h)]+[(i-k)%m]
  assert {tuple(sorted((word[h-1],word[h]))) for h in range(m+1)}==E
 for j in range(r+1,2*r+1):triples.append((j,j-r-1,2*j%k,(2*j-1,2*j)))
 for i,E in enumerate(cycles):
  repaired={e for e in E if set(e)<=V}
  ts=[t for t in triples if i in t[:2]]
  assert len(ts)==(i!=r)
  if ts:
   a,b,c,e=ts[0]; x=0 if i==a else k-1
   assert {v if u==x else u for u,v in E if x in (u,v)}==set(e)
   assert e not in repaired and source[e]==c and c not in (a,b)
   repaired.add(e)
  adj=defaultdict(set)
  for u,v in repaired:adj[u].add(v);adj[v].add(u)
  ends=[v for v in V if len(adj[v])==1]
  assert len(ends)==2
  p=[max(ends)]
  while len(p)<len(V):
   ns=adj[p[-1]]-set(p[-2:]);assert len(ns)==1;p.append(ns.pop())
  assert set(p)==V and len(set(p))==len(V); paths.append(p)
 count=Counter(tuple(sorted(e)) for p in paths for e in zip(p,p[1:])); M={t[3] for t in triples}
 assert len(set(v for e in M for v in e))==2*r
 assert count==Counter({e:3 if e in M else 1 for e in source if set(e)<=V})
 own={v:h for h,t in enumerate(triples) for v in t[:2]}; nxt={h:own.get(t[2],-1) for h,t in enumerate(triples)}
 cyclesF=[]; done=set()
 for h in nxt:
  trail=[];u=h
  while u!=-1 and u not in done and u not in trail:trail.append(u);u=nxt[u]
  if u in trail:cyclesF.append(trail[trail.index(u):])
  done.update(trail)
 active={nxt[h]:triples[h][2] for C in cyclesF for h in C}
 spare={h:next(x for x in triples[h][:2] if x!=v) for h,v in active.items()}
 order=sorted(set(range(k))-set(active.values()))+sorted(active.values())
 dirs=[{tuple(sorted((u,v))):int(u>v) for u,v in zip(p,p[1:])} for p in paths]
 def relation(h,u,v):
  t=triples[h];seq=sorted(t[:3],key=order.index)
  return dirs[u][t[3]]^dirs[v][t[3]]^(seq.index(u)%2)^(seq.index(v)%2)
 def obstruction(C):
  z=0
  for h in C:z^=relation(h,active[h],triples[h][2])
  return z
 before=[obstruction(C) for C in cyclesF]
 for ci,C in enumerate(cyclesF):
  if obstruction(C):
   h=C[-1];s=spare[h];u=active[h];v=triples[h][2]
   old=[obstruction(D) for D in cyclesF]
   order.remove(s);order.insert(min(order.index(u),order.index(v))+1,s)
   new=[obstruction(D) for D in cyclesF]
   assert [a^b for a,b in zip(old,new)]==[int(i==ci) for i in range(len(cyclesF))]
 assert not any(obstruction(C) for C in cyclesF)
 # All three pairwise constraints (not producer's two chosen constraints).
 adj=defaultdict(list)
 for h,t in enumerate(triples):
  for u,v in [(t[0],t[1]),(t[1],t[2]),(t[2],t[0])]:
   z=relation(h,u,v);adj[u].append((v,z));adj[v].append((u,z))
 sol={}
 for root in range(k):
  if root in sol:continue
  sol[root]=0; todo=[root]
  while todo:
   u=todo.pop()
   for v,z in adj[u]:
    if v not in sol:sol[v]=sol[u]^z;todo.append(v)
    assert sol[v]==sol[u]^z
 seq=[paths[i][::-1] if sol[i]^(j%2) else paths[i] for j,i in enumerate(order)]
 es=sorted(count);idx={e:i for i,e in enumerate(es)};ALL=(1<<len(es))-1; masks=[]; fwds=[]
 for p in seq:
  mask=fwd=0
  for u,v in zip(p,p[1:]):
   bit=1<<idx[tuple(sorted((u,v)))];mask|=bit
   if u>v:fwd|=bit
  masks.append(mask);fwds.append(fwd)
 initial=known=delta=0
 for mask,fwd in zip(masks,fwds):
  required=(fwd^delta)&mask
  assert (initial^required)&known&mask==0
  initial|=required & ~known;known|=mask;delta^=ALL^mask
 assert known==ALL and delta==0
 t=initial;seen=set()
 for mask,fwd in zip(masks,fwds):
  assert t not in seen;seen.add(t);assert t&mask==fwd;t^=ALL^mask
 assert t==initial
 return [len(C) for C in cyclesF],sum(before)
if __name__=='__main__':
 for r in list(range(1,61))+[75,100]:
  cs,flips=run(r);print(f'PASS r={r} n={4*r+1} cycles={cs} parity_flips={flips}',flush=True)
