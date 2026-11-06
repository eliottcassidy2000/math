import itertools,json
checks=0
for n in range(2,6):
 pairs=list(itertools.combinations(range(n),2));ix={p:k for k,p in enumerate(pairs)};full=(1<<n)-1;Hs={};deltas={}
 for code in range(1<<len(pairs)):
  def beat(a,b):return bool(code>>ix[tuple(sorted((a,b)))]&1)==(a<b)
  start={};end={}
  for S in range(1,1<<n):
   vs=[v for v in range(n) if S>>v&1]
   for v in vs:start[S,v]=end[S,v]=0
   for p in itertools.permutations(vs):
    if all(beat(a,b) for a,b in zip(p,p[1:])):start[S,p[0]]+=1;end[S,p[-1]]+=1
  Hs[code]=sum(end[full,v] for v in range(n));ds=[]
  for u,v in pairs:
   rest=full^(1<<u)^(1<<v);s=rest;uv=vu=0
   while True:
    A=s|(1<<u);uv+=end[A,u]*start[full^A,v]
    A=s|(1<<v);vu+=end[A,v]*start[full^A,u]
    if s==0:break
    s=(s-1)&rest
   ds.append(vu-uv if beat(u,v) else uv-vu)
  deltas[code]=ds
 for code,ds in deltas.items():
  for k,d in enumerate(ds):assert Hs[code^(1<<k)]-Hs[code]==d;checks+=1
print(json.dumps({'status':'PASS','universe':'all labelled tournaments orders2..5, every arc reversal','arc_checks':checks,'endpoint_method':'direct permutations of every induced subset'},indent=2))
