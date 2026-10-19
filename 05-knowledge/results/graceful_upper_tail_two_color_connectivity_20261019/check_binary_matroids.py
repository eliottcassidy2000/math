"""Independent finite binary-matroid controls for the rank-staircase lemma."""
import itertools, random, json, collections
rng=random.Random(20261019)
def span(vs):
    z={0}
    for v in vs: z|={x^v for x in list(z)}
    return frozenset(z)
def rank(vs): return len(span(vs)).bit_length()-1
cases=[]
for s in range(1,5):
  for trial in range(20):
    dimension=s+1
    classes=[]
    for i in range(s):
      while True:
        a=[rng.randrange(1<<dimension) for _ in range(max(2,s-i+1))]
        if rank(a)>=s-i: break
      classes.append(a)
    # Coordinates identify distinct matroid elements even for repeated vectors.
    states=[p for p in itertools.product(*(range(len(a)) for a in classes))
              if rank([classes[i][p[i]] for i in range(s)])==s]
    closures=[span([classes[i][p[i]] for i in range(s)]) for p in states]
    adj=[[] for _ in states]
    for i,p in enumerate(states):
      for j in range(i):
        d=sum(a!=b for a,b in zip(p,states[j]))
        if d==1 or (d==2 and closures[i]==closures[j]):
          adj[i].append(j); adj[j].append(i)
    diameter=0
    for root in range(len(states)):
      dist={root:0}; todo=collections.deque([root])
      while todo:
        v=todo.popleft()
        for w in adj[v]:
          if w not in dist: dist[w]=dist[v]+1; todo.append(w)
      assert len(dist)==len(states), (s,trial,classes)
      diameter=max(diameter,max(dist.values()))
    assert diameter<=2**s-1, (s,trial,diameter)
    cases.append(dict(s=s,trial=trial,classes=classes,states=len(states),diameter=diameter,bound=2**s-1))
summary=dict(seed=20261019,cases=len(cases),all_connected=True,all_diameter_bounds_pass=True,
             total_states=sum(c['states'] for c in cases),
             max_diameter_by_s={s:max(c['diameter'] for c in cases if c['s']==s) for s in range(1,5)},
             controls=cases)
print(json.dumps(summary,indent=2))
