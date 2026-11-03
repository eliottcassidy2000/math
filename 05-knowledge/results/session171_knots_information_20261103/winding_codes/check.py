import itertools,json,collections
SEEDS=(153,8345,8337,8369,177,185)
PAIRS=list(itertools.combinations(range(6),2))
def edge(m,a,b):
 i=PAIRS.index(tuple(sorted((a,b))))
 return bool((m>>i)&1)==(a<b)
def hp(m):
 return sum(all(edge(m,p[i],p[i+1]) for i in range(5)) for p in itertools.permutations(range(6)))
def need(p,msg):
 if not p:raise RuntimeError(msg)
def span(g):
 s={0}
 for v in g:s|={x^v for x in tuple(s)}
 return frozenset(s)
def syndrome(w,rows):
 return sum(((w&r).bit_count()%2)<<i for i,r in enumerate(rows))
def report(N,rows):
 C=[w for w in range(1<<N) if not syndrome(w,rows)]
 nonzero=[w for w in C if w]
 d=min(map(int.bit_count,nonzero)) if nonzero else None
 # Literal masks: independent arc Hamming distance, not a winding computation.
 measured=[]
 for w in C:
  distance=sum((SEEDS[0]^SEEDS[3 if w>>i&1 else 0]).bit_count() for i in range(N))
  need(distance==3*w.bit_count(),'literal metric mismatch')
  measured.append(distance)
 if d:need(min(x for x in measured if x)==3*d,'code systole')
 return dict(N=N,rows=rows,code_dimension=(len(C).bit_length()-1),distance=d,first_hidden_flips=3*d if d else None,codewords=C,H=29**N)
# Exact seed Hamilton counts by permutation enumeration.
need([hp(m) for m in SEEDS]==[29]*6,'seed H')
need(all((SEEDS[h]^SEEDS[(h+1)%6]).bit_count()==1 for h in range(6)),'edges')
need(all((SEEDS[h]^SEEDS[(h+3)%6]).bit_count()==3 for h in range(6)),'antipodes')
# All linear subspaces in F2^3; BFS computes first hidden physical return
# without using the formula 3d. At logical return all coordinates in {0,3}.
spaces={span(g) for k in range(4) for g in itertools.combinations(range(1,8),k)}
need(len(spaces)==16,'subspace census')
bfs=[]
for C in sorted(spaces,key=lambda s:(len(s),sorted(s))):
 start=(0,0,0); q=collections.deque([start]); dist={start:0}; hit=None
 while q:
  h=q.popleft()
  if all(x%3==0 for x in h):
   w=sum((x//3)<<i for i,x in enumerate(h))
   if w and w in C:hit=dist[h];break
  for i in range(3):
   for e in (-1,1):
    z=list(h);z[i]=(z[i]+e)%6;z=tuple(z)
    if z not in dist:dist[z]=dist[h]+1;q.append(z)
 expected=3*min(x.bit_count() for x in C if x) if len(C)>1 else None
 need(hit==expected,'BFS hidden-loop distance')
 bfs.append(dict(codewords=sorted(C),first_hidden_flips=hit))
# Hamming columns are all nonzero triples, hence every 1-error syndrome unique.
hamrows=[sum(((j>>i)&1)<<(j-1) for j in range(1,8)) for i in range(3)]
need(len({syndrome(1<<i,hamrows) for i in range(7)}|{0})==8,'one-error correction')
cases=[report(7,[]),report(7,hamrows),report(7,[127]),report(7,[(1<<i)|(1<<6) for i in range(6)]),report(7,[1<<i for i in range(7)])]
print(json.dumps(dict(seed_H=[29]*6,all_F2_3_subspace_BFS=bfs,cases=cases),indent=2,sort_keys=True))
