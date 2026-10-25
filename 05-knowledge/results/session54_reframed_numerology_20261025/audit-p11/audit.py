import itertools as I, json, pathlib, collections, hashlib
P=pathlib.Path(__file__).resolve().parent/'producer'
def ck(p,msg):
 if not p: raise RuntimeError(msg)
def solve(n,E,B):
 # Forward reachable endpoints, independently oriented recurrence.
 adj=[0]*n
 for i,(u,v) in enumerate(E):
  if not (B>>i)&1: adj[u]|=1<<v
 reach=[0]*(1<<n)
 for v in range(n): reach[1<<v]=1<<v
 for s in range(1,1<<n):
  ends=reach[s]
  while ends:
   bit=ends&-ends; ends-=bit;u=bit.bit_length()-1
   nxt=adj[u]&~s
   while nxt:
    b=nxt&-nxt;nxt-=b;reach[s|b]|=b
 return bool(reach[-1])
r=json.loads((P/'result.json').read_text());results=[]
for n in (7,11):
 Q={pow(x,2,n) for x in range(1,n)}
 E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q]
 R=next(x for x in r['censuses'] if x['n']==n)
 ck(E==[tuple(e) for e in R['edges']],'arc reconstruction')
 IN=[{u for u,v in E if v==w} for w in range(n)]
 ck(all(len(s)==(n-1)//2 for s in IN),'degree')
 ck(all(len(a&b)==(n-3)//4 for a,b in I.combinations(IN,2)),'pair intersections')
 # All S,T; cost computed as scalar matrix rectangle sum, no producer helper.
 hist={}; masks=set();tested=0;eq=collections.Counter()
 for k in range(2,n+1):
  H=collections.Counter()
  for S in I.combinations(range(n),k):
   SS=set(S);head=[sum(u in SS for u,v in E if v==w) for w in range(n)]
   for T in I.combinations(range(n),k-2):
    cost=k*((n-1)//2)-sum(head[w] for w in T);H[cost]+=1;tested+=1
    if cost==n-1:
     TT=set(T);masks.add(sum(1<<i for i,(u,v) in enumerate(E) if u in SS and v not in TT));eq[k-2]+=1
  old=next(x for x in R['rows'] if x['S_size']==k)
  ck(dict(H)=={int(a):b for a,b in old['cost_histogram'].items()},'full cost histogram')
  hist[k]=min(H)
 ck(masks=={b['mask'] for b in R['blockers']},'exact masks')
 stars={sum(1<<i for i,e in enumerate(E) if c in e) for c in range(n)}
 ck(not masks&stars,'star separation')
 restores=0
 for B in masks|stars:
  ck(not solve(n,E,B),'blocker')
  for i in range(len(E)):
   if B>>i&1:
    ck(solve(n,E,B^(1<<i)),'restoration');restores+=1
 for c,p in R['star_factor_cycles'].items():
  ck(set(p)==set(range(n))-{int(c)} and len(p)==n-1,'cycle vertices')
  ck(all(e in E for e in zip(p,p[1:]+p[:1])),'cycle arcs')
 results.append(dict(n=n,tested=tested,hall=len(masks),examples=len(masks|stars),restoration_checks=restores,minimum_by_k=hist,equality_certificate_count_by_t=dict(eq)))
 if n==11:
  EI={e:i for i,e in enumerate(E)};perm=[EI[(4*u%11,4*v%11)] for u,v in E]
  unseen=set(range(55));orbits=[]
  while unseen:
   a=min(unseen);orb={a};b=perm[a]
   while b!=a:orb.add(b);b=perm[b]
   unseen-=orb;orbits.append(sum(1<<i for i in orb))
  ck(len(orbits)==11 and all(o.bit_count()==5 for o in orbits),'orbit partition')
  F=json.loads((P/'fixed-points.json').read_text())
  ck(orbits==F['invariant_arc_orbits'],'fixed orbit data')
  rows=F['rows']+F['invariant_five_arc_survivors'];blocked=[]
  ck({x['mask'] for x in F['rows']}=={a|b for a,b in I.combinations(orbits,2)},'all invariant tens')
  ck({x['mask'] for x in F['invariant_five_arc_survivors']}==set(orbits),'all invariant fives')
  for row in rows:
   B=row['mask'];p=row['hamiltonian_path'];yes=solve(n,E,B)
   ck(yes==bool(p),'independent HP test')
   if yes:
    ck(len(p)==n and set(p)==set(range(n)),'path vertices')
    ck(all((u,v) in EI and not B>>EI[u,v]&1 for u,v in zip(p,p[1:])),'path arcs')
   else:blocked.append(B)
  star0=sum(1<<i for i,e in enumerate(E) if 0 in e)
  ck(blocked==[star0]==F['fixed_blockers'],'unique fixed blocker')
  results.append(dict(invariant_5_sets=11,invariant_10_sets=55,surviving_paths=65,unique_fixed_blocker=star0))
print(json.dumps({'status':'PASS','independent_checks':results},indent=2,sort_keys=True))
