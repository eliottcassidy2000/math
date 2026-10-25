import hashlib,json,pathlib,itertools,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parent
SOURCE=ROOT/'inputs'
P=SOURCE/'session54-cycle-packing-bound'
def ck(x,s):
 if not x: raise RuntimeError(s)
ck(hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest()=='5213de5e2193cd43e952d529277bb44eca5053fa4d86d9b9afcee944458bc116','manifest')
def setup(n):
 Q={a*a%n for a in range(1,n)}
 E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q]
 return E,[[i for i,(u,v) in enumerate(E) if (v-u)%n==d] for d in sorted(Q)]
def hp(n,E,B):
 # Independent forward subset DP, storing reachable last-vertex masks.
 adj=[0]*n
 for i,(u,v) in enumerate(E):
  if not(B>>i&1):adj[u]|=1<<v
 dp=[0]*(1<<n)
 for v in range(n):dp[1<<v]=1<<v
 for s in range(1,1<<n):
  ends=dp[s]
  while ends:
   b=ends&-ends;ends-=b;v=b.bit_length()-1
   nxt=adj[v]&~s
   while nxt:
    z=nxt&-nxt;nxt-=z;dp[s|z]|=z
 return bool(dp[-1])
result={'manifest':'PASS','inherited':{},'profile_only_hostiles':{}}
for n in (7,11):
 E,C=setup(n)
 if n==7:Bs=[x['blocker'] for x in json.loads((P/'inputs/P7-blockers.json').read_text())['bijection']]
 else:
  row=next(x for x in json.loads((P/'inputs/Hall-census.json').read_text())['censuses'] if x['n']==n)
  Bs=[x['mask'] for x in row['blockers']]+[sum(1<<i for i,(u,v) in enumerate(E) if c in (u,v)) for c in range(n)]
 ck(len(set(Bs))==len(Bs),'distinct blockers')
 for B in Bs:
  ck(B.bit_count()==n-1 and not hp(n,E,B),'independent blocker check')
  ck(all(sum((B>>i)&1 for i in Cj)==2 for Cj in C),'profile')
 result['inherited'][n]=len(Bs)
 cert=next(x for x in json.loads((P/'result.json').read_text())['paley_certificates'] if x['n']==n)['two_per_cycle_nonblocker']
 B=cert['mask'];p=cert['surviving_path']
 ck(len(p)==n and set(p)==set(range(n)),'hostile spanning')
 ck(all((u,v) in E and not(B>>E.index((u,v))&1) for u,v in zip(p,p[1:])),'hostile surviving')
 ck(hp(n,E,B),'hostile DP')
 result['profile_only_hostiles'][n]='PASS'
 if n==7:
  found=[]
  for pairs in itertools.product(*[list(itertools.combinations(c,2)) for c in C]):
   B=sum(1<<i for pair in pairs for i in pair)
   if not hp(n,E,B):found.append(B)
  ck(set(found)==set(Bs),'P7 completeness independent')
  result['P7_profile_universe']=21**3
  result['P7_all_minima']=len(found)
 else:
  idx={e:i for i,e in enumerate(E)};remaining=set(range(len(E)));orbits=[]
  while remaining:
   i=min(remaining);orbit=[];j=i
   while j not in orbit:
    orbit.append(j);u,v=E[j];j=idx[(3*u%n,3*v%n)]
   remaining-=set(orbit);orbits.append(orbit)
  ck(all(len(o)==5 for o in orbits),'multiplier orbit lengths')
  fixed=[sum(1<<i for o in pair for i in o) for pair in itertools.combinations(orbits,2)]
  fixedblock=[B for B in fixed if not hp(n,E,B)]
  ck(len(fixedblock)==1,'fixed minimum blockers')
  result['P11_multiplier3_fixed_minima']=len(fixedblock)
  result['P11_fixed_ten_arc_candidates']=len(fixed)
# Validate copied source bytes where current producer sources are available.
prov=json.loads((P/'provenance.json').read_text())
for dst,src in prov['source_files'].items():
 ck((P/dst).read_bytes()==(SOURCE/src).read_bytes(),'provenance '+src)
for src,digest in prov['source_manifests'].items():
 ck(hashlib.sha256((SOURCE/src).read_bytes()).hexdigest()==digest,'source manifest '+src)
result['source_provenance']='PASS'
print(json.dumps(result,indent=2,sort_keys=True))
