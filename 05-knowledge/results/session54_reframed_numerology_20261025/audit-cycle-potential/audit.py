import pathlib,json,itertools as it,hashlib
from fractions import Fraction as F
ROOT=pathlib.Path(__file__).resolve().parent
P=ROOT/'inputs/session54-cycle-potential-classification'
def ck(x,s):
 if not x:raise RuntimeError(s)
ck(hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest()=='6b05de37509f3a6f041362b21c8bb98f09ec8e27c4c9470df90bf0c3433a1693','frozen hash')
def det(A):
 a=[[F(v) for v in row] for row in A];d=F(1)
 for k in range(len(a)):
  p=next((i for i in range(k,len(a)) if a[i][k]),None)
  if p is None:return 0
  if p!=k:a[p],a[k]=a[k],a[p];d=-d
  z=a[k][k];d*=z
  for i in range(k+1,len(a)):
   c=a[i][k]/z
   for j in range(k,len(a)):a[i][j]-=c*a[k][j]
 return d
def hp(n,E,B):
 adj=[0]*n
 for i,(u,v) in enumerate(E):
  if not(B>>i&1):adj[u]|=1<<v
 dp=[0]*(1<<n)
 for v in range(n):dp[1<<v]=1<<v
 for s in range(1,1<<n):
  ends=dp[s]
  while ends:
   b=ends&-ends;ends-=b;v=b.bit_length()-1;nxt=adj[v]&~s
   while nxt:z=nxt&-nxt;nxt-=z;dp[s|z]|=z
 return bool(dp[-1])
def potential(n,E,B):
 # Independently propagate integral potentials through all support edges.
 val={0:0}
 while len(val)<2*n:
  old=len(val)
  for i,(u,v) in enumerate(E):
   v+=n;b=(B>>i)&1
   if u in val:val[v]=val[u]-b if v not in val else val[v]
   if v in val:val[u]=val[v]+b if u not in val else val[u]
  ck(len(val)>old,'support connected')
 ck(all(val[u]-val[n+v]==((B>>i)&1) for i,(u,v) in enumerate(E)),'potential consistency')
 lo=min(val.values());a=[val[u]-lo for u in range(n)];d=[val[n+v]-lo for v in range(n)];M=max(a+d)
 cuts=[];deltas=[]
 for k in range(1,M+1):
  S={u for u in range(n) if a[u]>=k};T={v for v in range(n) if d[v]>=k}
  ck(all(u in S for u,v in E if v in T),'compatible level')
  delta=len(S)-len(T);ck(delta>0,'positive deficiency');deltas.append(delta)
  cuts.append(sum(1<<i for i,(u,v) in enumerate(E) if u in S and v not in T))
 ck(sum(deltas)==2 and sum(cuts)==B,'integral level decomposition')
 ck(all(not x&y for x,y in it.combinations(cuts,2)),'level disjointness')
 return deltas
out=[]
for n in (7,11):
 r=(n-1)//2;lam=(n-3)//4;V=set(range(n));Q={x*x%n for x in range(1,n)};E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];idx={e:i for i,e in enumerate(E)}
 cert=json.loads((P/f'rank-{n}.json').read_text());row=next(x for x in json.loads((P/'result.json').read_text())['classification'] if x['n']==n)
 ck(cert['edges']==[list(e) for e in E],'edge indexing')
 Cs=[]
 for seq in [cert['base_cycle']]+cert['independent_cycles']:
  ck(len(seq)==n and set(seq)==V,'cycle spans');arcs=list(zip(seq,seq[1:]+seq[:1]));ck(all(e in idx for e in arcs),'cycle arcs');Cs.append([int(e in arcs) for e in E])
 A=[[x-y for x,y in zip(c,Cs[0])] for c in Cs[1:]];piv=cert['pivot_columns'];z=det([[v[i] for i in piv] for v in A]);ck(z==(1 if n==7 else -2),'independent Fraction determinant')
 ck(len(A)==len(E)-2*n+1,'rank bound reached')
 for v in V:
  ck(all(sum(a[i] for i,e in enumerate(E) if e[0]==v)==0 and sum(a[i] for i,e in enumerate(E) if e[1]==v)==0 for a in A),'row column orthogonality')
 tree=row['support_spanning_tree'];reach={0}
 for u,v in tree:
  ck(u in reach and v not in reach,'tree ordering')
  ck((u,v-n) in idx if u<n else (v,u-n) in idx,'tree edge');reach.add(v)
 ck(len(reach)==2*n,'support tree')
 IN={v:{u for u,w in E if w==v} for v in V}
 ck(all(len(IN[v])==r for v in V) and all(len(IN[u]&IN[v])==lam for u,v in it.combinations(V,2)),'double regular')
 orbit=set()
 for a in Q:
  for b in V:
   f=lambda x:(a*x+b)%n
   ck({(f(u),f(v)) for u,v in E}==set(E),'affine automorphism')
   orbit.add((f(0),f(1)))
 ck(orbit==set(E),'arc transitivity')
 for path in cert['endpoint_paths']:
  ck(len(path)==n and set(path)==V and all(e in idx for e in zip(path,path[1:])),'endpoint path')
 ck(cert['endpoint_paths'][0][0]==0 and cert['endpoint_paths'][0][-1]==1,'forward witness')
 # Independent all S,T mask enumeration (not the producer neighbor-union enumeration).
 cuts={1:set(),2:set()};sizecounts={1:{},2:{}}
 for sm in range(1<<n):
  for tm in range(1<<n):
   delta=sm.bit_count()-tm.bit_count()
   if delta not in cuts:continue
   if any((tm>>v&1) and not(sm>>u&1) for u,v in E):continue
   B=sum(1<<i for i,(u,v) in enumerate(E) if (sm>>u&1) and not(tm>>v&1));cuts[delta].add(B)
   t=tm.bit_count();sizecounts[delta][t]=sizecounts[delta].get(t,0)+1
 rows=[sum(1<<i for i,(u,v) in enumerate(E) if u==c) for c in V];cols=[sum(1<<i for i,(u,v) in enumerate(E) if v==c) for c in V]
 ck(cuts[1]==set(rows+cols),'deficiency 1 classification')
 two={a|b for a,b in it.combinations(cuts[1],2) if not a&b};candidates=cuts[2]|two
 ck(len(candidates)==(84 if n==7 else 176),'candidate count')
 levels={}
 for B in candidates:
  ds=potential(n,E,B);levels[len(ds)]=levels.get(len(ds),0)+1
  ck(all(sum(v for i,v in enumerate(c) if B>>i&1)==2 for c in Cs),'cycle intersection')
 blocked={B for B in candidates if not hp(n,E,B)}
 ck(blocked==set(row['all_minimum_blockers']),'all blocker masks independent')
 ck(blocked==cuts[2]|{rows[v]|cols[v] for v in V},'Hall plus stars')
 for item in row['mixed_surviving_paths']:
  B=item['mask'];p=item['path'];ck(B in candidates-blocked and set(p)==V and len(p)==n,'mixed path domain')
  ck(all(e in idx and not(B>>idx[e]&1) for e in zip(p,p[1:])),'mixed path arcs')
 ck({x['mask'] for x in row['mixed_surviving_paths']}==candidates-blocked,'all finite paths covered')
 for t in range(1,2*r):
  ck(r*r*t-(t+1)*(lam*t+r-lam)==r*(r-2)+lam*(t-1)*(2*r-1-t)>0,'delta1 polynomial')
 for t in range(1,2*r-1):
  z=r*r*t-(t+2)*(lam*t+r-lam);ck(z==r*(r-3)+lam*(t-1)*(2*r-2-t),'delta2 polynomial');ck(z>0 if n==11 else (z==0)==(t in (1,4)),'delta2 positivity')
 out.append({'n':n,'Fraction_minor':str(det([[v[i] for i in piv] for v in A])),'rank':len(A),'candidate_count':len(candidates),'minimum_count':len(blocked),'cut_certificate_size_counts':sizecounts,'levels':levels,'mixed_paths_verified':len(candidates-blocked)})
print(json.dumps({'status':'PASS','independent':out},indent=2,sort_keys=True))
