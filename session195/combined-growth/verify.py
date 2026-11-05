#!/usr/bin/env python3
"""Self-contained exact verifier: no imports from repository or peer packets."""
import argparse,hashlib,itertools,json,math
from collections import deque
from pathlib import Path
HERE=Path(__file__).resolve().parent

def check(ok,why):
 if not ok:raise RuntimeError(why)

def graphs(n):
 for m in range(1<<(n*(n-1)//2)):
  a=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if m>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
    k+=1
  yield a

def dp(a):
 n=len(a);N=1<<n;H=[0]*N;H[0]=1;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):
  check(0<=a[i]<N and not(a[i]>>i&1),'invalid row')
  for j in range(i):check(((a[i]>>j)^(a[j]>>i))&1,'invalid pair')
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    if not rest:E[mask][v]=S[mask][v]=1
    else:
     E[mask][v]=sum(E[rest][u] for u in range(n) if a[u]>>v&1)
     S[mask][v]=sum(S[rest][u] for u in range(n) if a[v]>>u&1)
  H[mask]=sum(E[mask]);check(H[mask]==sum(S[mask]),'start/end mismatch')
 full=N-1;out=[]
 for r in range(n):
  deleted=full^(1<<r)
  L=sum(E[u][r]*H[full^u] for u in range(N) if u>>r&1)
  R=sum(S[u][r]*H[full^u] for u in range(N) if u>>r&1)
  B=sum(H[u]*H[deleted^u] for u in range(N) if u&deleted==u)
  out.append((H[-1],S[-1][r],L,E[-1][r],R,B,int(n==1),H[deleted]))
 return out

def literal(a,r):
 n=len(a);h=s=e=L=R=B=d=0
 def path(p):return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  if path(p):h+=1;s+=p[0]==r;e+=p[-1]==r
  for k in range(n+1):
   if path(p[:k]) and path(p[k:]):
    L+=bool(k and p[k-1]==r);R+=bool(k<n and p[k]==r)
 for p in itertools.permutations([i for i in range(n) if i!=r]):
  d+=path(p)
  for k in range(n):B+=path(p[:k]) and path(p[k:])
 return h,s,L,e,R,B,int(n==1),d

def state(p):return p[0],p[3],p[4],p[7]
def cone(x):
 h,e,R,d=x;return e,d-e,h-d,R+e-h-d

def tr(x,op,bank):
 h,e,R,d=x
 if op=='A':return h+R-e,h-e,h+R,h
 if op=='K':return h+R-e,d,2*R-d,d
 if op=='D':return h,0,2*R,d
 g,s,L,E,Q,B,delta,dG=bank[op]['sig']
 return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,h*dG

def positive(x,op,bank):
 a,b,c,f=cone(x)
 if op=='A':return b+c,a,2*b+c+f,0
 if op=='K':return a+b,0,2*b+2*c+f,f
 if op=='D':return 0,a+b,c,2*b+c+2*f
 g,s,L,E,Q,B,delta,dG=bank[op]['sig'];U=Q+E-g-dG;V=B-L;W=s-delta
 return ((E-delta)*a+E*b+E*c,(dG-E+delta)*a+(dG-E)*(b+c),
 (g-s+L-dG)*a+(g+2*L-dG)*b+(g+L-dG)*c+L*f,
 (U+V+W)*a+(U+2*V)*b+(U+V)*c+V*f)

def build(rows,r,op,bank):
 rows=list(rows);n=len(rows)
 if op=='D':return [x|(1<<n) for x in rows]+[0],r
 if op in ('A','K'):
  return [x|((i!=r)<<n) for i,x in enumerate(rows)]+[1<<r],n if op=='A' else r
 y=bank[op];G=y['rows'];b=y['root'];mask=((1<<len(G))-1)<<n
 out=[x|mask for x in rows]+[x<<n for x in G]
 out[r]^=1<<(n+b);out[n+b]|=1<<r
 return out,n+b

def saturation(X,bank):
 seen={(1,1,1,1)};q=deque(seen);values={1};ops=['A','K','D']+list(range(len(bank)))
 while q:
  x=q.popleft()
  for op in ops:
   z=tr(x,op,bank)
   check(min(cone(z))>=0 and z[0]>=x[0] and z[2]>=x[2] and z[0]%2==1,'state invariant')
   if z[0]>X:continue
   values.add(z[0])
   if z[2]<=X and z not in seen:seen.add(z);q.append(z)
 return len(seen),len(values),sorted(set(range(1,X+1,2))-values)

def determinant(M):
 total=0;n=len(M)
 for p in itertools.permutations(range(n)):
  v=math.prod(M[i][p[i]] for i in range(n));inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  total+=(-1)**inv*v
 return total

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=1000);args=ap.parse_args()
 raw=json.loads((HERE/'inventory.json').read_text());bank=[]
 for y in raw:
  sig=tuple(y['signature'])+(y['H_deleted_root'],)
  check(dp(y['rows'])[y['root']]==sig==literal(y['rows'],y['root']),'physical gadget profile')
  bank.append({'sig':sig,'rows':y['rows'],'root':y['root']})
 observed=set();small=[];ng=nr=0
 for n in range(1,6):
  for a in graphs(n):
   ps=dp(a);observed.update(ps);ng+=1;nr+=n
   if n<=4:small += [(a,r,state(p)) for r,p in enumerate(ps)]
 check((ng,nr,len(observed))==(1099,5405,65),'complete inventory cardinalities')
 check(observed=={y['sig'] for y in bank},'complete eight-coordinate inventory equality')
 (HERE/'inventory.txt').write_text(''.join(' '.join(map(str,y['sig']))+'\n' for y in bank))
 checks=0
 for a,r,x in small:
  check(min(cone(x))>=0,'input cone')
  ops=['A','K','D']+(list(range(len(bank))) if len(a)<=3 else [])
  for op in ops:
   b,t=build(a,r,op,bank);z=state(dp(b)[t]);pred=tr(x,op,bank)
   check(z==pred,'physical transition');check(cone(z)==positive(x,op,bank),'positive transform');checks+=1
 # Triple collision and subsequent unequal continuations.
 C=[2,12,9,17,7];p=dp(C);collision=[]
 for r in (2,3):
  x=state(p[r]);a=C
  for op in ['K','A']:a,r=build(a,r,op,bank);x=tr(x,op,bank);check(state(dp(a)[r])==x,'collision graph')
  collision.append(x[0])
 check(collision==[55,59],'triple insufficiency hostile')
 # Rank-four minor from actual graph continuations.
 columns=[]
 for word in ('','A','AA','AAD'):
  a=[0];r=0;x=(1,1,1,1)
  for op in word:a,r=build(a,r,op,bank);x=tr(x,op,bank)
  col=[]
  for suffix in ('','A','AA','KA'):
   b=a;t=r;y=x
   for op in suffix:b,t=build(b,t,op,bank);y=tr(y,op,bank)
   check(state(dp(b)[t])==y,'rank witness graph');col.append(y[0])
  columns.append(col)
 M=[list(x) for x in zip(*columns)];check(determinant(M)==-16,'rank minor')
 # Explicit interleaved witnesses, using the frozen audited gadget ordering.
 outputs=[]
 for target,word in [(163,[53,'K','A']),(247,[41,'K','K','A'])]:
  a=[0];r=0;x=(1,1,1,1);steps=[]
  for op in word:
   a,r=build(a,r,op,bank);x=tr(x,op,bank);check(state(dp(a)[r])==x,'explicit witness graph');steps.append({'operator':op,'state':x})
  check(x[0]==target,'explicit witness target')
  outputs.append({'H':target,'rows':a,'root':r,'state':x,'steps':steps})
 # Strongness of the prime witness follows directly, independently checked.
 for w in outputs:
  a=w['rows'];n=len(a)
  for root in range(n):
   reached=1<<root
   while True:
    nxt=reached
    for u in range(n):
     if reached>>u&1:nxt|=a[u]
    if nxt==reached:break
    reached=nxt
   check(reached==(1<<n)-1,'witness not strong')
 # Parameter families and their small claimed prime/7p values.
 for seed,root,coef,off in [([54,48,2,7,12,28],3,88,-13),([54,36,0,7,14,28],3,66,-17)]:
  initial=state(dp(seed)[root])
  for k in range(6):
   x=initial;a=seed;r=root
   for op in ['K']*k+['A']:a,r=build(a,r,op,bank);x=tr(x,op,bank)
   check(state(dp(a)[r])==x and x[0]==coef*2**k+off,'spine family')
 for p in (73,163,691,1039):check(p>1 and all(p%d for d in range(2,math.isqrt(p)+1)),'small prime')
 hostile=state(dp([28,17,26,2,8])[3]);check(hostile==(9,3,9,3),'gcd hostile graph')
 for op in ['A','K','D']+list(range(len(bank))):check(all(v%3==0 for v in tr(hostile,op,bank)),'gcd hostile')
 ns,nv,miss=saturation(args.bound,bank)
 if args.bound>=173:check([x for x in miss if x<=173]==[7,21,67,173],'exact small grammar exclusions')
 (HERE/'explicit_witnesses.json').write_text(json.dumps(outputs,indent=2)+'\n')
 print(json.dumps({'status':'PASS','gadget_graphs':ng,'rooted_gadgets':nr,'eight_entry_signatures':len(bank),'physical_transition_checks':checks,'target_bound':args.bound,'retained_states':ns,'attained_odd_counts':nv,'missing_odd_counts':miss,'rank_minor':M,'rank_determinant':-16,'triple_collision_continuations':collision,'gcd_hostile':hostile,'inventory_sha256':hashlib.sha256((HERE/'inventory.json').read_bytes()).hexdigest(),'scope':'Complete arbitrary-depth combined grammar support at the selected target bound, not all tournaments.'},indent=2))
if __name__=='__main__':main()
