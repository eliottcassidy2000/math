#!/usr/bin/env python3
import json,random,time,hashlib
from pathlib import Path
BASE=Path(__file__).resolve().parent
n=11;FULL=(1<<n)-1;rng=random.Random(7307311)
squares={x*x%n for x in range(1,n)}
T=tuple(sum(1<<j for j in range(n) if (j-i)%n in squares) for i in range(n))
def need(x,msg):
 if not x:raise RuntimeError(msg)
def fas_dp(out):
 d=[0]*(1<<n)
 for S in range(1,1<<n):
  z=S;best=999
  while z:
   v=(z&-z).bit_length()-1;z&=z-1;R=S^(1<<v);c=d[R]+(out[v]&R).bit_count()
   if c<best:best=c
  d[S]=best
 return d
def sample_order(out,d,slack=0):
 S=FULL;rev=[];remaining=slack
 while S:
  opts=[];z=S
  while z:
   v=(z&-z).bit_length()-1;z&=z-1;R=S^(1<<v);extra=d[R]+(out[v]&R).bit_count()-d[S]
   if extra<=remaining:opts.append((v,extra))
  v,extra=rng.choice(opts);remaining-=extra;rev.append(v);S^=1<<v
 return tuple(reversed(rev))
def defects(out,order):
 pos={v:i for i,v in enumerate(order)}
 return [(i,j) for i in range(n) for j in range(i+1,n) if bool(out[i]>>j&1)!=(pos[i]<pos[j])]
def beta(edges):
 p=list(range(n))
 def root(x):
  while p[x]!=x:p[x]=p[p[x]];x=p[x]
  return x
 for a,b in edges:
  x,y=root(a),root(b)
  if x!=y:p[x]=y
 return len(edges)-n+len({root(x) for x in range(n)})
def max_path_dp(out,order):
 pos={v:i for i,v in enumerate(order)}
 w=[[int(pos[i]>pos[j]) for j in range(n)] for i in range(n)]
 d=[[-999]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=0
 for S in range(1,1<<n):
  z=S
  while z:
   v=(z&-z).bit_length()-1;z&=z-1;cur=d[S][v]
   if cur<0:continue
   opts=out[v] & (FULL^S)
   while opts:
    x=(opts&-opts).bit_length()-1;opts&=opts-1;R=S|(1<<x);value=cur+w[v][x]
    if value>d[R][x]:d[R][x]=value
 return d,w
def sample_path(out,d,w):
 best=max(d[FULL]);v=rng.choice([v for v in range(n) if d[FULL][v]==best]);rev=[v];S=FULL
 while S!=(1<<v):
  R=S^(1<<v);opts=[u for u in range(n) if R>>u&1 and out[u]>>v&1 and d[R][u]+w[u][v]==d[S][v]]
  v=rng.choice(opts);rev.append(v);S=R
 return tuple(reversed(rev)),best
def tau(out,P):
 need(all(out[a]>>b&1 for a,b in zip(P,P[1:])),'illegal path')
 pe={tuple(sorted(x)) for x in zip(P,P[1:])};res=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   arc=bool(out[i]>>j&1)
   if (i,j) not in pe:arc=not arc
   a,b=(i,j) if arc else (j,i);res[a]|=1<<b
 return tuple(res)

rng=random.Random(73073311)
input_receipt=json.loads((BASE/'input_two_step_receipt.json').read_text())
best=input_receipt['best_literal_candidate'];start=tuple(best['final_rows']);ds=fas_dp(start)
need(ds[-1]==13,'input best state cost')
stats={'seed':73073311,'budget':{'optimal_rank_draws':32,'maximum_overlap_paths_per_rank':4,'final_rank_draws_if_FAS_at_most10':32,'final_slack_limit':2},'rank_draws':0,'path_draws':0,'distinct_targets':0,'final_rank_draws':0,'minimum_FAS':999,'minimum_observed_beta':None}
seen=set();found=None;best_seen=None
for a in range(32):
 stats['rank_draws']+=1;r=sample_order(start,ds);pd,pw=max_path_dp(start,r)
 for j in range(4):
  stats['path_draws']+=1;R,overlap=sample_path(start,pd,pw);F=tau(start,R)
  if F in seen:continue
  seen.add(F);stats['distinct_targets']+=1;df=fas_dp(F)
  if df[-1]<stats['minimum_FAS']:
   stats['minimum_FAS']=df[-1];best_seen={'R':R,'rows':F,'FAS':df[-1],'overlap':overlap,'source_rank':r}
  if df[-1]>10:continue
  for h in range(32):
   stats['final_rank_draws']+=1;slack=0 if h<16 else min(2,10-df[-1]);rf=sample_order(F,df,slack);D=defects(F,rf);be=beta(D)
   stats['minimum_observed_beta']=be if stats['minimum_observed_beta'] is None else min(be,stats['minimum_observed_beta'])
   if be==0:
    found={'P':best['P'],'Q':best['Q'],'R':R,'forest_rows':F,'forest_rank':rf,'forest_edges':D,'forest_FAS':df[-1],'forest_beta':be,'intermediate_after_P':best['intermediate_rows'],'intermediate_after_PQ':best['final_rows'],'forward_complements':[R,best['Q'],best['P']],'whole_complement_clock':3,'remarking_note':'Actual named paths are legal at their stages. No elementary relocation clock claim.'};break
  if found:break
 if found:break
stats['status']='FOUND_LITERAL_THREE_COMPLEMENT_FOREST_CERTIFICATE' if found else 'NO_WITNESS_IN_DECLARED_BEST_STATE_BATCH; THREE_COMPLEMENT_EXISTENCE_NOT_DECIDED'
if found:stats['certificate']=found
if best_seen:stats['best_target']=best_seen
print(json.dumps(stats,indent=2,sort_keys=True))
