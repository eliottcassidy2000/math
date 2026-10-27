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
# Fully declared finite search budget. Randomness selects DP ties, not optimal costs.
budget={'first_optimal_rank_draws':64,'first_maximum_overlap_paths_per_rank':8,'second_optimal_rank_draws_per_distinct_first_state':8,'second_maximum_overlap_paths_per_rank':4,'final_rank_draws_per_second_state':32,'max_final_rank_slack':2}
stats={'seed':7307311,'budget':budget,'first_rank_draws':0,'first_path_draws':0,'distinct_first_states':0,'second_rank_draws':0,'second_path_draws':0,'distinct_second_states':0,'final_rank_draws':0,'first_min_fas':999,'second_min_fas':999,'minimum_observed_final_beta':None}
seen1=set();seen2=set();dt=fas_dp(T);need(dt[-1]==20,'P11 FAS')
found=None
best_candidate=None
def canonical_order(out,d):
 S=FULL;rev=[]
 while S:
  v=next(v for v in range(n) if S>>v&1 and d[S^(1<<v)]+(out[v]&(S^(1<<v))).bit_count()==d[S])
  rev.append(v);S^=1<<v
 return tuple(reversed(rev))
for a in range(budget['first_optimal_rank_draws']):
 stats['first_rank_draws']+=1;r=sample_order(T,dt);pd,pw=max_path_dp(T,r)
 for b in range(budget['first_maximum_overlap_paths_per_rank']):
  stats['first_path_draws']+=1;P,s1=sample_path(T,pd,pw);A=tau(T,P)
  if A in seen1:continue
  seen1.add(A);stats['distinct_first_states']+=1;da=fas_dp(A);stats['first_min_fas']=min(stats['first_min_fas'],da[-1])
  for c in range(budget['second_optimal_rank_draws_per_distinct_first_state']):
   stats['second_rank_draws']+=1;r2=sample_order(A,da);qd,qw=max_path_dp(A,r2)
   for e in range(budget['second_maximum_overlap_paths_per_rank']):
    stats['second_path_draws']+=1;Q,s2=sample_path(A,qd,qw);F=tau(A,Q)
    if F in seen2:continue
    seen2.add(F);stats['distinct_second_states']+=1;df=fas_dp(F)
    if df[-1]<stats['second_min_fas']:
     rr=canonical_order(F,df);ee=defects(F,rr)
     best_candidate={'P':P,'Q':Q,'intermediate_rows':A,'intermediate_FAS':da[-1],'final_rows':F,'final_FAS':df[-1],'final_rank':rr,'final_edges':ee,'final_rank_beta':beta(ee),'whole_complement_clock':2}
    stats['second_min_fas']=min(stats['second_min_fas'],df[-1])
    if df[-1]>n-1:continue
    for h in range(budget['final_rank_draws_per_second_state']):
     stats['final_rank_draws']+=1;slack=min(budget['max_final_rank_slack'],n-1-df[-1]) if h>=16 else 0
     rf=sample_order(F,df,slack);D=defects(F,rf);be=beta(D);stats['minimum_observed_final_beta']=be if stats['minimum_observed_final_beta'] is None else min(stats['minimum_observed_final_beta'],be)
     if be==0:
      found={'P11_rows':T,'P':P,'first_rank':r,'first_path_overlap':s1,'intermediate_rows':A,'intermediate_FAS':da[-1],'Q':Q,'second_rank':r2,'second_path_overlap':s2,'forest_rows':F,'forest_rank':rf,'forest_edges':D,'forest_FAS':df[-1],'forest_beta':be,'forward_complements':[Q,P],'whole_complement_clock':2,'required_remarking':'after Q complement, replace active Q by literal legal P; no elementary relocation count claimed'}
      break
    if found:break
   if found:break
  if found:break
 if found:break
stats['status']='FOUND_LITERAL_TWO_COMPLEMENT_FOREST_CERTIFICATE' if found else 'NO_WITNESS_IN_DECLARED_FINITE_SEARCH; TWO_COMPLEMENT_EXISTENCE_OPEN'
if found:stats['certificate']=found
if best_candidate:stats['best_literal_candidate']=best_candidate
(BASE/'search_result.json').write_text(json.dumps(stats,indent=2,sort_keys=True)+'\n')
print(json.dumps(stats,indent=2,sort_keys=True))
