#!/usr/bin/env python3
"""Hostile controls for the all-order full binary-return theorem. No scalar coverage census."""
import hashlib,itertools,json,math
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


def ready(x):
 h,e,R,d=x
 return d>=3 and e==d and abs(R-h-d)==1

def gap(x):
 h,e,R,d=x;c=R-h;f=R+e-h-d
 return h>=9 and d>=3 and all(t==0 or t>=3 for t in (e,c,f))

def main():
 raw=json.loads((HERE/'inventory67.json').read_text());bank=[]
 for z in raw:
  sig=tuple(z['signature'])+(z['H_deleted_root'],)
  check(dp(z['rows'])[z['root']]==sig==literal(z['rows'],z['root']),'gadget profile')
  bank.append({'sig':sig,'rows':z['rows'],'root':z['root']})
 check(len(bank)==67,'bank size')
 profiles=set();graph_count=rooted_count=0
 for n in range(1,7):
  for a in graphs(n):
   ps=dp(a);graph_count+=1;rooted_count+=n
   profiles.update(state(p) for p in ps)
 check((graph_count,rooted_count)==(33867,202013),'control enumeration size')
 extras=[([52,57,58,17,32,8],3),([2,12,9,17,7],1)]
 for a,r in extras:profiles.add(state(dp(a)[r]))
 ops=['K','D']+list(range(len(bank)));gap_inputs=ready_inputs=endpoint_tests=gap_checks=entry_checks=return_checks=0
 seeds=set()
 for x in profiles:
  h,e,R,d=x
  if h>=9:
   check(h-e==0 or h-e>=3,'endpoint gap');endpoint_tests+=1
  if ready(x):
   seeds.add(x);ready_inputs+=1
   check(h>=3*d,'ready weak-deletion ratio')
   for op in ops:
    z=tr(x,op,bank)
    if op=='K':check(ready(z) and z[3]==d and z[2]-z[0]-z[3]==R-h-d,'K ready return')
    else:check(gap(z),'entry into gap')
    entry_checks+=1
  if gap(x):
   gap_inputs+=1
   for op in ops:
    z=tr(x,op,bank);check(gap(z),'gap invariant');gap_checks+=1
   for op in range(len(bank)):
    z=tr(x,op,bank)
    for suffix in ('','K','DK','DDK'):
     y=z
     for s in suffix:y=tr(y,s,bank)
     check(not ready(y),'last-root-switch short return');return_checks+=1
 # Exact saturation of marked histories up to R bound, not scalar spectrum coverage.
 bound=3000;q=deque((x,False) for x in seeds);seen=set(q);hist_checks=0;dirty_states=0;neutral_transitions=0
 while q:
  x,dirty=q.popleft()
  for op in ops:
   z=tr(x,op,bank);newdirty=dirty or op!='K';hist_checks+=1
   check(z[2]>=x[2] and z[0]>=x[0],'history monotonicity')
   if z[0]==x[0]:neutral_transitions+=1
   if newdirty:
    check(gap(z),'history gap invariant')
    check(not ready(z),'non-K history returned')
   else:check(ready(z),'pure K failed to return')
   item=(z,newdirty)
   if z[2]<=bound and item not in seen:seen.add(item);q.append(item)
 dirty_states=sum(d for x,d in seen)
 print(json.dumps({'status':'PASS','labelled_graphs_through6':graph_count,'rooted_profiles_through6':rooted_count,
  'distinct_control_states':len(profiles),'gadget_profiles':len(bank),'endpoint_gap_inputs':endpoint_tests,
  'gap_inputs':gap_inputs,'mature_ready_inputs':ready_inputs,'gap_invariance_checks':gap_checks,
  'mature_entry_checks':entry_checks,'last_switch_short_suffix_checks':return_checks,
  'history_R_bound':bound,'history_seed_states':sorted(seeds),'retained_history_states':len(seen),
  'retained_dirty_states':dirty_states,'history_transitions_tested':hist_checks,'neutral_count_transitions':neutral_transitions,
  'inventory_sha256':hashlib.sha256((HERE/'inventory67.json').read_bytes()).hexdigest(),
  'scope':'Exhaustive finite hostile controls for the structural proof; no scalar spectrum census and no inference of universal truth from the finite bound.'},indent=2))
if __name__=='__main__':main()
