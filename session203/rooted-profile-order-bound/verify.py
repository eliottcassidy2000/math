#!/usr/bin/env python3
"""Independent literal controls for the rooted (B,s) order bound."""
import hashlib,itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
def need(ok,why):
 if not ok:raise RuntimeError(why)
def graphs(n):
 es=list(itertools.combinations(range(n),2))
 for code in range(1<<len(es)):
  a=[0]*n
  for k,(i,j)in enumerate(es):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a
def induced(a,keep):
 return [sum(1<<j for j,v in enumerate(keep)if a[u]>>v&1)for u in keep]
def profile(a):
 n=len(a);N=1<<n;F=N-1;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for i in range(n):end[1<<i][i]=start[1<<i][i]=1
 for S in range(1,N):
  for i in range(n):
   if not(S>>i&1):continue
   rem=S^(1<<i)
   if rem:
    for j in range(n):
     if rem>>j&1:
      if a[j]>>i&1:end[S][i]+=end[rem][j]
      else:start[S][i]+=start[rem][j]
  H[S]=sum(end[S])
 h=H[F];s=start[F][0];L=sum(end[S][0]*H[F^S]for S in range(N))
 return(h,s,L,h+L-s)
def scc(a):
 n=len(a);r=[row|(1<<i)for i,row in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if r[i]>>k&1:r[i]|=r[k]
 groups=[];left=set(range(n))
 while left:
  i=min(left);g=sorted(j for j in left if r[i]>>j&1 and r[j]>>i&1);groups.append(g);left.difference_update(g)
 original_groups=groups[:]
 groups.sort(key=lambda g:-sum(bool(a[g[0]]>>h[0]&1)for h in original_groups if h!=g))
 for i,g in enumerate(groups):
  for h in groups[i+1:]:need(all(a[u]>>v&1 for u in g for v in h),'SCC order')
 return groups
def ham(a,vertices,end=None):
 for w in itertools.permutations(vertices):
  if end is not None and(not w or w[-1]!=end):continue
  if all(a[u]>>v&1 for u,v in zip(w,w[1:])):return w
 raise RuntimeError(('no path',vertices,end))
def test_injection(a,L,s):
 n=len(a);T=induced(a,list(range(1,n)));cs=[[x+1 for x in c]for c in scc(T)]
 if len(cs)==1:
  free=[u for u in range(1,n)if a[u]&1];fixed=();case='strong_deletion'
  need(n-1-len(free)<=s,'outdegree bound')
 else:
  D=cs[-1];free=[u for u in range(1,n)if u not in D]
  w=next(u for u in D if a[u]&1);fixed=ham(a,D,w);case='split_deletion'
  HD=profile(induced(a,D))[0];need(len(D)<=HD<=s,'terminal SCC budget')
 seen=set()
 for bits in range(1<<len(free)):
  S=[v for i,v in enumerate(free)if bits>>i&1]
  first=ham(a,S)+fixed+(0,);rest=[v for v in range(n)if v not in first];second=ham(a,rest)
  need(all(a[u]>>v&1 for path in(first,second)for u,v in zip(path,path[1:])),'illegal cover')
  need(len(first)+len(second)==n and len(set(first+second))==n,'not spanning cover')
  seen.add((first,second))
 need(len(seen)==1<<len(free) and L>=len(seen),'cover injection size')
 return case,len(seen)
def join(parts):
 a=[];offset=0;N=sum(map(len,parts))
 for part in parts:
  m=len(part);later=((1<<N)-1)^((1<<(offset+m))-1)
  a.extend((row<<offset)|later for row in part);offset+=m
 return a
stats={'graphs':0,'positive_roots':0,'strong_positive':0,'weak_positive':0,'deleted_singleton_components':0,'strong_deletion_injections':0,'split_deletion_injections':0,'literal_injected_covers':0};by_order={};boundary=[]
def check(a,inject=False):
 h,s,L,B=profile(a)
 if not s:return
 stats['positive_roots']+=1;need(B>=L>=1,'B positivity')
 cs=scc(a);need(0 in cs[0],'positive root not initial')
 c0=cs[0];p0=profile(induced(a,c0));q=1
 for c in cs[1:]:q*=profile(induced(a,c))[0]
 need((h,s,L,B)==tuple(q*x for x in p0),'component product h,s,L,B')
 keep=sorted(c0+[v for c in cs[1:]if len(c)>1 for v in c]);b=induced(a,keep);pp=profile(b)
 need(pp==(h,s,L,B),'normalization changed profile')
 need(len(b)<=s+1+B.bit_length()-1,'normalized order bound')
 stats['deleted_singleton_components']+=len(a)-len(b)
 if len(cs)==1:
  stats['strong_positive']+=1;need(L*(1<<s)>=1<<(len(a)-1),'strong exponential bound')
  if inject and len(a)>1:
   case,num=test_injection(a,L,s);stats[case+'_injections']+=1;stats['literal_injected_covers']+=num
 else:stats['weak_positive']+=1
 if len(a)<=2:boundary.append({'rows':a,'profile_h_s_L_B':[h,s,L,B],'normalized_rows':b})
for n in range(1,6):
 count=0
 for a in graphs(n):stats['graphs']+=1;count+=1;check(a,True)
 by_order[n]=count
tri=[2,4,1];single=[0]
examples=[join([single]*12),join([tri,single,tri,single,tri]),join([single,tri,single,tri,single,tri]),join([single,single,tri,single,tri])]
example_rows=[]
for a in examples:
 check(a);h,s,L,B=profile(a);cs=scc(a);keep=sorted(cs[0]+[v for c in cs[1:]if len(c)>1 for v in c]);example_rows.append({'rows':a,'order':len(a),'profile_h_s_L_B':[h,s,L,B],'normalized_rows':induced(a,keep),'order_cap':s+B.bit_length()})
need(profile([0])==(1,1,1,1),'singleton convention')
need(profile([2,0])==(1,1,1,1),'source edge convention')
need(profile([0,1])[1]==0,'sink edge excluded')
payload={'scope':'All labeled root-zero tournaments through order5 plus four explicit larger weak-tail graphs; exact controls, not an occupancy proof.','status':'PASS','labeled_universe_by_order':by_order,'statistics':stats,'small_boundaries':boundary,'weak_tail_examples':example_rows}
data=json.dumps(payload,sort_keys=True,indent=2)+'\n';print(data,end='')
