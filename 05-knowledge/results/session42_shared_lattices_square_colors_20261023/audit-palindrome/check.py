#!/usr/bin/env python3
import json,math,itertools,functools
from pathlib import Path
checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def square(x):return x>0 and math.isqrt(x)**2==x
def Aset(N,d):return [v for v in range(1,N+1) if (v-1)//d%2==0]
def quotient(N,d):
 A=Aset(N,d);adj=[tuple(j for j,b in enumerate(A) if a!=b and square(a+b+d)) for a in A];marks=[i for i,a in enumerate(A) if square(2*a+d)]
 return A,adj,marks
cert=[];nontrivial=[]
for N in range(2,48,2):
 for d in range(1,N//2+1):
  if N%(2*d):continue
  side=lambda v:(v-1)//d%2
  full={a:[b for b in range(1,N+1) if side(a)!=side(b) and square(a+b)] for a in range(1,N+1)}
  isolated=[a for a in full if not full[a]];leaves=[a for a in full if len(full[a])==1]
  reached={1};todo=[1]
  for a in todo:
   for b in full[a]:
    if b not in reached:reached.add(b);todo.append(b)
  if isolated:record={'kind':'isolated','witness':isolated}
  elif len(leaves)>2:record={'kind':'leaves','witness':leaves}
  elif len(reached)<N:record={'kind':'disconnected','witness':sorted(reached)}
  else:
   A,adj,marks=quotient(N,d);nontrivial.append((N,d));allmask=(1<<len(A))-1
   @functools.lru_cache(None)
   def count(last,mask):
    if mask==allmask:return 1
    return sum(count(j,mask|(1<<j)) for j in adj[last] if not (mask>>j)&1)
   total=sum(count(i,1<<i) for i in marks)
   need(total==0,('earlier palindrome',N,d,total))
   record={'kind':'quotient_exact_dp','representatives':A,'marked_representatives':[A[i] for i in marks],'paths_starting_at_mark':total,'memoized_states':count.cache_info().currsize}
  cert.append({'N':N,'d':d,**record})
need(nontrivial==[(36,6),(40,4)],'only two nontrivial candidates')
# A smaller structural certificate supersedes DP for N36.
A36,adj36,marks36=quotient(36,6);cut={1,14};seen=set(cut);components=[]
for a in A36:
 if a in seen:continue
 todo=[a];seen.add(a)
 for x in todo:
  for j in adj36[A36.index(x)]:
   y=A36[j]
   if y not in seen:seen.add(y);todo.append(y)
 components.append(sorted(todo))
need(sorted(components)==sorted([[2,4,6,13,15,17,26,28,30],[3,16,27],[5,18,25],[29]]),'structural separator certificate')
need(len(components)>len(cut)+1,'Hamilton separator obstruction')
# Positive witness is literal, independent of mutable root task files.
p=(3,46,18,31,5,44,37,27,22,42,7,9,40,41,8,28,21,43,6,10,39,25,24,12,4,32,17,47,2,14,35,29,20,16,33,48,1,15,34,30,19,45,36,13,23,26,38,11)
N=48;d=8
need(sorted(p)==list(range(1,49)),'witness permutation')
s=tuple(a+b for a,b in zip(p,p[1:]));need(all(map(square,s)) and s==s[::-1],'palindromic square word')
A,adj,marks=quotient(N,d);tau=lambda v:v+d if (v-1)//d%2==0 else v-d
need(tuple(map(tau,p))==p[::-1],'reversal is block twin')
rep=lambda v:min(v,tau(v))
q=tuple(rep(v) for v in p[:N//2]);witness_q=q;need(sorted(q)==A,'quotient permutation')
need(all(square(a+b+d) for a,b in zip(q,q[1:])) and square(2*q[-1]+d),'quotient witness')
# Exact reduction checked independently for every pair-order and start side through N8.
reduction_checks=0
for N in range(2,10,2):
 m=N//2
 for d in range(1,m+1):
  if m%d:continue
  A=Aset(N,d)
  for q in itertools.permutations(A):
   predicted=all(square(a+b+d) for a,b in zip(q,q[1:])) and square(2*q[-1]+d)
   for start in (0,1):
    half=tuple(a+((i+start)%2)*d for i,a in enumerate(q))
    t=lambda v:v+d if (v-1)//d%2==0 else v-d
    path=half+tuple(t(v) for v in half[::-1]);sums=tuple(a+b for a,b in zip(path,path[1:]))
    actual=all(map(square,sums))
    need(sorted(path)==list(range(1,N+1)) and sums==sums[::-1],'all pair-order constructions')
    need(actual==predicted,'quotient reduction')
    reduction_checks+=1
print(json.dumps({'status':'PASS','checks':checks,'first_nontrivial_palindromic_square_sum_order':48,'lower_order_certificates':cert,'nontrivial_lower_candidates':nontrivial,'N36_separator':{'cut':sorted(cut),'components':components},'finite_reduction_controls':reduction_checks,'witness':{'N':48,'d':8,'path':p,'quotient_path':witness_q,'central_representative':witness_q[-1],'square_clock':[math.isqrt(a) for a in s]}},indent=2))
