#!/usr/bin/env python3
import sys as _sys
if "--fail" in _sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Independent auditor: digit-pair Q, pumping certificates, partition refinement."""
from collections import defaultdict
from itertools import product,combinations,permutations
import json

def check(p,msg):
 if not p: raise RuntimeError(msg)

def q_direct(n):
 # Evaluate the digit-pair identity; never call the producer's Q recurrence.
 ones=[i for i in range(n.bit_length()) if (n>>i)&1]
 Q=sum(1<<(j-i-1) for i in ones for j in ones if i<j)+sum(j==i+1 for i in ones for j in ones)
 check(Q%2==0,('Q not even',n)); return Q//2

def witness(M,r,z):
 u=M;k=1
 while u%2==0: u//=2;k+=1
 b=max(k,2)
 while ((1<<b)-1)%u: b+=1
 c=(1<<(b-1))-1
 # t in 1..M avoids relying on t=0 and supports the large-witness statement.
 t=next(t for t in range(1,M+1) if (q_direct(r)+t*c-z)%M==0)
 A=r.bit_length()+k
 n=r+sum(((1<<b)-1)<<(A+j*(b+k)) for j in range(t))
 return n

def graph(M):
 states=list(product(range(2*M),range(M)))
 def tr(s,d):
  r,z=s
  return ((2*r+d)%(2*M),(z+d*((r+1)//2))%M)
 def output(s):return (s[0]-s[1])%M
 # Reachability as graph closure, then actual equivalence partition refinement.
 reached={(0,0)}
 while True:
  enlarged=reached|{tr(s,d) for s in reached for d in (0,1)}
  if enlarged==reached:break
  reached=enlarged
 groups={s:output(s) for s in states};rounds=0; counts=[len(set(groups.values()))]
 while True:
  sig={s:(output(s),groups[tr(s,0)],groups[tr(s,1)]) for s in states}
  ids={v:i for i,v in enumerate(sorted(set(sig.values())))}
  refined={s:ids[sig[s]] for s in states}; rounds+=1
  count=len(ids);counts.append(count)
  if count==len(set(groups.values())):break
  groups=refined
 check(len(reached)==2*M*M,('reach',M))
 check(count==(1 if M==1 else 2*M*M),('minimize',M,count))
 check(counts[1]==count,('longer suffix needed',M,counts))
 return dict(M=M,reachable=len(reached),minimal=count,partition_counts=counts)

mods=[1,2,3,4,5,6,8,9,10,12,15,16,18,21,25,27,30,32,35]
automata=[graph(M) for M in mods]
seeds=[(1,2,0,1),(5,6,2,3),(9,11,3,5)]
witnesses=0;maxbits=0;joint=[]
for M in range(1,18):
 observed=set(); parity_observed={0:set(),1:set()}
 for r,z in product(range(2*M),range(M)):
  n=witness(M,r,z);q=q_direct(n);witnesses+=1;maxbits=max(maxbits,n.bit_length())
  check(n%(2*M)==r and q%M==z,('certificate',M,r,z))
  vals=[]
  for h,R,e,d in seeds:
   H=h+R*n-2*d*q-e*(n%2)
   check(H%2==1,('not odd',H))
   vals.append(((H-1)//2)%M)
  observed.add(tuple(vals));parity_observed[n%2].add((vals[0],vals[2]))
 predicted={tuple(((h-1+(R-e)*b)//2+R*j-d*z)%M for h,R,e,d in seeds) for b,j,z in product(range(2),range(M),range(M))}
 check(observed==predicted,('joint equality',M))
 check(all(len(s)==M*M for s in parity_observed.values()),('fixed parity unimodular',M))
 pair={(a,b) for a,b,c in observed}
 lines={(a,(2+3*a-b)%M) for a,b in product(range(M),range(2))}
 check(pair==lines,('two exact lines',M))
 joint.append(dict(M=M,joint_triples=len(observed),root_and_rank_one=len(pair),root_and_unimodular=M*M))
# Recurrence independently checked only after generating digit-pair values.
for n in range(10000):
 q=q_direct(n)
 check(q_direct(2*n)==q,('even recurrence',n))
 check(q_direct(2*n+1)==q+(n+1)//2,('odd recurrence',n))
def brute(rows, vertices=None, ending=None):
 vs=list(range(len(rows))) if vertices is None else vertices
 return sum(all(rows[a]>>b&1 for a,b in zip(P,P[1:])) and (ending is None or P[-1]==ending) for P in permutations(vs))
def append(rows,b):
 n=len(rows);out=rows[:]
 for v in range(n):
  if v!=0 or not b:out[v]|=1<<n
 out.append(1 if b else 0)
 return out
literal_checks=0;found=[]
for target,n in [((1,2,0,1),2),((5,6,2,3),4),((9,11,3,5),5)]:
 h,R,e,d=target;pairs=list(combinations(range(n),2));base=None
 for mask in range(1<<len(pairs)):
  rows=[0]*n
  for j,(a,b) in enumerate(pairs):
   u,v=(a,b) if mask>>j&1 else (b,a);rows[u]|=1<<v
  if brute(rows)!=h or brute(rows,ending=0)!=e or brute(rows,list(range(1,n)))!=d:continue
  if brute(append(rows,1))-h+e==R:base=rows;break
 check(base is not None,('physical seed missing',target))
 found.append(dict(seed=target,rows=base))
 for m in range(32 if n==2 else 8):
  rows=base[:]
  for i in range(m.bit_length()):rows=append(rows,m>>i&1)
  H=brute(rows)
  check(H==h+R*m-2*d*q_direct(m)-e*(m%2),('literal formula',target,m))
  for M in (1,2,3,6):
   check(((H-1)//2)%M==((h-1+(R-e)*(m%2))//2+R*(m//2)-d*q_direct(m))%M,('normalization',target,m,M))
  literal_checks+=1
print(json.dumps(dict(status='PASS',physical_seeds=found,literal_permutation_checks=literal_checks,automata=automata,actual_pumping_witnesses=witnesses,max_witness_bits=maxbits,joint_images=joint,digit_recurrence_checks=20000),indent=2))
