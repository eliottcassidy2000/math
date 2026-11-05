#!/usr/bin/env python3
"""Exact finite controls for the all-order first symmetric 7p theorem.
No third-party packages, no assert statements, no external data.
"""
import itertools as it
import json
import math
from collections import Counter
from pathlib import Path

checks=0

def check(p,msg):
 global checks
 checks+=1
 if not p:raise RuntimeError(msg)

def rows_from_mask(n,mask):
 r=[0]*n
 for k,(i,j) in enumerate(it.combinations(range(n),2)):
  if mask>>k&1:r[i]|=1<<j
  else:r[j]|=1<<i
 return r

def valid(r,p):return all(r[a]>>b&1 for a,b in zip(p,p[1:]))

def allpaths(r):
 n=len(r);paths=[]
 def go(p,unused):
  if not unused:paths.append(tuple(p));return
  for v in range(n):
   if unused>>v&1 and (not p or r[p[-1]]>>v&1):go(p+[v],unused^(1<<v))
 go([], (1<<n)-1)
 return paths

def hp_table(r):
 n=len(r);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for m in range(1,1<<n):
  for v in range(n):
   if not(m>>v&1):continue
   prev=m^(1<<v)
   dp[m][v]=1 if not prev else sum(dp[prev][u] for u in range(n) if prev>>u&1 and r[u]>>v&1)
  h[m]=sum(dp[m])
 return h

def hp(r):return hp_table(r)[-1]
def f2(r):
 h=hp_table(r);full=(1<<len(r))-1
 return sum(h[m]*h[full^m] for m in range(1,full))
def induce(r,vs):return [sum(1<<j for j,b in enumerate(vs) if r[a]>>b&1)for a in vs]
def strong(r):
 n=len(r)
 for root in range(n):
  seen=1<<root;old=-1
  while old!=seen:
   old=seen
   for v in range(n):
    if seen>>v&1:seen|=r[v]
  if seen!=(1<<n)-1:return False
 return True

def circ(q,ds):return [sum(1<<((i+d)%q) for d in ds)for i in range(q)]
def skewsets(q):
 for b in it.product([0,1],repeat=(q-1)//2):yield tuple(i if x else q-i for i,x in enumerate(b,1))
def ear(r):
 n=len(r);return [x|1<<n for x in r]+[1<<(n+1),(1<<n)-1]
def substitute(q,root,r):
 groups=[];pos=0
 for i in range(len(q)):
  sz=len(r) if i==root else 1
  groups.append(list(range(pos,pos+sz)));pos+=sz
 out=[0]*pos
 for a in range(len(q)):
  if a==root:
   for i,u in enumerate(groups[a]):
    for j,v in enumerate(groups[a]):
     if r[i]>>j&1:out[u]|=1<<v
  for b in range(len(q)):
   if q[a]>>b&1:
    for u in groups[a]:
     for v in groups[b]:out[u]|=1<<v
 return out

def merge(r,p,q):
 p=list(p);q=list(q);o=[]
 while p and q:
  o.append(p.pop(0) if r[p[0]]>>q[0]&1 else q.pop(0))
 return tuple(o+p+q)

def aut(r):
 n=len(r);cells={}
 for i,row in enumerate(r):cells.setdefault(row.bit_count(),[]).append(i)
 cells=list(cells.values());ans=[]
 for z in it.product(*(list(it.permutations(c))for c in cells)):
  perm=[None]*n
  for c,d in zip(cells,z):
   for u,v in zip(c,d):perm[u]=v
  if all(((r[u]>>v)&1)==((r[perm[u]]>>perm[v])&1)for u in range(n)for v in range(n)):
   ans.append(tuple(perm))
 return ans

def canonical(r):
 n=len(r);best=None
 # Only used n<=4, or degree-partition-canonical n10 with seven equal scores.
 cells={}
 for i,row in enumerate(r):cells.setdefault(row.bit_count(),[]).append(i)
 cells=[cells[s] for s in sorted(cells)]
 for chunks in it.product(*(list(it.permutations(c))for c in cells)):
  pp=sum((list(c)for c in chunks),[])
  val=tuple(sum(1<<j for j in range(n) if r[pp[i]]>>pp[j]&1) for i in range(n))
  if best is None or val<best:best=val
 return best

def reverse(r):
 full=(1<<len(r))-1
 return [(full^(1<<i))^row for i,row in enumerate(r)]

# Full labelled controls: independent permutation/path enumeration, shuffle,
# induced marked strong4, and strict proper strong inclusion.
labelled=0;shuffles=0;marked=0;strong5_hist=Counter()
for n in range(1,6):
 for mask in range(1<<(n*(n-1)//2)):
  r=rows_from_mask(n,mask);labelled+=1;h=hp_table(r)
  paths=allpaths(r)
  check(len(paths)==h[-1],('hp',n,mask))
  for m in range(1,(1<<n)-1):
   aa=[v for v in range(n)if m>>v&1];bb=[v for v in range(n)if not(m>>v&1)]
   ap=[tuple(aa[i]for i in p)for p in allpaths(induce(r,aa))]
   bp=[tuple(bb[i]for i in p)for p in allpaths(induce(r,bb))]
   seen=set()
   for p in ap:
    for q in bp:
     z=merge(r,p,q);shuffles+=1
     check(valid(r,z) and len(z)==n,('shuffle directed',r,p,q,z))
     check(tuple(v for v in z if v in aa)==p and tuple(v for v in z if v in bb)==q,'shuffle restrictions')
     check(z not in seen,'shuffle injectivity');seen.add(z)
   if strong(r) and len(aa)>1 and strong(induce(r,aa)):
    check(h[m]<h[-1],'strict strong induced monotonicity')
  if n==5 and strong(r):strong5_hist[h[-1]]+=1
  if n>=4 and strong(r):
   for root in range(n):
    check(any(root in vs and strong(induce(r,vs))for vs in it.combinations(range(n),4)),'marked strong4')
    marked+=1

check(dict(strong5_hist)=={9:240,11:120,13:120,15:64},'complete strong5 histogram')

# Complete synchronized pairs of order-three automorphism orbits.
c3_pair_hist=Counter()
for direction in [1,2]:
 for mask in range(8):
  r=circ(3,[1])+[x<<3 for x in circ(3,[direction])]
  for i in range(3):
   for j in range(3):
    if mask>>((j-i)%3)&1:r[i]|=1<<(3+j)
    else:r[3+j]|=1<<i
  value=hp(r);c3_pair_hist[value]+=1
  check(value==(9 if mask in [0,7] else 45),'synchronized C3 pair')
  check(len(allpaths(r))==value,'independent C3 pair DFS')
check(dict(c3_pair_hist)=={9:4,45:12},'complete synchronized C3 histogram')

orbit_rows=[]
for q in [5,7]:
 vals=set()
 for ds in skewsets(q):
  r=circ(q,ds);triple=(hp(r),f2(r),hp(ear(r)))
  check(triple[2]==3*triple[0]+triple[1],'ear identity')
  vals.add(triple);orbit_rows.append({'q':q,'steps':list(ds),'H':triple[0],'F2':triple[1],'ear_H':triple[2]})
 check(vals==({(15,90,135)} if q==5 else {(175,1638,2163),(189,1638,2205)}),('orbit table',q,vals))

# An elementary, source-independent replacement for the cubic orbit floor:
# p>=17 has >=8 constant-step cycles. At p=11,13 enumerate all skew sets and
# exhibit eight Hamilton paths with first vertex 0; translates are disjoint.
small_prime_cert=[]
for q in [11,13]:
 for ds in skewsets(q):
  r=circ(q,ds);found=[]
  def dfs(p,unused):
   if len(found)>=8:return
   if not unused:found.append(tuple(p));return
   v=p[-1]
   for w in range(q):
    if unused>>w&1 and r[v]>>w&1:dfs(p+[w],unused^(1<<w))
  dfs([0],((1<<q)-1)^1)
  check(len(found)==8 and len(set(found))==8 and all(p[0]==0 and len(p)==q and len(set(p))==q and valid(r,p)for p in found),'prime symmetry certificate')
  small_prime_cert.append({'q':q,'steps':list(ds),'anchored_paths':[list(p)for p in found]})

# Every labelled marked strong4 presentation, both 7-orbit types.
# Count the universal run patterns independently, rather than inferring
# profile coefficients from the two tested module types.
def run_patterns(q,k):
 words=set(it.permutations([0]*k+[1,2,3]))
 return [p for p in words if all(a!=b and q[a]>>b&1 for a,b in zip(p,p[1:]))]
q4_count=0;sub_values=set();coefs=set();strong4_classes=set()
for mask in range(64):
 q=rows_from_mask(4,mask)
 if not strong(q):continue
 q4_count+=1;strong4_classes.add(canonical(q));check(hp(q)==5,'strong4 H')
 check(len(run_patterns(q,1))==5 and len(run_patterns(q,2))in[2,3] and len(run_patterns(q,3))==0,'universal rooted4 run patterns')
 for ds in skewsets(7):
  r=circ(7,ds);t=substitute(q,0,r);v=hp(t);sub_values.add(v)
  delta=v-5*hp(r);check(delta%f2(r)==0,'rooted4 linear profile')
  coefs.add(delta//f2(r));check(delta//f2(r)==len(run_patterns(q,2)),'rooted4 coefficient')
check(q4_count==24 and len(strong4_classes)==1,'complete strong4 inventory')
check(sub_values=={4151,4221,5789,5859}and coefs=={2,3},'rooted4 substitution values')

q=[2,8,3,5];c=circ(7,[1,2,3]);witness=substitute(q,0,c)
check(witness==[142,156,184,240,225,195,135,512,255,383],'witness rows')
check(strong(witness)and hp(witness)==4151,'witness H and strong')
check(len(allpaths(witness))==4151,'independent DFS witness')
check(len(aut(witness))==7,'witness automorphisms')
check(all(593%d for d in range(2,math.isqrt(593)+1)),'593 prime')
roots=[hp(substitute(q,i,c))for i in range(4)]
check(roots==[4151,5789,4151,5789],'four marked positions')
w2=substitute(q,2,c)
check(canonical(witness)!=canonical(w2),'two equality cores distinct')
check(canonical(reverse(witness))==canonical(w2),'two equality cores converse')
check(canonical(reverse(witness))!=canonical(witness),'equality core not self-converse')
check(len(aut(q))==1,'strong4 rigid')

receipt={'status':'PASS','checks':checks,'all_labelled_orders':'1..5','labelled_tournaments':labelled,
 'shuffle_pairs':shuffles,'marked_strong4_controls':marked,'q5_q7_presentations':len(orbit_rows),'universal_q4_run_patterns_verified':True,'strong5_H_histogram':dict(strong5_hist),'synchronized_C3_pair_H_histogram':dict(c3_pair_hist),
 'q11_q13_anchored_certificates':len(small_prime_cert),'q4_labelled_strong_presentations':q4_count,
 'q4_q7_substitutions':q4_count*8,'rooted4_H_values':sorted(sub_values),'witness':witness,
 'witness_H':4151,'witness_Aut':7,'prime_factor':593,
 'scope':'Finite exact controls and certificates; all-order classification follows the accompanying structural proof. No inherited spectrum omission is needed: Appendix A proves H7 absence and hypothetical H21 rigidity.'}
Path('orbit_table.json').write_text(json.dumps(orbit_rows,indent=2)+'\n')
Path('small_prime_certificates.json').write_text(json.dumps(small_prime_cert,indent=2)+'\n')
Path('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,sort_keys=True))
