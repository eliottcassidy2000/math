#!/usr/bin/env python3
"""Standalone literal controls for C1_COMPLETE.md; explicit -O-safe checks."""
from itertools import combinations,permutations
from pathlib import Path
import json,math
P=Path(__file__).resolve().parent

def require(ok,why):
 if not ok:raise RuntimeError(why)

def graphs(n):
 pairs=list(combinations(range(n),2))
 for code in range(1<<len(pairs)):
  rows=[0]*n
  for i,(u,v) in enumerate(pairs):
   if code>>i&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
  yield rows

def cyclic(rows,word):
 return len(word)==len(rows) and len(set(word))==len(rows) and all(rows[u]>>v&1 for u,v in zip(word,word[1:]+word[:1]))
def cycles(rows):
 if len(rows)<3:return []
 return [p for rest in permutations(range(1,len(rows))) for p in [(0,)+rest] if cyclic(rows,p)]
def induced(rows,vs):return [sum(1<<j for j,v in enumerate(vs) if rows[u]>>v&1) for u in vs]
def components(rows):
 n=len(rows);reach=[r|1<<u for u,r in enumerate(rows)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 left=set(range(n));parts=[]
 while left:
  u=min(left);part=tuple(v for v in left if reach[u]>>v&1 and reach[v]>>u&1);left.difference_update(part);parts.append(part)
 return sorted(parts,key=lambda D:-sum(bool(rows[D[0]]>>E[0]&1) for E in parts if E!=D))
def strong(rows):return len(components(rows))==1

def path_counts(rows,vs):
 h=0;s={u:0 for u in vs};e=s.copy()
 for p in permutations(vs):
  if all(rows[u]>>v&1 for u,v in zip(p,p[1:])):h+=1;s[p[0]]+=1;e[p[-1]]+=1
 return h,s,e

def main():
 labelled=unique=two_strong=deletion_checks=products=bridges=0
 for n in range(1,6):
  for rows in graphs(n):
   labelled+=1;cs=cycles(rows)
   if not cs:continue
   deletions=[]
   for root in range(n):
    keep=[u for u in range(n) if u!=root];r=induced(rows,keep);parts=components(r)
    if len(parts)==1:
     if len(cs)==1:
      smaller=tuple(keep.index(u) for u in cs[0] if u!=root)
      require(cyclic(r,smaller),'deletion lemma')
      require(len(cycles(r))==1,'smaller uniqueness');deletion_checks+=1
     continue
    Ds=[tuple(keep[u] for u in D) for D in parts];stats=[path_counts(rows,D) for D in Ds]
    A=sum(stats[0][1][a] for a in Ds[0] if rows[root]>>a&1)
    B=sum(stats[-1][2][b] for b in Ds[-1] if rows[b]>>root&1)
    expected=A*B*math.prod(v[0] for v in stats[1:-1])
    require(expected==len(cs),'separating-vertex product');products+=1;deletions.append(root)
    if len(cs)==1:
     aa=[a for a in Ds[0] if rows[root]>>a&1];require(len(aa)==1,'unique initial outneighbor')
     copy=rows[:];copy[root]^=1<<aa[0]
     require(not strong(copy),'exposed arc not bridge');bridges+=1
   if not deletions:
    require(len(cs)>=2,'2-strong lower bound');two_strong+=1
   if len(cs)==1:
    require(deletions,'unique cycle without separating vertex');unique+=1
 surgery=0
 for n in range(5,501):
  forced=[(1<<((i+1)%n))|(1<<((i+2)%n)) for i in range(n)]
  original=tuple(range(n))
  if n%2:
   word=tuple(range(0,n,2))+tuple(range(1,n,2))
   require(cyclic(forced,word) and word!=original,'odd step2');surgery+=1
  else:
   a=forced[:];a[4]|=1<<1
   word=(0,2,4,1,3)+tuple(range(5,n))
   require(cyclic(a,word) and word!=original,'local even surgery');surgery+=1
   b=forced[:];b[1]|=1<<4
   word=tuple(range(3,n,2))+(1,)+tuple(range(4,n,2))+(0,2)
   require(cyclic(b,word),'even cycle merge');j=word.index(0);require(word[j:]+word[:j]!=original,'distinct merge');surgery+=1
 result=dict(status='PASS',labelled_tournaments=labelled,maximum_labelled_order=5,uniquely_Hamiltonian_labelled=unique,two_strong_labelled=two_strong,deletion_lemma_checks=deletion_checks,separating_vertex_cycle_products=products,exposed_strong_bridges=bridges,forced_arc_surgeries=surgery,surgery_orders=[5,500],all_order_claim='Every strong C=1 tournament has a strong bridge.')
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'c1_complete_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
