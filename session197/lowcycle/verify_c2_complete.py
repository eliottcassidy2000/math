#!/usr/bin/env python3
"""Complete small base and literal insertion controls for C2_COMPLETE.md."""
import itertools as it,json,math
from pathlib import Path
P=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)
def graphs(n):
 pairs=list(it.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  r=[0]*n
  for k,(u,v) in enumerate(pairs):
   if code>>k&1:r[u]|=1<<v
   else:r[v]|=1<<u
  yield tuple(r)
def induced(r,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if r[u]>>v&1) for u in vs)
def strong(r):
 n=len(r);reach=[r[u]|1<<u for u in range(n)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 return all(x==(1<<n)-1 for x in reach)
def cycles(r):
 n=len(r)
 if n<3:return []
 return [(0,)+p for p in it.permutations(range(1,n)) if all(r[u]>>v&1 for u,v in zip((0,)+p,p+(0,)))]
def cyckey(w):
 j=w.index(min(w));return w[j:]+w[:j]
def reg_iso(r):
 R=(6,12,24,17,3)
 return any(induced(r,p)==R for p in it.permutations(range(5)))

def main():
 total=insertions=degreecovers=0;stats=[];labelledC2=[]
 for n in range(1,7):
  hist={};allc2=0;two=0
  for r in graphs(n):
   total+=1
   if n<3 or not strong(r):continue
   cs=cycles(r);C=len(cs)
   if C==2:allc2+=1
   deletions=[]
   for v in range(n):
    vs=[u for u in range(n) if u!=v];small=induced(r,vs)
    if strong(small):deletions.append((v,vs,small))
   if len(deletions)!=n:continue
   two+=1;hist[C]=hist.get(C,0)+1
   if n==5:need(C==2 and reg_iso(r),'regular5 classification')
   if n==6:need(C in(4,5),'order6 base')
   all_degree_one=True
   for v,vs,small in deletions:
    Ds=cycles(small);images=set()
    for D in Ds:
     old=tuple(vs[u] for u in D);j=next((j for j in range(n-1) if r[old[j]]>>v&1 and r[v]>>old[(j+1)%(n-1)]&1),None)
     need(j is not None,'missing cycle insertion');new=old[:j+1]+(v,)+old[j+1:]
     key=cyckey(new);need(key in cs,'insertion not a cycle');images.add(key)
    need(len(images)==len(Ds)<=C,'insertion injectivity');insertions+=1
    low=[u for u,x in enumerate(small) if x.bit_count()==1]
    if not low:all_degree_one=False
    for u in low:
     old=vs[u];need(r[old].bit_count()==2 and r[old]>>v&1,'degree-two covering edge')
   if all_degree_one:
    L=[u for u,x in enumerate(r) if x.bit_count()==2];ell=len(L)
    need(all(any(r[u]>>v&1 for u in L) for v in range(n)),'degree-two set does not cover')
    out=sum(bool(r[u]>>v&1) for u in L for v in range(n) if v not in L)
    need(out==2*ell-ell*(ell-1)//2 and n-ell<=out and 1<=ell<=5 and n<=6,'degree cover bound');degreecovers+=1
  stats.append(dict(order=n,two_strong_labelled=two,cycle_histogram=hist));labelledC2.append(allc2)
 need(stats[4]['cycle_histogram']=={2:24},'n5 histogram')
 need(stats[5]['cycle_histogram']=={4:720,5:1920},'n6 histogram')
 R=(6,12,24,17,3);aut=sum(induced(R,p)==R for p in it.permutations(range(5)))
 need(aut==5,'regular5 automorphisms')
 starts=[]
 for u in range(5):starts.append(sum(all(R[a]>>b&1 for a,b in zip(w,w[1:])) for rest in it.permutations([v for v in range(5) if v!=u]) for w in [(u,)+rest]))
 need(starts==[3]*5,'regular5 start counts')
 need(labelledC2==[0,0,0,0,144,3600],'labelled C2 coefficients')
 bounds=[ell*(7-ell)//2 for ell in range(1,6)];need(bounds==[3,5,6,6,5],'integer cover maxima')
 # Explicit elementary order-six controls, independent of the census filter.
 regular_insertions=0
 for k in (2,3):
  for outs in it.combinations(range(5),k):
   r=list(R)+[sum(1<<v for v in outs)]
   for u in range(5):
    if u not in outs:r[u]|=1<<5
   images=set()
   for old in cycles(R):
    for j in range(5):
     if r[old[j]]>>5&1 and r[5]>>old[(j+1)%5]&1:
      images.add(cyckey(old[:j+1]+(5,)+old[j+1:]))
   need(len(images)==3 and all(w in cycles(r) for w in images),'regular-five insertion base')
   regular_insertions+=1
 triangle_cases=0
 for reverse in (False,True):
  r=[0]*6
  for u,v in ((0,1),(1,2),(2,0)):r[u]|=1<<v
  for u,v in ((3,4),(4,5),(5,3)):
   if reverse:u,v=v,u
   r[u]|=1<<v
  for u in range(3):
   for v in range(3,6):
    if v==u+3:r[u]|=1<<v
    else:r[v]|=1<<u
  words=[(0,3,1,4,2,5),(0,3,2,5,1,4),((0,1,4,3,2,5) if reverse else (0,1,2,5,3,4))]
  need(len(set(words))==3 and all(w in cycles(r) for w in words),'matched-triangle base cycles')
  triangle_cases+=1
 result=dict(status='PASS',regular_five_insertion_cases=regular_insertions,matched_triangle_cases=triangle_cases,labelled_tournaments=total,complete_max_order=6,two_strong_inventory=stats,literal_cycle_insertion_checks=insertions,physical_degree_cover_cases=degreecovers,degree_cover_bounds=bounds,labelled_C2_counts_by_order=labelledC2,R5_automorphism_order=aut,R5_root_start_counts=starts,all_order_classification='Every strong C2 tournament is regular5 or has a strong bridge.')
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'c2_complete_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
