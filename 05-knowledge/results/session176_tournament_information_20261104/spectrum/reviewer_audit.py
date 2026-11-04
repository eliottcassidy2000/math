#!/usr/bin/env python3
"""Independent reviewer: literal path/cut enumeration, plus tested target circuit."""
import itertools as it,json,random
from check import streaming

def need(v,m):
 if not v:raise RuntimeError(m)
def literal(a):
 return sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in it.permutations(range(len(a))))
def inflate(a,r):
 origins=[0]*r+list(range(1,5));g=[0]*len(origins)
 for i,x in enumerate(origins):
  for j,y in enumerate(origins):
   if i!=j and ((x==y and i<j) or (x!=y and a[x]>>y&1)):g[i]|=1<<j
 return g
def patterns(g):
 # Enumerate words with indistinguishable module symbols; independent of gap mask implementation.
 table={}
 for p in it.permutations(range(1,5)):
  counts=[0]*5
  for k in range(1,6):
   for places in it.combinations(range(4+k),k):
    seq=[];idx=0
    for j in range(4+k):
     if j in places:seq.append(0)
     else:seq.append(p[idx]);idx+=1
    if all(x!=y and g[x]>>y&1 for x,y in zip(seq,seq[1:])):counts[k-1]+=1
  if any(counts):table[''.join(map(str,p))]=counts
 return table

def main():
 L=[16,17,3,7,12];U=[16,25,3,5,12];tables=[patterns(L),patterns(U)]
 for table,first in zip(tables,[11,13]):
  need([sum(row[j] for row in table.values()) for j in range(5)]==[first,6,0,0,0],'symbolic profiles')
 gates=[]
 for r in range(1,5):
  vals=[literal(inflate(seed,r)) for seed in [L,U]];need(vals==[6*2**r-1,6*2**r+1],'literal inflation');gates.append(vals)
 # Exhaustive tournaments on <=5 vertices, selecting every transitive initial backbone.
 tests=0
 for N in range(1,6):
  edges=list(it.combinations(range(N),2))
  for bits in range(1<<len(edges)):
   g=[0]*N
   for z,(i,j) in enumerate(edges):
    if bits>>z&1:g[i]|=1<<j
    else:g[j]|=1<<i
   h=literal(g)
   for n in range(N+1):
    if all(g[i]>>j&1 for i in range(n) for j in range(i+1,n)):
     need(streaming(g,n)==h,('stream literal',N,bits,n));tests+=1
 need(streaming([],0)==1,'empty convention')
 print(json.dumps({'status':'PASS','literal_inflation_pairs':gates,'complete_skeleton_tables':tables,'literal_streaming_comparisons':tests,'review_note':'The n=0 complexity bound was corrected to n+1 during review. No substantive identity defect found.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
