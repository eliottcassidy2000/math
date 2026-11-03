#!/usr/bin/env python3
"""Independent rank reconstruction from exterior contexts, with hostile controls."""
import hashlib,itertools as it,json,math,sys
from pathlib import Path
import numpy as np
from audit_carrier import groups
HERE=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def validate(N,G):
 if N<6:raise ValueError('too few vertices')
 if len(set(tuple(sorted(e)) for e in G))!=len(G):raise ValueError('repeated edge')
 for a,b in G:
  if not(0<=a<N and 0<=b<N) or a==b:raise ValueError('bad edge')
  if max(a,b)<6:raise ValueError('old-old edge in exterior context')
def coefficient(N,edges):
 adjacency=[[] for _ in range(N)]
 for a,b in edges:
  adjacency[a].append(b);adjacency[b].append(a)
 if any(len(a)>2 for a in adjacency):return None
 remaining=set(range(N));sgn=1;components=0;nontrivial=0
 while remaining:
  start=next((v for v in remaining if len(adjacency[v])<=1),None)
  if start is None:return None
  components+=1
  if not adjacency[start]:remaining.remove(start);continue
  nontrivial+=1;prev=None;at=start;length=0
  while True:
   if at not in remaining:return None
   remaining.remove(at)
   onward=[v for v in adjacency[at] if v!=prev]
   if not onward:break
   v=onward[0];sgn*=1 if at<v else -1
   prev,at=at,v;length+=1
  if length%2:return None
 return sgn,nontrivial,components

def rank_mod(rows,p):
 A=np.array(rows,dtype=np.int64)%p;r=0;det=1
 for c in range(A.shape[1]):
  piv=next((i for i in range(r,len(A)) if A[i,c]),None)
  if piv is None:continue
  if piv!=r:A[[r,piv]]=A[[piv,r]];det=-det
  z=int(A[r,c]);det=det*z%p;A[r,c:]=A[r,c:]*pow(z,-1,p)%p
  # Independent full rectangular elimination, not producer's online insertion.
  for i in range(r+1,len(A)):
   if A[i,c]:A[i,c:]=(A[i,c:]-A[i,c]*A[r,c:])%p
  r+=1
  if r==len(A):break
 return r,det%p if r==len(A)==A.shape[1] else None

def main():
 path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'rank_contexts.txt'
 contexts=[]
 for line in path.read_text().splitlines():
  a=list(map(int,line.split()));N,k=a[:2];need(len(a)==2+2*k,'context parse')
  G=list(zip(a[2::2],a[3::2]));validate(N,G);contexts.append((N,G))
 _,gg=groups(6);es=list(it.combinations(range(6),2))
 reps=[min(g) for key,g in sorted(gg.items())]
 need(len(reps)==len(contexts)==571,'square dimension')
 raw=[]
 for N,G in contexts:
  row=[]
  for m,s in reps:
   F=[e for j,e in enumerate(es) if m>>j&1]
   q=coefficient(N,F+G)
   row.append(None if q is None else (q[0]*s,q[1],q[2]))
  raw.append(row)
 results=[]
 for p in (101,103):
  rows=[[0 if q is None else q[0]*pow(2,q[1],p)*math.factorial(q[2])*pow(pow(2,q[2]-1,p),-1,p)%p for q in row] for row in raw]
  rank,det=rank_mod(rows,p);need(rank==571 and det!=0,'full modular rank')
  deletion,_=rank_mod(rows[:-1],p);need(deletion==570,'lost-row hostile')
  corrupted=rows[:-1]+[rows[0]];dup,_=rank_mod(corrupted,p);need(dup==570,'duplicate-row hostile')
  results.append({'prime':p,'rank':rank,'determinant':det,'remove_last_rank':deletion,'replace_last_by_first_rank':dup})
 for bad in [(6,[(0,1)]),(7,[(0,6),(6,0)]),(7,[(6,6)])]:
  rejected=False
  try:validate(*bad)
  except ValueError:rejected=True
  need(rejected,'invalid exterior context rejected')
 # Direct permutation coefficient controls, all edge subsets through N=5.
 for N in range(1,6):
  ee=list(it.combinations(range(N),2));idx={e:j for j,e in enumerate(ee)};numer=[0]*(1<<len(ee))
  for word in it.permutations(range(N)):
   directed=[(idx[tuple(sorted((a,b)))],1 if a<b else -1) for a,b in zip(word,word[1:])]
   for select in range(1<<(N-1)):
    m=0;s=1
    for j,(e,z) in enumerate(directed):
     if select>>j&1:m|=1<<e;s*=z
    numer[m]+=s*(1<<m.bit_count())
  for m,v in enumerate(numer):
   q=coefficient(N,[e for j,e in enumerate(ee) if m>>j&1])
   expect=0 if q is None else q[0]*2**q[1]*math.factorial(q[2])*2**(N-q[2])
   need(v==expect,'full Hamilton coefficient formula permutation control')
 out={'source_contexts_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'n':6,'coordinates':571,'queries':len(contexts),'new_vertex_range':[min(N-6 for N,G in contexts),max(N-6 for N,G in contexts)],'exterior_edge_range':[min(len(G) for N,G in contexts),max(len(G) for N,G in contexts)],'ranks':results,'formula_controls':'all Hamilton permutations and all edge monomials N1..5','checks':checks}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
