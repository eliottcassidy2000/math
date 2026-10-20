#!/usr/bin/env python3
"""Independently classify every accepted normalizer representative and audit signs
by a full vertex stabilizer-chain search, not the census twin quotient.
Only masks/normalizer representatives are input from the first run.
"""
from collections import defaultdict,Counter
from itertools import combinations
import json,time
from pathlib import Path
import independent_engine as q
ROOT=Path(__file__).resolve().parent
CHECKS=0

def require(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)

def matrix(mask,n,T):
 A=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(combinations(range(n),2)):
  b=(mask>>k)&1
  if T:A[u][v],A[v][u]=1-b,b
  else:A[u][v]=A[v][u]=b
 return tuple(map(tuple,A))

def transporter(A,B,ac,bc,forced=None):
 n=len(A);f=dict(forced or {});used=set(f.values())
 if len(used)!=len(f) or any(ac[i]!=bc[j] for i,j in f.items()):return None
 if any(A[i][k]!=B[j][v] for i,j in f.items() for k,v in f.items()):return None
 candidates={i:[j for j in range(n) if ac[i]==bc[j]] for i in range(n)}
 def rec():
  if len(f)==n:return tuple(f[i] for i in range(n))
  best=None;choices=None
  for i in range(n):
   if i in f:continue
   opts=[j for j in candidates[i] if j not in used and all(A[i][u]==B[j][v] and A[u][i]==B[v][j] for u,v in f.items())]
   if not opts:return None
   if choices is None or len(opts)<len(choices):best,choices=i,opts
  for j in choices:
   f[best]=j;used.add(j);r=rec()
   if r is not None:return r
   used.remove(j);del f[best]
  return None
 return rec()

def types(masks,n,T):
 buckets=defaultdict(list);out=[]
 for mask in masks:
  A=matrix(mask,n,T);ac=q.vertex_colors(A);key=tuple(sorted(ac));found=None
  for i in buckets[key]:
   B,bc,count,rep=out[i]
   if transporter(A,B,ac,bc) is not None:found=i;break
  if found is None:
   found=len(out);out.append([A,ac,0,mask]);buckets[key].append(found)
  out[found][2]+=1
 return out

def sign_audit(A,colors):
 n=len(A);tested=0;order=1
 for i in range(n):
  orbit=1;fixed={j:j for j in range(i)}
  for j in range(i+1,n):
   if colors[i]!=colors[j]:continue
   # Directly known twin transpositions are full-vertex permutations, not quotient reductions.
   if all(A[i][v]==A[j][v] for v in range(n) if v not in (i,j)):
    p=list(range(n));p[i],p[j]=j,i
   else:p=transporter(A,A,colors,colors,fixed|{i:j})
   if p is None:continue
   require(all(A[u][v]==A[p[u]][p[v]] for u in range(n) for v in range(n)),'valid full automorphism')
   tested+=1;orbit+=1
   if sum(A[u][v] for u,v in combinations(range(n),2) if p[u]>p[v])%2:return False,tested,None
  order*=orbit
 return True,tested,order

def readmasks(row,T):
 n=row['n'];gens=row['generators'];H=q.generated_group(gens,n)
 unseen=set(combinations(range(n),2));pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)};orbits=[];base=0
 while unseen:
  u,v=min(unseen);arcs={(h[u],h[v]) for h in H};edges={tuple(sorted(e)) for e in arcs}
  require(len(arcs)==len(edges),'pair orientation consistency')
  orbits.append(sum(1<<index[e] for e in edges));unseen-=edges
  base|=sum(1<<index[tuple(sorted((u,v)))] for u,v in arcs if u>v)
 tag=f'{n}_{len(H)}_{len(orbits)}_{"T" if T else "G"}.orbits'
 out=[]
 for line in (ROOT/tag).read_text().splitlines():
  s,_,_=map(int,line.split());out.append((base if T else 0)^sum(o for i,o in enumerate(orbits) if s>>i&1))
 return out

def main():
 receipts=[]
 for row in json.loads((ROOT/'large13.json').read_text())['rows']:
  start=time.monotonic();n=row['n'];T=types(readmasks(row,True),n,True);G=types(readmasks(row,False),n,False)
  require(len(T)==row['tau_T'],'independent tournament classification')
  print(row['class'],'types',len(T),len(G),'seconds',time.monotonic()-start,flush=True)
  even=[];odd=[];tested=0;maxorder=0
  for k,(A,c,mult,rep) in enumerate(G):
   yes,tests,order=sign_audit(A,c);tested+=tests
   if order:maxorder=max(order,maxorder)
   (even if yes else odd).append(mult)
   if k%2000==0:print('sign audit',k,'seconds',time.monotonic()-start,flush=True)
  require(len(even)==row['tau_E'],'independent global sign count')
  require(dict(Counter(even))=={int(k):v for k,v in row['even_graph_multiplicity_histogram'].items()},'even histogram')
  require(dict(Counter(odd))=={int(k):v for k,v in row['escape_graph_multiplicity_histogram'].items()},'odd histogram')
  r={'class':row['class'],'tau_T':len(T),'tau_E':len(even),'difference':len(T)-len(even),'even_multiplicities':dict(Counter(even)),'odd_multiplicities':dict(Counter(odd)),'full_vertex_automorphisms_tested':tested,'max_even_automorphism_order':maxorder,'seconds':time.monotonic()-start};receipts.append(r);print(r,flush=True)
 (ROOT/'independent_large_audit.json').write_text(json.dumps({'status':'PASS','checks':CHECKS,'rows':receipts},indent=2)+'\n')
if __name__=='__main__':main()
