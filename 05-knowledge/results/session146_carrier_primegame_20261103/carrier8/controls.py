#!/usr/bin/env python3
"""Hostile context-validation and direct permutation/Fourier controls."""
import itertools as it, math, subprocess, tempfile
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x: raise RuntimeError(msg)

def direct_numerator(n,edges):
 edges={tuple(sorted(e)) for e in edges}; total=0
 for word in it.permutations(range(n)):
  adj={tuple(sorted((a,b))): (1 if a<b else -1) for a,b in zip(word,word[1:])}
  if edges<=adj.keys():
   sign=math.prod(adj[e] for e in edges)
   total+=sign*(2**len(edges))
 return total

def formula_numerator(n,edges):
 adj=[[] for _ in range(n)]
 for a,b in edges:adj[a].append(b);adj[b].append(a)
 if any(len(x)>2 for x in adj):return 0
 seen=set();components=0;nontrivial=0;sgn=1
 for start in range(n):
  if start in seen or len(adj[start])>1:continue
  components+=1;seen.add(start)
  if not adj[start]:continue
  nontrivial+=1;prev=None;at=start;length=0
  while True:
   onward=[x for x in adj[at] if x!=prev]
   if not onward:break
   nxt=onward[0]
   if nxt in seen:return 0
   sgn*=1 if at<nxt else -1
   seen.add(nxt);length+=1;prev,at=at,nxt
  if length%2:return 0
 if len(seen)<n:return 0
 return sgn*2**nontrivial*math.factorial(components)*2**(n-components)

cases=[[],[(0,1)],[(0,3),(3,1)],[(0,7),(7,1),(2,6),(6,3)],[(0,2),(2,3),(3,4),(4,7)],[(0,1),(1,2),(2,0)],[(0,1),(0,2),(0,3)],[(0,3),(3,1),(1,6),(6,2),(2,5),(5,4)],[(0,2),(2,4),(4,6),(6,7),(7,5),(5,3),(3,1)]]
for es in cases:
 need(direct_numerator(8,es)==formula_numerator(8,es),f'permutation mismatch {es}')
# Exact Walsh extraction for a nontrivial exterior monomial from actual fixed cuts.
# Old vertices 0,1; exterior vertices 2,3; old edge orientation is fixed at +1/2.
n=4;old={(0,1)};ext=[e for e in it.combinations(range(n),2) if e not in old];G={(0,2),(1,2)}
total=Fraction(0)
for signs in it.product((-1,1),repeat=len(ext)):
 y={(0,1):Fraction(1,2)};y.update({e:Fraction(s,2) for e,s in zip(ext,signs)})
 H=Fraction(0)
 for word in it.permutations(range(n)):
  H+=math.prod(Fraction(1,2)+(y[tuple(sorted((a,b)))] if a<b else -y[tuple(sorted((a,b)))]) for a,b in zip(word,word[1:]))
 total+=H*math.prod(signs[ext.index(e)] for e in G)
walsh=total*Fraction(2**len(G),2**len(ext))
# Sum the independently checked full coefficients over old-edge choices.
predicted=Fraction(0)
for F in (set(),old):
 W=F|G; predicted+=Fraction(formula_numerator(n,W),2**(n-1))*Fraction(1,2)**len(F)
need(walsh==predicted,'Walsh extraction to fixed cut evaluations')
invalid={'old-old':'7 1 0 1\n','loop':'7 1 6 6\n','duplicate-reversed':'7 2 0 6 6 0\n','out-of-range':'7 1 0 7\n','truncated':'7 1 0\n','trailing':'7 0 99\n'}
for name,data in invalid.items():
 with tempfile.NamedTemporaryFile(mode='w') as f:
  f.write(data);f.flush()
  result=subprocess.run([str(HERE/'verify'),'6',f.name,str(HERE/'reps6.txt')],capture_output=True,text=True)
  need(result.returncode!=0 and 'FAIL:' in result.stderr, f'{name} accepted or crashed: {result.stderr}')
print('direct Hamilton permutation controls N=8:',len(cases),'passed')
print('exact fixed-cut Walsh extraction N=4:',walsh,'passed')
print('invalid contexts rejected:',', '.join(invalid))
