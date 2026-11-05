#!/usr/bin/env python3
"""Independent graph/certificate audit; no source-checker imports."""
import hashlib
import itertools
import json
import math
from pathlib import Path

SOURCE=Path('/workspace/scratch/07d5532f2fce/session197-orbit-wildcard')
EXPECTED={'PROOF.md':'69497d610f4ffdb3530ed260cec3618af3d973aadaadbbddd23aaff343424f4a',
          'verify.py':'9f354c23dfe8616f73f58cc15f56b6054a5f197bdb51fc944c6550b04d6ceba9',
          'receipt.json':'50f7e00c1375e818bdc01f419d0c13b2340d59caee197cd2e2aaeae3e5c31aef',
          'witness.json':'b4358dec35824a42a5756bd73a2cd180a3b7fd242870628fb1e54ae5430d36d7'}

def check(ok,message):
 if not ok:raise RuntimeError(message)

def graph():
 n=14;out=[0]*n
 for x in range(n):
  for y in range(x+1,n):
   bx,i=divmod(x,7);by,j=divmod(y,7)
   win=((j-i)%7 in {1,2,4}) if bx==by else ((j-i)%7 in {1,3})
   a,b=(x,y) if win else (y,x)
   out[a]|=1<<b
 return out

def h_walk_inclusion_exclusion(rows):
 """Length n-1 walks using all n vertices are exactly Hamilton paths."""
 n=len(rows);answer=0;terms=0
 for mask in range(1,1<<n):
  vertices=[i for i in range(n) if mask>>i&1]
  arcs=[(i,j) for i in vertices for j in vertices if rows[i]>>j&1]
  counts=[int(mask>>i&1) for i in range(n)]
  for _ in range(n-1):
   nxt=[0]*n
   for i,j in arcs:nxt[j]+=counts[i]
   counts=nxt
  answer+=(-1 if (n-len(vertices))&1 else 1)*sum(counts)
  terms+=1
 return answer,terms

def literal_h(rows):
 return sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(rows)-1))
            for p in itertools.permutations(range(len(rows))))

def modules(rows):
 n=len(rows);full=(1<<n)-1;answer=[]
 for M in range(1,full+1):
  if M.bit_count()<2 or M==full:continue
  exterior=full^M
  first=(M&-M).bit_length()-1
  pattern=rows[first]&exterior
  if all(rows[u]&exterior==pattern for u in range(n) if M>>u&1):answer.append(M)
 return answer

def automorphisms(rows):
 n=len(rows);degrees=[r.bit_count() for r in rows]
 images=[-1]*n;used=0;output=[];nodes=0
 def visit(used):
  nonlocal nodes
  nodes+=1
  remaining=[v for v in range(n) if images[v]<0]
  if not remaining:
   output.append(images[:]);return
  option=[]
  for v in remaining:
   possible=[]
   for w in range(n):
    if used>>w&1 or degrees[w]!=degrees[v]:continue
    if all(((rows[v]>>u)&1)==((rows[w]>>images[u])&1)
           for u in range(n) if images[u]>=0):possible.append(w)
   if not possible:return
   option.append((len(possible),v,possible))
  _,v,possible=min(option)
  for w in possible:
   images[v]=w;visit(used|(1<<w));images[v]=-1
 visit(0)
 return output,nodes

def trial_prime(n):
 return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def trial_factor(n):
 out={};d=2
 while d*d<=n:
  while n%d==0:out[d]=out.get(d,0)+1;n//=d
  d+=1
 if n>1:out[n]=out.get(n,0)+1
 return out

def verify_certificate(c,seen):
 p=c['p'];check(trial_prime(p),('independent certificate-node primality',p))
 if p==2:
  seen[p]={'trial_bound':1};return
 factors={}
 for z in c['factors']:
  q=z['certificate']['p']
  check(q not in factors and q<p and isinstance(z['e'],int) and z['e']>0,'factor node')
  factors[q]=z['e'];verify_certificate(z['certificate'],seen)
 check(factors==trial_factor(p-1),('independent complete factorization',p))
 a=c['a'];check(pow(a,p-1,p)==1,'Lucas Fermat condition')
 check(all(math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in factors),'Lucas full-order conditions')
 seen[p]={'trial_bound':math.isqrt(p),'factorization_minus_one':factors,'witness':a}

def main():
 for name,expected in EXPECTED.items():
  check(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==expected,('frozen hash',name))
 witness=json.loads((SOURCE/'witness.json').read_text())
 rows=graph();check(rows==witness['rows'],'independent modular construction')
 n=14;full=(1<<n)-1
 check(all(not(rows[i]>>i&1) and rows[i]<1<<n for i in range(n)),'row range')
 check(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)),'tournament')
 # Validate the independent walk formula against literal permutations first.
 small_checks=0
 for size in range(1,6):
  pairs=list(itertools.combinations(range(size),2))
  for bits in range(1<<len(pairs)):
   small=[0]*size
   for k,(i,j) in enumerate(pairs):
    a,b=(i,j) if bits>>k&1 else (j,i);small[a]|=1<<b
   check(h_walk_inclusion_exclusion(small)[0]==literal_h(small),'walk/literal control')
   small_checks+=1
 H,terms=h_walk_inclusion_exclusion(rows)
 check(H==15144269==witness['H'],'independent walk inclusion-exclusion count')
 p=H//7;check(H%7==0 and p==2163467==witness['prime_quotient'],'prime factor')
 check(trial_prime(p),'trial prime')
 certificate_nodes={};verify_certificate(witness['primality_certificate'],certificate_nodes)
 nontrivial=modules(rows);check(not nontrivial,'all-subset module exclusion')
 reach=[rows[i]|1<<i for i in range(n)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 check(all(z==full for z in reach),'strongness by transitive closure')
 auto,nodes=automorphisms(rows)
 rotations=[[7*(i//7)+(i%7+k)%7 for i in range(n)] for k in range(7)]
 check(sorted(auto)==sorted(rotations)==sorted(witness['all_automorphisms']),'exact full automorphism group')
 orbits=sorted({tuple(sorted({perm[i] for perm in auto})) for i in range(n)})
 check(orbits==[tuple(range(7)),tuple(range(7,14))],'two exact orbits')
 degrees=[r.bit_count() for r in rows];check(degrees==[5]*7+[8]*7,'degree orbit split')
 # Each block is exactly Paley, has H189, and is mixed to every opposite vertex.
 block_H=[]
 for block in [list(range(7)),list(range(7,14))]:
  sub=[sum(1<<j for j,v in enumerate(block) if rows[u]>>v&1) for u in block]
  h,_=h_walk_inclusion_exclusion(sub);block_H.append(h);check(h==189,'Paley block H')
  M=sum(1<<u for u in block)
  for v in range(n):
   if v not in block:check(rows[v]&M not in (0,M),'opposite vertex mixed to full orbit')
 check(math.factorial(14)%7==0 and math.factorial(14)//7==12454041600,'labelled mass')
 result={'status':'PASS','independent_algorithm':'vertex-cover inclusion-exclusion over length13 directed walks',
         'small_literal_walk_controls':small_checks,'walk_inclusion_exclusion_terms':terms,
         'H':H,'prime_quotient':p,'prime_trial_bound':math.isqrt(p),'verified_Lucas_nodes':certificate_nodes,
         'nontrivial_modules':nontrivial,'module_universe':1<<14,
         'automorphism_backtracking_nodes':nodes,'automorphism_order':len(auto),'automorphisms':auto,
         'orbits':orbits,'outdegrees':degrees,'block_H_values':block_H,
         'path_orbit_mass':p,'labelled_mass':math.factorial(14)//7,
         'frozen_source_hashes':EXPECTED,
         'scope':'Independent audit of the explicit two-orbit graph and certificate, not minimality or global spectrum.'}
 print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
