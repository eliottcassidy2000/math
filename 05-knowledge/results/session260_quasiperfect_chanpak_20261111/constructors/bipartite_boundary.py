#!/usr/bin/env python3
"""Rooted-tree boundary compiler, bipartite moves, and exact small spectrum."""
from pathlib import Path
from itertools import product
import json
HERE=Path(__file__).resolve().parent

def boundary(G,v):
 n=len(G);full=(1<<n)-1;dom=[0]*(1<<n);out=[0,0,0]
 for S in range(1<<n):
  if S:
   b=S&-S;dom[S]=dom[S^b]|b|G[b.bit_length()-1]
  if dom[S]&(full^(1<<v))!=full^(1<<v):continue
  out[0 if S>>v&1 else 1 if dom[S]>>v&1 else 2]+=1
 return tuple(out)
def count(G):return sum(boundary(G,0)[:2])
def leaf(G,v):
 H=G+[1<<v];H[v]|=1<<len(G);return H
def twins(G,v,m):
 old=[i for i in range(len(G)) if i!=v];k=len(old);H=[0]*(k+m)
 for i,x in enumerate(old):
  for j,y in enumerate(old):
   if G[x]>>y&1:H[i]|=1<<j
  if G[x]>>v&1:
   H[i]|=((1<<m)-1)<<k
   for a in range(m):H[k+a]|=1<<i
 return H

def shore_state(G,L):
 n=len(G);R=[i for i in range(n) if i not in L]
 def F(shore):
  ans=0
  for mask in range(1<<len(shore)):
   hit=0
   for i,v in enumerate(shore):
    if mask>>i&1:hit|=G[v]
   ans+=1<<hit.bit_count()
  return ans
 return dict(d=count(G),FL=F(L),FR=F(R),P=1<<n,eL=int(all(G[v] for v in L)),eR=int(all(G[v] for v in R)))
def hub_left(G,L):
 n=len(G);R=[i for i in range(n) if i not in L];H=G+[sum(1<<v for v in R)]
 for v in R:H[v]|=1<<n
 return H,L+[n]

def main():
 states={1:{(1,0,1)}};mins=[]
 for n in range(1,11):
  if n>1:
   states[n]={(A*(a+b+c),B*(a+b)+C*a,C*b) for i in range(1,n) for A,B,C in states[i] for a,b,c in states[n-i]}
  mins.append(dict(n=n,rooted_boundary_states=len(states[n]),minimum=min(A+B for A,B,C in states[n])))
 if mins[6]['minimum']!=45 or mins[9]['minimum']!=225:raise RuntimeError('tree cutoff')
 census=json.loads((HERE/'bipartite_spectrum_receipt.json').read_text());values=set().union(*(set(r['values']) for r in census))
 closed={1}
 for x in range(3,225,2):
  if x in values or any(x%d==0 and d in closed and x//d in closed for d in range(3,x,2)):closed.add(x)
 holes=[x for x in range(1,225,2) if x not in closed]
 if holes!=[7,13,21,43]:raise RuntimeError('spectrum')
 if not closed<=values:raise RuntimeError('direct witnesses needed')
 # Literal non-chain and chain interfaces, all roots, not a separate census.
 samples=[];checks=0
 for row in census:
  if row['n']>7:continue
  for d,p,code in row['witnesses']:
   if d not in [3,5,9,11,15,17,19,23,29,31,37,41,45,47]:continue
   n=row['n'];q=n-p;G=[0]*n
   for a,b in product(range(p),range(q)):
    if code>>(a*q+b)&1:G[a]|=1<<(p+b);G[p+b]|=1<<a
   L=list(range(p));z=shore_state(G,L);H,LL=hub_left(G,L);zz=shore_state(H,LL)
   expect=dict(d=z['d']+z['FR']-z['eR'],FL=z['FL']+z['P'],FR=2*z['FR']-1,P=2*z['P'],eL=z['eL'],eR=1)
   if zz!=expect:raise RuntimeError('shore hub')
   checks+=1
   if n<=6 and d in [5,11,29]:
    for aa,bb in product(range(1,3),repeat=2):
     HH=G+[0]*(aa+bb);LL=set(L)|set(range(n,n+aa));nn=n+aa+bb
     for i in range(nn):
      for j in range(max(i+1,n),nn):
       if (i in LL)!=(j in LL):HH[i]|=1<<j;HH[j]|=1<<i
     pred=z['P']*(2**aa-1)*(2**bb-1)+(z['FR']-1)*(2**aa-1)+(z['FL']-1)*(2**bb-1)+d+2-z['eL']-z['eR']
     if count(HH)!=pred:raise RuntimeError('two-shore formula')
     checks+=1
   for v in range(n):
    A,B,C=boundary(G,v)
    if boundary(leaf(G,v),v)!=(2*A,B+C,0):raise RuntimeError('leaf')
    if G[v] and C>B:raise RuntimeError('C injection')
    # alpha counts outside sets meeting N(v), valid when v selected.
    alpha=0
    for S in range(1<<n):
     if S>>v&1 or not S&G[v]:continue
     hit=G[v]|(1<<v)
     for u in range(n):
      if S>>u&1:hit|=G[u]|(1<<u)
     if hit==(1<<n)-1:alpha+=1
    beta=d-2*alpha
    for m in (1,2,3):
     if count(twins(G,v,m))!=(1<<m)*alpha+beta:raise RuntimeError('false twins')
    checks+=4
   if d in [5,11,29]:samples.append(dict(d=d,rows=G,shore_state=z))
 out=dict(status='PASS',tree_minima=mins,global_bipartite_exceptions_below225=holes,scope='Complete below225 using connected order10 lower bound225 and component products; not a global four-exception theorem',boundary_controls=checks,sample_interfaces=samples)
 (HERE/'bipartite_boundary_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
