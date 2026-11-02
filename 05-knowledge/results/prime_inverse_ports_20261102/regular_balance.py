#!/usr/bin/env python3
"""Literal cycles versus a trace/degree proof, including hostile controls."""
from pathlib import Path
from itertools import combinations, permutations, product
from fractions import Fraction
from math import comb
import json, random
from research import circ, decode, need, ordered_covers, all_h, skeleton, substitute
HERE=Path(__file__).resolve().parent

def literal(rows):
 n=len(rows);tri=[];five=0
 for sub in combinations(range(n),3):
  if all(sum(rows[u]>>v&1 for v in sub)==1 for u in sub):tri.append(sub)
 for sub in combinations(range(n),5):
  for tail in permutations(sub[1:]):
   w=(sub[0],)+tail
   five+=all(rows[u]>>v&1 for u,v in zip(w,w[1:]+w[:1]))
 masks=[sum(1<<u for u in t) for t in tri]
 pairs=sum(not(a&b) for i,a in enumerate(masks) for b in masks[i+1:])
 edge={(u,v):0 for u in range(n) for v in range(n) if rows[u]>>v&1}
 vertices=[0]*n
 for t in tri:
  for u in t:
   vertices[u]+=1
   for v in t:
    if (u,v) in edge:edge[u,v]+=1
 return dict(C3=len(tri),C5=five,M33=pairs,vertex_triangles=vertices,edge_triangles=list(edge.values()))

def check(rows):
 n=len(rows);need(all(r.bit_count()==(n-1)//2 for r in rows),'regular input')
 a=literal(rows);c=n*(n*n-1)//24;t=(n*n-1)//8
 need(a['C3']==c and a['vertex_triangles']==[t]*n,'triangle count')
 square=sum(x*x for x in a['edge_triangles'])
 K=Fraction(n*(n*n-1)*(n*n+11),160)
 balance=Fraction(n*(n*n-1)*(n-3)*(5*n*n-12*n-41),2880)
 need(a['C5']+square==K,'trace law checked by literal cycles')
 need(a['C5']+2*a['M33']==balance,'disjoint cycle balance')
 need(a['M33']==comb(c,2)-sum(comb(x,2) for x in a['vertex_triangles'])+sum(comb(x,2) for x in a['edge_triangles']),'pair inclusion-exclusion')
 upper=Fraction(n*(n*n-1)*(n-2)*(n-3),160)
 if n%4==1:upper-=Fraction(n*(n-1),8)
 need(a['C5']<=upper,'known C5 bound and integer correction')
 # Independently multiply the skew matrix; do not use edge counts here.
 S=[[0 if u==v else 1 if rows[u]>>v&1 else -1 for v in range(n)] for u in range(n)]
 S2=[[sum(S[u][w]*S[w][v] for w in range(n)) for v in range(n)] for u in range(n)]
 tr4=sum(x*x for row in S2 for x in row)
 need(tr4==32*square-n**4+n,'skew trace bridge')
 return dict(n=n,H=None,C3=c,C5=a['C5'],M33=a['M33'],balance=int(balance),edge_triangle_square_sum=square)

def main():
 results=[];counts={}
 for n in [1,3,5,7]:
  codes=[''] if n==1 else (HERE/'inputs'/f'quotients_n{n}.txt').read_text().split()
  for code in codes:
   rows=decode(code,n)
   if all(r.bit_count()==(n-1)//2 for r in rows):results.append(check(rows));counts[str(n)]=counts.get(str(n),0)+1
 for n in [9,11]:
  for bits in product([0,1],repeat=(n-1)//2):results.append(check(circ(n,bits)));counts[str(n)]=counts.get(str(n),0)+1
 # A degree-preserving random walk is a finite hostile sample, not a class census.
 rng=random.Random(20261102)
 for n in [13,15]:
  rows=circ(n,[0]*((n-1)//2))
  for sample in range(16):
   for step in range(20):
    ts=[t for t in combinations(range(n),3) if all(sum(rows[u]>>v&1 for v in t)==1 for u in t)]
    t=rng.choice(ts)
    for u,v in combinations(t,2):rows[u]^=1<<v;rows[v]^=1<<u
   results.append(check(rows));counts[str(n)]=counts.get(str(n),0)+1
 # Same score list is insufficient outside regularity. Search all small classes.
 hostile=None
 for n in range(3,7):
  seen={}
  for code in (HERE/'inputs'/f'quotients_n{n}.txt').read_text().split():
   rows=decode(code,n);a=literal(rows);key=tuple(sorted(r.bit_count() for r in rows));value=a['C5']+2*a['M33']
   if key in seen and seen[key]['balance']!=value:
    hostile=dict(n=n,scores=key,first=seen[key],second=dict(rows=rows,balance=value,H=all_h(rows)[-1]));break
   seen[key]=dict(rows=rows,balance=value,H=all_h(rows)[-1])
  if hostile:break
 need(hostile is not None,'same-score hostile exists')
 higher=[];modules=[]
 for bits in [(1,1,1,1),(1,1,1,0),(1,0,1,1)]:
  rows=circ(9,bits);C=ordered_covers(rows);alpha=Fraction(C[2]-31110,6);beta=C[1]-1-60-768-alpha
  for r in range(12):
   need(sum(comb(r,j)*C[j] for j in range(1,min(9,r)+1))==r**9+60*r**7+768*r**5+alpha*r**3+beta*r,'order9 polynomial')
  higher.append(dict(bits=bits,H=C[1],C2=C[2],alpha=int(alpha),beta=int(beta)))
  modules.append(rows)
 contexts=0
 for n in range(1,5):
  codes=[''] if n==1 else (HERE/'inputs'/f'quotients_n{n}.txt').read_text().split()
  for code in codes:
   q=decode(code,n)
   for v in range(n):
    b=skeleton(q,v)+[0,0,0]
    actual=[all_h(substitute(q,v,B))[-1] for B in modules]
    for i in [1,2]:
     dh=higher[i]['H']-higher[0]['H'];da=higher[i]['alpha']-higher[0]['alpha']
     need(actual[i]-actual[0]==dh*b[1]+6*da*(b[2]+b[3]),'nine module repair in full adjacency')
     contexts+=1
 out=dict(status='PASS',literal_regular_checks=len(results),universes=counts,
          scope='All regular classes through7; all cyclic orientations at9,11; seeded triangle-reversal samples at13,15',
          results=results,same_score_hostile=hostile,order9_polynomials=higher,order9_full_adjacency_exchange_checks=contexts)
 (HERE/'regular_balance.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print('PASS',len(results),'regular controls; hostile order',hostile['n'])
if __name__=='__main__':main()
