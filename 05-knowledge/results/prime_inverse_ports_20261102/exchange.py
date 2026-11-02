#!/usr/bin/env python3
from pathlib import Path
from itertools import permutations,combinations
from math import comb,gcd
from functools import reduce
import json,subprocess
from research import all_h,ordered_covers,circ,substitute,decode,need,GEN
HERE=Path(__file__).resolve().parent

def cycle_data(r):
 n=len(r);cs={k:[] for k in [3,5,7]}
 for k in cs:
  for sub in combinations(range(n),k):
   for t in permutations(sub[1:]):
    w=(sub[0],)+t
    if all(r[a]>>b&1 for a,b in zip(w,w[1:]+w[:1])):cs[k].append(frozenset(sub))
 m=sum(not(a&b) for i,a in enumerate(cs[3]) for b in cs[3][i+1:])
 return {**{str(k):len(v) for k,v in cs.items()},'33':m}

def codes(n):
 f=HERE/'inputs'/f'quotients_n{n}.txt'
 return f.read_text().split() if f.exists() else subprocess.check_output([GEN,'-q',str(n)],text=True).split()

def main():
 regular=[]
 for code in codes(7):
  r=decode(code,7)
  if not all(x.bit_count()==3 for x in r):continue
  C=ordered_covers(r);cyc=cycle_data(r)
  need(cyc['3']==14 and cyc['5']+2*cyc['33']==56,'regular seven balance')
  for x in range(13):
   need(sum(comb(x,k)*C[k] for k in range(1,min(7,x)+1))==x**7+28*x**5+112*x**3+2*cyc['7']*x,'cycle colouring polynomial')
  regular.append(dict(rows=r,covers=C,cycles=cyc))
 need(len(regular)==3,'complete regular seven class bank')
 A=circ(7,(0,0,0));B=circ(7,(0,0,1));contexts=0
 for n in range(1,6):
  bank=[''] if n==1 else codes(n)
  for code in bank:
   Q=decode(code,n);Hq=all_h(Q)[-1]
   for v in range(n):
    ha=all_h(substitute(Q,v,A))[-1];hb=all_h(substitute(Q,v,B))[-1]
    need(hb-ha==14*Hq,'homogeneous exchange');need(ha%7==0 and (ha//7-7*Hq)%18==0,'mod18');contexts+=1
 d=json.loads((HERE/'module_census.json').read_text());literal=[]
 for k,r in sorted(d['threshold_30625_symmetric_7p'].items(),key=lambda kv:int(kv[0])):
  G=substitute(r['quotient_rows'],r['port'],r['block_rows']);need(all_h(G)[-1]==int(k),'full prime witness')
  literal.append(dict(H=int(k),p=r['p'],vertices=len(G),rows=G))
 hostile=[]
 for bits in [(1,1,1,1),(1,1,1,0),(1,0,1,1)]:
  R=circ(9,bits);h=all_h(R);E=sum(h[s]*h[511^s] for s in range(512));hostile.append(dict(bits=bits,H=h[-1],C2=E-2*h[-1]))
 need(len({r['C2'] for r in hostile})==3,'order nine counterexample')
 def addmark(R,mask):
  n=len(R);r=list(R)+[mask]
  for u in range(n):
   if not(mask>>u&1):r[u]|=1<<n
  return r
 cut=[all_h(addmark(R,3))[-1] for R in [A,B]];need(cut==[561,579],'nonhomogeneous boundary')
 out=dict(status='PASS',regular_seven=regular,homogeneous_context_checks=contexts,
          higher_cover_content=reduce(gcd,regular[0]['covers'][2:]),actual_7p_witnesses=literal,
          modulo18_cycles=[[1,7,13],[5,17,11]],order9_hostiles=hostile,nonhomogeneous_cut3_counts=cut)
 (HERE/'exchange.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print('PASS',contexts,'contexts',len(literal),'full 7p witnesses')
if __name__=='__main__':main()
