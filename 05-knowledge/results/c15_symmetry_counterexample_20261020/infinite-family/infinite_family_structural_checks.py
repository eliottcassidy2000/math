#!/usr/bin/env python3
"""Bounded structural controls for the infinite affine-Paley obstruction family.
No full graph/tournament census at orders22 or38 is claimed or attempted.
"""
from itertools import combinations
from fractions import Fraction
from math import gcd,factorial
from collections import Counter
import json
checks=0
def ck(p,msg):
 global checks
 checks+=1
 if not p:raise AssertionError(msg)

def prime(p):return p>=2 and all(p%d for d in range(2,int(p**.5)+1))
def paley(p):
 squares={x*x%p for x in range(1,p)}
 return [sum(1<<j for j in range(p) if (j-i)%p in squares) for i in range(p)]
def nontrivial_modules(rows):
 n=len(rows);allmask=(1<<n)-1;found=[]
 for M in range(1,allmask):
  if M.bit_count()<2:continue
  outside=allmask^M;v=outside;good=True
  while v:
   bit=v&-v;i=bit.bit_length()-1;v-=bit;z=rows[i]&M
   if z and z!=M:good=False;break
  if good:found.append(M)
 return found

def pair_orbits(n,gens):
 remaining=set(combinations(range(n),2));out=[]
 while remaining:
  e=min(remaining);O={e};todo=[e]
  while todo:
   u,v=todo.pop()
   for g in gens:
    q=tuple(sorted((g[u],g[v])))
    if q not in O:O.add(q);todo.append(q)
  remaining-=O;out.append(O)
 return out

def affine(n,offset,p,mult,shift):
 g=list(range(n))
 for i in range(p):g[offset+i]=offset+(mult*i+shift)%p
 return g

def primitive(p):
 for r in range(2,p):
  if len({pow(r,j,p) for j in range(p-1)})==p-1:return r
 raise ValueError(p)

def odd_partitions(n,m=None):
 if n==0:yield ();return
 if m is None:m=n
 for x in range(min(n,m),0,-1):
  if x%2:
   for tail in odd_partitions(n-x,x):yield (x,)+tail

def marked_tournaments(k):
 ans=Fraction(0)
 for lam in odd_partitions(k):
  count=Counter(lam);z=1
  for l,m in count.items():z*=l**m*factorial(m)
  ls=lam+(1,1)
  exponent=sum((l-1)//2 for l in ls)+sum(gcd(ls[i],ls[j]) for i in range(len(ls)) for j in range(i+1,len(ls)))
  ans+=Fraction(2**exponent,z)
 ck(ans.denominator==1,'Burnside quotient integrality')
 return int(ans)

def main():
 prime_controls=[]
 for p in (3,7,11,19):
  ck(prime(p) and p%4==3,'prime modulus')
  rows=paley(p);squares={x*x%p for x in range(1,p)}
  ck(all(row.bit_count()==(p-1)//2 for row in rows),'regularity')
  ck(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i,j in combinations(range(p),2)),'tournament')
  mods=nontrivial_modules(rows);ck(not mods,('prime modules',p,mods))
  hgens=[affine(p,0,p,1,1),affine(p,0,p,pow(primitive(p),2,p),0)]
  ck(len(pair_orbits(p,hgens))==1,'unordered 2-homogeneity')
  neg=affine(p,0,p,p-1,0)
  ck(all(((rows[neg[i]]>>neg[j])&1)==1-((rows[i]>>j)&1) for i,j in combinations(range(p),2)),'converse multiplier')
  complete_sign=(-1)**sum(neg[i]>neg[j] for i,j in combinations(range(p),2))
  ck(complete_sign==-1,'complete graph odd normalizer witness')
  prime_controls.append({'p':p,'H_order':p*(p-1)//2,'proper_nontrivial_modules':len(mods),'subsets_tested':2**p-p-2,'pair_orbits':1,'negation_sign_complete_graph':complete_sign})
 # Hostile: regular cyclic order9 does have proper strong modules.
 rows9=[sum(1<<j for j in range(9) if (j-i)%9 in {1,3,4,7}) for i in range(9)]
 mods9=nontrivial_modules(rows9);ck(len(mods9)==3,'composite hostile modules')
 controls=[]
 for a,b in ((7,11),(11,19)):
  k=b-a;n=2*b;ck(a<b<2*a,'root budget')
  N=[];H=[]
  for off,p in ((0,a),(a,b)):
   r=primitive(p)
   H.extend([affine(n,off,p,1,1),affine(n,off,p,r*r%p,0)])
   N.extend([affine(n,off,p,1,1),affine(n,off,p,r,0)])
  for i in range(a+b,n-1):
   g=list(range(n));g[i],g[i+1]=i+1,i;N.append(g)
  O=pair_orbits(n,H);expected=3+2*k+k*(k-1)//2;ck(len(O)==expected,'global H-pair orbits')
  L=list(range(a))+list(range(a+b,n));R=list(range(a,a+b));ls=set(L)
  edge=lambda u,v:(u in ls)!=(v in ls)
  for g in N:
   ck(all(edge(u,v)==edge(g[u],g[v]) for u,v in combinations(range(n),2)),'N fixes Kbb')
   ck(sum(edge(u,v) and g[u]>g[v] for u,v in combinations(range(n),2))%2==0,'N generator even')
  flip=list(range(n))
  for u,v in zip(L,R):flip[u]=v;flip[v]=u
  ck(all(edge(u,v)==edge(flip[u],flip[v]) for u,v in combinations(range(n),2)),'part swap automorphism')
  ck(sum(edge(u,v) and flip[u]>flip[v] for u,v in combinations(range(n),2))%2==1,'part swap odd')
  t=marked_tournaments(k)
  controls.append({'a':a,'b':b,'k':k,'n':n,'H_order':(a*(a-1)//2)*(b*(b-1)//2),'normalizer_order':a*(a-1)*b*(b-1)*factorial(k),'H_pair_orbits':len(O),'tournament_types_proved_by_formula':t,'even_graph_types_proved_by_exact_gap':t-1,'Kbb_N_generator_signs':'all +1','Kbb_full_part_swap_sign':-1})
 print(json.dumps({'status':'PASS','checks':checks,'prime_block_controls':prime_controls,'composite_order9_hostile_nontrivial_modules':[ [i for i in range(9) if M>>i&1] for M in mods9],'family_structural_controls':controls,'scope':'structural checks only; no full census at n22 or n38'},indent=2,sort_keys=True))
if __name__=='__main__':main()
