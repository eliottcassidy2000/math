#!/usr/bin/env python3
"""Exact controls for uniform single-cycle D4 bounds and actual Euler obstructions.
Finite pair banks do not replace the analytic all-length theorem.
"""
from itertools import combinations,permutations
import json,hashlib
GATES=0

def check(value,label):
 global GATES
 GATES+=1
 if not value:raise RuntimeError(label)
def mul(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inv(p):return tuple(p.index(i) for i in range(len(p)))
def conj(p,q):return mul(mul(p,q),inv(p))
def support(p):return set(i for i,j in enumerate(p) if i!=j)
def fixed(p):return set(range(len(p)))-support(p)
def braided(p,q):return mul(mul(p,q),p)==mul(mul(q,p),q)
def commute(p,q):return mul(p,q)==mul(q,p)
def perm(n,*cs):
 p=list(range(n))
 for c in cs:
  for i,j in zip(c,c[1:]+c[:1]):p[i]=j
 return tuple(p)
def all_cycles(n,m):
 for cs in combinations(range(n),m):
  for tail in permutations(cs[1:]):yield perm(n,(cs[0],)+tail)
def group(gens):
 identity=tuple(range(len(gens[0])));seen={identity};todo=[identity]
 for p in todo:
  for g in gens:
   q=mul(p,g)
   if q not in seen:seen.add(q);todo.append(q)
 return seen
WORDS=[[-1,-2,1,1,1,2,1],[2,1,3,3,3,-1,-2],[2,1,1,1,-2],
       [2,-1,2,2,1,-2],[2,1,2,3,3,-2,-1,-2],[2,3,2,2,-3,-2]]
def action(word,row):
 row=list(row)
 for letter in word:
  k=abs(letter)-1;p,q=row[k:k+2]
  row[k:k+2]=[conj(p,q),p] if letter>0 else[q,conj(inv(q),p)]
 return row

def pair_banks():
 rows=[]
 for m in [2,3,4,5]:
  n=2*m;sigma=perm(n,tuple(range(m)));counts={};total=0;commuting=0
  for tau in all_cycles(n,m):
   total+=1;j=len(support(sigma)&support(tau))
   if commute(sigma,tau):
    commuting+=1;check(j in [0,m],'commuting single cycles have disjoint or equal supports')
   if not braided(sigma,tau):continue
   counts[j]=counts.get(j,0)+1
   A=fixed(sigma);B=fixed(tau);g=mul(sigma,tau)
   check({g[i] for i in A}==B,'full fixed sets obey actual abstract reaccess')
   n0=len(A&B);check(2*n0>=3*len(A)-n,'abstract cusp injection control')
   check(2*j>=m,'braided single-cycle support overlap')
  rows.append(dict(m=m,ambient=n,total_cycles=total,braid_overlap_counts=counts,commuting_cycles=commuting))
 return rows

def named(n,gens,expected_order,expected_chi):
 a,b,c,d=gens;e=conj(b,c);f=conj(inv(e),a)
 for w in WORDS:check(action(w,gens)==list(gens),'named full six-word control')
 GG=group(gens);check(len(GG)==expected_order,'named exact image order')
 check({g[0] for g in GG}==set(range(n)),'named action really transitive')
 cusp=[(b,c),(b,d),(a,e)];node=[(f,b),(a,d),(c,d)]
 for p,q in cusp:check(braided(p,q),'named cusp braid')
 for p,q in node:check(commute(p,q),'named node commutation')
 check(conj(e,f)==a and conj(e,b)==c,'first node simultaneous conjugation')
 check(conj(c,a)==a and conj(c,e)==b,'third cusp simultaneous conjugation')
 k=len(fixed(b));nc=[len(fixed(p)&fixed(q)) for p,q in cusp];ww=[len(support(p)&support(q)) for p,q in node]
 for p,q in cusp:check({mul(p,q)[i] for i in fixed(p)}==fixed(q),'named full-fixed actual subset control')
 chi=-2*k+sum(nc)+sum(ww);check(chi==expected_chi,'named exact Euler value')
 return dict(degree=n,kept=k,cycle_length=len(support(b)),group_order=len(GG),cusp_counts=nc,node_overlaps=ww,chi=chi)

def mixed_node_bank():
 n=9;sigma=perm(n,(0,1,2),(3,4));counts={};total=0
 for cs in combinations(range(n),3):
  for tail in permutations(cs[1:]):
   for trans in combinations(sorted(set(range(n))-set(cs)),2):
    tau=perm(n,(cs[0],)+tail,trans);total+=1
    if not commute(sigma,tau):continue
    j=len(support(sigma)&support(tau));counts[j]=counts.get(j,0)+1
    check(j in [0,2,3,5],'commuting 32 supports are unions of the two nontrivial cycles')
    check(j!=1,'positive actual node overlap cannot be one')
 check(total==2520,'complete degree9 mixed-cycle universe')
 check(2 in counts and 3 in counts,'both primitive positive overlap controls retained')
 return dict(degree=n,total=total,commuting_overlap_counts=counts)

def main():
 pairs=pair_banks()
 s4=named(4,[perm(4,(2,3)),perm(4,(1,2)),perm(4,(0,1)),perm(4,(2,3))],24,1)
 a4=named(4,[perm(4,(0,2,3)),perm(4,(0,1,2)),perm(4,(0,2,3)),perm(4,(0,2,3))],12,7)
 c4=perm(8,(0,4,2,5));single4=named(8,[c4,perm(8,(0,1,2,3)),c4,perm(8,(1,6,3,7))],192,2)
 mixed=mixed_node_bank()
 data=dict(pair_banks=pairs,sharp_transposition_ledger=s4,sharp_three_cycle_support=a4,sharp_four_cycle_support=single4,mixed32_node_bank=mixed)
 print('FINITE-EXACT SINGLE-CYCLE SUPPORT AND ACTUAL EULER CONTROLS PASS')
 print('Uniform all-cycle proof is analytic; mapping degree>=10 additionally needs the independent involution supplier')
 print(json.dumps(data,sort_keys=True,indent=2))
 print('Always-active gates:',GATES)
 print('Semantic SHA256:',hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest())
if __name__=='__main__':main()
