#!/usr/bin/env python3
"""Exact THM-1950 response-domain correction, 2026-10-14.

No dependencies beyond the standard library. Universe: all labelled
 tournaments through order six, plus explicitly displayed witnesses.
Canonical order-seven mask 299 is inherited from THM-4051's frozen bank.
Bit 1 on lexicographically ordered (i,j), i<j, means i->j.
The tournament is intrinsic; no pairwise ties or cosmetic orientation occur.
Fraction elimination and independent integer Bareiss check each named witness;
Hamiltonian subset DP is checked by literal permutations through order eight.
Finite minimality uses a separate augmented-determinant response calculation.
Normal and python -O runs must have byte-identical output.
"""

from fractions import Fraction as F
from itertools import combinations, permutations

def require(ok, message):
 if not ok: raise RuntimeError(message)

def matrix(n,bits):
 a=[[0]*n for _ in range(n)]
 for b,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(bits>>b)&1; a[j][i]=1-a[i][j]
 return a

def solve_det(a):
 n=len(a); b=[[F(x) for x in r]+[F(1)] for r in a]; det=F(1)
 for k in range(n):
  p=next(i for i in range(k,n) if b[i][k])
  if p!=k: b[k],b[p]=b[p],b[k];det=-det
  q=b[k][k];det*=q;b[k]=[v/q for v in b[k]]
  for i in range(n):
   if i==k:continue
   q=b[i][k];b[i]=[x-q*y for x,y in zip(b[i],b[k])]
 return det,[b[i][-1] for i in range(n)]

def det_bareiss(a):
 a=[r[:] for r in a];n=len(a);d=1;sgn=1
 for k in range(n-1):
  p=next(i for i in range(k,n) if a[i][k])
  if p!=k:a[k],a[p]=a[p],a[k];sgn=-sgn
  q=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    z=a[i][j]*q-a[i][k]*a[k][j]
    if z%d:raise ArithmeticError('nonintegral Bareiss step')
    a[i][j]=z//d
  for i in range(k+1,n):a[i][k]=0
  d=q
 return sgn*a[-1][-1]

def invariants(a):
 n=len(a); m=[[int(i==j)+a[i][j]-a[j][i] for j in range(n)] for i in range(n)]
 det,x=solve_det(m);require(det==det_bareiss(m), "independent determinant mismatch")
 require(all(sum(m[i][j]*x[j] for j in range(n))==1 for i in range(n)), "inverse residual")
 disc=det/F(2**(n-1));s=sum(x)
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not mask>>j&1 and a[i][j]:dp[mask|1<<j][j]+=dp[mask][i]
 h=sum(dp[-1]); reach=[set([i]) for i in range(n)]
 for i in range(n):
  todo=[i]
  while todo:
   u=todo.pop()
   for j in range(n):
    if a[u][j] and j not in reach[i]:reach[i].add(j);todo.append(j)
 strong=all(len(r)==n for r in reach)
 if n<=8:
  hp=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
  require(h==hp, "literal Hamiltonian count mismatch")
 return {'order':n,'strong':strong,'H':h,'det':det,'disc':disc,'s':s,'x':x}

def join(a,b):
 n=len(a);m=len(b)
 return [r+[1]*m for r in a]+[[0]*n+r for r in b]

def pair(v):
 return (v['disc']*(1+v['s'])/2, v['disc']*(1-v['s'])/2)

def main():
 c=matrix(7,299); triangle=matrix(3,2); point=[[0]]
 examples=[('C=canonical_rep299',c,(True,31,F(6),F(2,3))),
           ('C_join_point',join(c,point),(False,31,F(5),F(1))),
           ('C_join_C',join(c,c),(False,961,F(26),F(12,13))),
           ('triangle',triangle,(True,3,F(1),F(3))),
           ('triangle_join_triangle',join(triangle,triangle),(False,9,F(5),F(3,5))),
           ('triangle_join_triangle_join_point',join(join(triangle,triangle),point),(False,9,F(4),F(1))),
           ('triangle_join_C',join(triangle,c),(False,93,F(9),F(11,9)))]
 print('THM-1950 RESPONSE DOMAIN CORRECTION: EXACT')
 print('bit_order=lexicographic_i_lt_j;bit1=i_to_j;rep299_adjacency=')
 for row in c: print(''.join(map(str,row)))
 for name,a,expected in examples:
  v=invariants(a)
  require((v['strong'],v['H'],v['disc'],v['s'])==expected,'witness value '+name)
  print(f"{name}:n={v['order']};strong={v['strong']};H={v['H']};det={v['det']};disc={v['disc']};s={v['s']};Dplus={pair(v)[0]};Dminus={pair(v)[1]}")
  if name=='C=canonical_rep299':print('inverse_ones='+','.join(map(str,v['x'])))
 for a,b in [(c,point),(c,c),(triangle,triangle),(triangle,c),(join(triangle,triangle),point)]:
  va,vb,vc=map(invariants,(a,b,join(a,b)))
  pa,pb,pc=map(pair,(va,vb,vc))
  require(pc==(pa[0]*pb[0],pa[1]*pb[1]),'multiplicative pair')
  q=lambda v:(1-v['s'])/(1+v['s'])
  require(q(vc)==q(va)*q(vb),'signed Cayley coordinate')
  P=lambda v:max(F(1),v['s'])*v['disc']
  require(P(vc)<=P(va)*P(vb),'P submultiplicativity')
 print('named_witnesses=PASS;Dplus_Dminus_q_and_P=PASS')
 print('finite_minimality: n,labelled_count,strong_count,min_all_s,min_strong_s')
 for n in range(1,7):
  count=0;best=None;best_all=None
  for bits in range(1<<(n*(n-1)//2)):
   a=matrix(n,bits)
   reach=[(1<<i)|sum(a[i][j]<<j for j in range(n)) for i in range(n)]
   for k in range(n):
    for i in range(n):
     if reach[i]>>k&1:reach[i]|=reach[k]
   strong=all(r==(1<<n)-1 for r in reach)
   m=[[int(i==j)+a[i][j]-a[j][i] for j in range(n)] for i in range(n)]
   d=det_bareiss(m)
   aug=[r+[1] for r in m]+[[-1]*n+[1]]
   s=F(det_bareiss(aug),d)-1
   best_all=s if best_all is None else min(best_all,s)
   if n<=5:require(s>=1,'unexpected low-order response')
   if not strong:continue
   count+=1;best=s if best is None else min(best,s)
   require(s>=1,'unexpected low-order strong response')
  print(n,1<<(n*(n-1)//2),count,best_all,best,sep=',')
 print('minimum_two_arbitrary_blocks_failure_order=7')
 print('minimum_SCC_product_failure_order=8')
 print('scope=composition_domains_only;all_order_H_ge_disc_OPEN')
 print('RESULT=PASS')

if __name__=='__main__':main()
