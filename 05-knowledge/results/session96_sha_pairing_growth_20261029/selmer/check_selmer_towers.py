#!/usr/bin/env python3
"""Exact model checks, not an arithmetic descent computation for a curve."""
from itertools import product
from math import prod
from fractions import Fraction
import json
checks=0

def check(ok,context):
 global checks
 checks+=1
 if not ok:raise RuntimeError(context)

def b(n,es,r=0):return r*n+sum(min(n,e) for e in es)

finite_groups=0;literal_elements=0
for p in (2,3,5):
 for length in range(4):
  for es in product(range(1,4),repeat=length):
   if p**sum(es)>1500:continue
   mods=[p**e for e in es];elements=list(product(*[range(m) for m in mods]))
   finite_groups+=1;literal_elements+=len(elements)
   counts=[]
   for n in range(6):
    actual=sum(all(p**n*x%m==0 for x,m in zip(xs,mods)) for xs in elements)
    check(actual==p**b(n,es),('literal torsion',p,es,n))
    counts.append(b(n,es))
   inc=[counts[n]-counts[n-1] for n in range(1,6)]
   check(all(inc[i]>=inc[i+1] for i in range(4)),('monotonic',es))
   for n in range(5):
    if counts[n+1]==counts[n]:check(all(e<=n for e in es),('plateau',es,n))
   for e in range(1,4):
    check(inc[e-1]-inc[e]==es.count(e),('second difference',es,e))
# Injections on quotient representatives: each active cyclic factor contributes F_p.
for es in [(),(1,),(2,),(1,3),(2,2,3)]:
 for n in range(1,5):
  active_next=[i for i,e in enumerate(es) if e>=n+1]
  active_now=[i for i,e in enumerate(es) if e>=n]
  check(set(active_next)<=set(active_now),('quotient injection',es,n))
# Eventual count recurrence and exact logarithmic coefficient formula.
for p in (2,3,5):
 for es in [(),(1,),(2,),(1,3),(3,3),(1,2,4)]:
  for r in range(4):
   E=max(es,default=0)
   for n in range(10):
    check(b(n,es,r)==r*n+sum(es.count(e)*min(n,e) for e in set(es)),('profile',es,r,n))
    # Coefficient of [rz + sum m_e z(1-z^e)]/(1-z)^2.
    gf=r*n+sum(es.count(e)*(n-max(0,n-e)) for e in set(es))
    check(gf==b(n,es,r),('generating function',es,r,n))
    if n>=E:check(p**b(n+1,es,r)==p**r*p**b(n,es,r),('tail ratio',es,r,n))
# Selmer-ratio certificate across exact model towers and deliberately slack bounds.
certificate_hits=0
for p in (2,3):
 for R in range(4):
  for tes in [(),(1,),(2,),(1,2)]:
   for ges in [(),(1,),(2,),(1,1),(3,3)]:
    for rdiv in range(3):
     for n in range(5):
      sn=p**(R*n+b(n,tes)+b(n,ges,rdiv))
      snext=p**(R*(n+1)+b(n+1,tes)+b(n+1,ges,rdiv))
      tau=p**(b(n+1,tes)-b(n,tes))
      for r0 in range(R+1):
       for lower_loss,upper_slack in product(range(2),repeat=2):
        L=Fraction(sn,p**lower_loss);U=snext*p**upper_slack
        condition=U<=L*p**r0*tau
        if condition:
         certificate_hits+=1
         check(R==r0 and rdiv==0 and all(e<=n for e in ges),('certificate',p,R,tes,ges,rdiv,n,r0,lower_loss,upper_slack))
       exact=snext<=sn*p**r0*tau
       check(exact==(R==r0 and rdiv==0 and all(e<=n for e in ges)),('exact iff',p,R,tes,ges,rdiv,n,r0))
check(certificate_hits>0,'positive certificate coverage')
# Wrong-rank hostile: actual R=0 and divisible corank1 gives fake rank1 residual.
for p in (2,3,5):
 for n in range(5):check(p**(n+1)==p**n*p,('wrong rank hostile',p,n))
# Natural inclusions scale each perfect alternating pairing by p.
pairing_pairs=0
for p,n in [(2,1),(2,2),(2,3),(3,1),(3,2),(5,1)]:
 q=p**n;qn=p*q
 for x1,x2,y1,y2 in product(range(q),repeat=4):
  det=x1*y2-x2*y1
  upper_num=((p*x1)*(p*y2)-(p*x2)*(p*y1))%qn
  lower_scaled_num=p*p*(det%q)
  check((upper_num-lower_scaled_num)%qn==0,('pairing scale',p,n))
  pairing_pairs+=1
 for u in range(1,q):
  if u%p==0:continue
  check(all((p*v-u)%q for v in range(1,qn) if v%p),('no unit repair',p,n,u))
 check((p*p-p)%qn!=0,('not isometric on basis',p,n))
# Restricted radical computed using actual standard-basis tests in ambient H_E.
radical_cases=0
for p in (2,3,5):
 for E in range(1,5):
  for n in range(1,E+1):
   q=p**n;ambient=p**E;scale=p**(E-n)
   if q*q>20000:continue
   radical=0
   for x,y in product(range(q),repeat=2):
    if (scale*scale*x)%ambient==0 and (scale*scale*y)%ambient==0:radical+=1
   check(radical==p**(2*min(n,E-n)),('radical',p,E,n,radical))
   check((radical==q*q)==(E>=2*n),('zero restriction',p,E,n))
   radical_cases+=1
# Any fixed cutoff is shared by a finite large block and the zero-paired divisible model.
for p in (2,3,5):
 for N in range(1,6):
  E=2*N
  for n in range(N+1):
   check(b(n,(E,E))==2*n,('same finite prefix',p,N,n))
   if n:check(p**(2*(E-n))%p**E==0,('same zero pairing',p,N,n))
print(json.dumps({'status':'PASS','checks':checks,'finite_groups':finite_groups,'literal_elements':literal_elements,'certificate_hits':certificate_hits,'pairing_pairs':pairing_pairs,'radical_cases':radical_cases,'arithmetic_curve_descent_performed':False},sort_keys=True))
