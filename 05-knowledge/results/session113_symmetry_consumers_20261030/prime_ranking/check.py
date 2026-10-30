#!/usr/bin/env python3
from prime_rank import PrimeBank
from bank_rank import digits, number, period, least
from literal_ring import Ring,DSU,grouped_sum
from itertools import product
import json,random
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def audit_ring(p,k,discunit):
 B=PrimeBank(p,k);R=Ring(p,k,discunit);q=R.q
 reps=[];offset=[]
 for d,(_,j) in zip(B.depths,B.blocks):
  # Each valuation has h copies, named by its occurrence among that valuation.
  u=1+sum(z==d for z in B.depths[:len(offset)])
  z=(u,0)
  for _ in range(d):z=R.mul(z,(0,1))
  offset.append(len(reps))
  for i in range(p**j):
   sc=pow(1+p,i,q);reps.append((sc*z[0]%q,sc*z[1]%q))
 mapping=[R.lookup[z] for z in reps]
 need(len(set(i for i,e in mapping))==R.bits,'raw representative bijection')
 inverse={old:(new,e) for new,(old,e) in enumerate(mapping)}
 for r,eqs in enumerate(R.cuts,1):
  literal=DSU(R.bits);pred=DSU(R.bits)
  for a,b,e in eqs:
   x,ea=inverse[a];y,eb=inverse[b]
   need(e^ea^eb==0,'common cut origin')
   literal.add(x,y,0)
  for idx,(d,(_,j)) in enumerate(zip(B.depths,B.blocks)):
   if d>=r:continue
   cap=(r-d-1)//2;step=p**cap;off=offset[idx]
   for t in range(p**j):pred.add(off+t,off+t%step,0)
  lp={};pp={}
  for x in range(R.bits):
   a,b=literal.find(x),pred.find(x)
   need(lp.setdefault(a,b)==b and pp.setdefault(b,a)==a,'literal cut iff word-period caps')
 for a,b in ((1,0),(0,1)):
  action=R.action(a,b)
  for idx,((kind,j),off) in enumerate(zip(B.blocks,offset)):
   shift=(a+(p**(k-2)*b if kind=='B' else 0))%p**j
   for t in range(p**j):
    old,e=mapping[off+t];dest,sgn=action[old];new,ee=inverse[dest]
    need(e^sgn^ee==0 and new==off+(t+shift)%p**j,'literal named action in common gauge')
 need(B.prime_total()==grouped_sum(p,k)[0],'DP versus independent prime convolution')
 rng=random.Random(113+p+k);total=B.prime_total()
 ranks=sorted({0,1,total-2,total-1,total//2}|{rng.randrange(total) for _ in range(8)})
 for rank in ranks:
  w=B.unrank(rank)
  need(B.rank(w)==rank,'prime rank inverse')
  need(B.canonical(w)==w,'literal full group lex canonical')
  raw=[0]*R.bits
  for new,(old,e) in enumerate(mapping):raw[old]=B.flat(w)[new]^e
  need(not any(all(raw[a]^raw[b]==e for a,b,e in eqs) for eqs in R.cuts),'ranked bank literally prime')
 try:B.rank([[0]*p**j for _,j in B.blocks])
 except ValueError:pass
 else:raise RuntimeError('composite bank accepted')
 return dict(p=p,k=k,discriminant=p*discunit,coordinates=R.bits,prime_orbits=str(total),rank_controls=len(ranks),literal_cuts=len(R.cuts))
def audit_thresholds():
 # Exhaust every binary word and every threshold at block length5, all strides.
 B=PrimeBank(5,2)
 for j in (0,1):
  for q in range(j+1):
   reps=[];L=5**j;packet=5**q
   for w in product((0,1),repeat=L):
    ws=[number(w[z:z+packet],2) for z in range(0,L,packet)]
    if ws!=least(ws):continue
    reps.append((number(w,2),period(w)))
   for t in range(2**L+1):
    need(B.ge_periods(j,q,t)==[sum(x>=t and e==5**v for x,e in reps) for v in range(j+1)],'all-threshold exact ORIGINAL period counts')
if __name__=='__main__':
 audit_thresholds()
 rows=[audit_ring(p,k,d) for p,k,d in ((5,2,1),(5,2,2),(7,2,1),(5,3,1))]
 print(json.dumps(dict(status='PASS',checks=checks,families=rows),indent=2))
