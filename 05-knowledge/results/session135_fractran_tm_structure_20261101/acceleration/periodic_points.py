"""Exact bounded-period point densities from a finite trace-preserving cap."""
from fractions import Fraction
from itertools import product
from math import prod
import json
from accelerate import Program,PRIMEGAME

def need(x,m):
 if not x:raise RuntimeError(m)
def least_return(p,x,L):
 start=tuple(x);cur=list(x)
 for t in range(1,L+1):
  step=p.step(cur)
  if step is None:return None
  _,cur=step
  if tuple(cur)==start:return t
 return None

def density(p,L):
 if type(L)!=int or L<1:raise ValueError('positive period cutoff required')
 caps=[L*max((b[j] for b in p.need),default=0) for j in range(len(p.primes))]
 den=prod(q**c for q,c in zip(p.primes,caps));numerators={};counts={};total=0
 weights=[[((q-1)*q**(c-e-1) if e<c else 1) for e in range(c+1)] for q,c in zip(p.primes,caps)]
 for x in product(*(range(c+1) for c in caps)):
  total+=1;period=least_return(p,x,L)
  if period is not None:
   w=prod(ws[e] for ws,e in zip(weights,x));numerators[period]=numerators.get(period,0)+w;counts[period]=counts.get(period,0)+1
 return {'period_cutoff':L,'caps':dict(zip(p.primes,caps)),'capped_states':total,
         'representative_counts':counts,'density_by_least_period':{k:str(Fraction(v,den)) for k,v in sorted(numerators.items())},
         'density_period_at_most_L':str(Fraction(sum(numerators.values()),den))}

def observations(p,x,T):
 result=[];x=list(x)
 for _ in range(T):
  st=p.step(x)
  if st is None:result.append('halt');break
  i,x=st;result.append(i)
 return result

def main():
 horizon=0
 for rows in [[(1,2)],[(3,2),(2,3)],[(1,10),(5,2)],[(2,1)],[(1,1)],[]]:
  p=Program(rows)
  for T in range(1,6):
   caps=[T*max((b[j] for b in p.need),default=0) for j in range(len(p.primes))]
   for x in product(*(range(2*c+3) for c in caps)):
    capped=[min(v,c) for v,c in zip(x,caps)]
    need(observations(p,x,T)==observations(p,capped,T),'finite horizon cap')
    need(least_return(p,x,T)==least_return(p,capped,T),'return cap');horizon+=1
 toy=density(Program([(3,2),(2,3)]),2);need(toy['density_period_at_most_L']=='5/12','toy density')
 reverse=density(Program([(2,3),(3,2)]),2);need(reverse['density_period_at_most_L']=='5/9','priority changes point density')
 fixed=density(Program([(1,1)]),3);need(fixed['density_period_at_most_L']=='1','fixed points')
 drift=density(Program([(2,1)]),3);need(drift['density_period_at_most_L']=='0','drift not periodic')
 prime=density(Program(PRIMEGAME),2);need(prime['density_period_at_most_L']=='2654208/36409529','inherited PRIMEGAME period2 density')
 # Final halt observation at time T requires a larger cap; the first T means times0..T-1.
 p=Program([(1,2)]);need(observations(p,[2],3)!=observations(p,[3],3),'horizon endpoint hostile')
 return {'status':'PASS','finite_horizon_and_return_checks':horizon,'toy':toy,'reversed_priority_toy':reverse,
  'fixed':fixed,'drift':drift,'primegame_period2':prime,
  'scope':'Periodic POINTS with least period <= specified finite L; not eventual basins or a universal bound on periods.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
