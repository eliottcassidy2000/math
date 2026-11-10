"""Small exact sanity probe of NEW arithmetic gluing only, not inherited banks."""
from fractions import Fraction
from itertools import product

def test(p,alpha,beta,seeds,levels=4):
    seams=[alpha*b<=beta*a for a,b in zip(seeds,seeds[1:])]
    seams.append(alpha*p*seeds[0]<=beta*seeds[-1])
    intervals=sorted((alpha*p**t*a,beta*p**t*a) for t in range(levels) for a in seeds)
    predicted=all(seams)
    # Exact finite union test without enumerating the (potentially long) ranges.
    observed=all(nxt[0]<=cur[1] for cur,nxt in zip(intervals,intervals[1:]))
    if predicted != observed: raise RuntimeError((p,alpha,beta,seeds,seams))
    return intervals

narrow=test(3,3,5,[237,343,495])
wide=test(3,3,9,[365])
count=0
for p in [3,5,7]:
 for alpha,beta in [(1,3),(3,5),(3,9)]:
  for a in range(1,12,2):
   for b in range(a+2,p*a,2):
    test(p,alpha,beta,[a,b]);count+=1
if not Fraction(5,3)**2 < 3 < Fraction(5,3)**3: raise RuntimeError('phase bound')
# Critical wide rail: one +2 error creates exactly three missing odds.
h=365
old_upper=9*h
new_lower=3*(3*h+2)
if list(range(old_upper,new_lower,2)) != [3285,3287,3289]: raise RuntimeError('seam debt')
print('PASS: exact interval union agrees with cyclic seam criterion on',count,'small cases')
print('NARROW_FIRST_WINDOWS',narrow[:4])
print('WIDE_FIRST_WINDOWS',wide[:3])
print('NARROW_MIN_PHASES',3)
print('WIDE_PLUS_TWO_SEAM_HOLES',list(range(old_upper,new_lower,2)))
