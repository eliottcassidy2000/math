#!/usr/bin/env python3
"""Independent cycle-partition Burnside check of the two large source counts."""
from fractions import Fraction
from math import factorial,gcd
from collections import Counter
from pathlib import Path
import json

def oddparts(n,m=1):
 if n==0:yield ();return
 for a in range(m,n+1,2):
  for b in oddparts(n-a,a):yield(a,)+b

def b(n,c):
 total=Fraction(0);terms=[]
 for p in oddparts(n):
  z=1
  for a,m in Counter(p).items():z*=a**m*factorial(m)
  e=sum((a-1)//2 for a in p)+sum(gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
  term=Fraction(c**len(p)*2**e,z);total+=term;terms.append({'partition':p,'numerator':term.numerator,'denominator':term.denominator})
 if total.denominator!=1:raise RuntimeError('nonintegral orbit count')
 return total.numerator,terms

def main():
 a,at=b(5,4);c,ct=b(6,2)
 if (2*a,c)!=(17824,3040):raise RuntimeError('Burnside mismatch')
 Path(__file__).with_suffix('.json').write_text(json.dumps({'status':'PASS','C15_relative_orbits':2*a,'F21_relative_orbits':c,'four_colored_tournaments_on_five':a,'C15_terms_before_factor_two':at,'F21_terms':ct},indent=2)+'\n')
if __name__=='__main__':main()
