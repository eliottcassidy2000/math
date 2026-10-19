"""Literal coefficient certificate for the degree-42 exceptional class at p=29.
Uses only standard-library polynomial routines from the independent verifier.
No production search or matrix code is imported.
"""
from itertools import permutations
from functools import reduce
from math import gcd
from pathlib import Path
import json
from verify_periods_independent import plus,var,substitution,require
BASE=Path(__file__).resolve().parent
TERMS=[(2,8,32,4),(2,12,28,1),(2,14,26,2),(2,16,24,-13),
(2,18,22,-8),(4,6,32,-12),(4,10,28,11),(4,12,26,1),
(4,14,24,-5),(4,16,22,-10),(4,18,20,-11),(6,8,28,-9),
(6,10,26,-8),(6,12,24,7),(6,14,22,-5),(6,16,20,11),
(8,10,24,8),(8,12,22,14),(8,14,20,12),(8,16,18,2),
(10,12,20,-14),(10,14,18,-1),(12,14,16,-1)]

def polynomial():
 f={}
 for a,b,c,v in TERMS:
  for pm in permutations(range(3)):
   sign=(-1)**sum(pm[i]>pm[j] for i in range(3) for j in range(i+1,3))
   f=plus(f,{tuple((a,b,c)[i] for i in pm):v*sign})
 return f

if __name__=='__main__':
 f=polynomial();v=[var(3,i) for i in range(3)]
 A=[var(3,1,-1),plus(v[0],v[1],-1),v[2]]
 A2=[plus(v[1],v[0],-1),var(3,0,-1),v[2]]
 N=plus(plus(f,substitution(f,A)),substitution(f,A2))
 require(len(f)==138,'term count')
 require(all(sum(t)==42 and min(t)>0 and all(e%2==0 for e in t) for t in f),'support')
 require(all(f.get((t[1],t[2],t[0]))==c for t,c in f.items()),'cyclic')
 require(all(f.get((t[1],t[0],t[2]))==-c for t,c in f.items()),'alternating')
 require(f.get((2,12,28))==1,'nonzero coefficient')
 require(all(c%29==0 for c in N.values()),'norm modulo 29')
 receipt={'status':'PASS','prime':29,'degree':42,'alternants':len(TERMS),
 'monomials':len(f),'cyclic':True,'alternating':True,'all_interior':True,
 'all_exponents_even':True,'literal_norm_zero_mod_29':True,
 'integer_norm_coefficient_gcd':abs(reduce(gcd,N.values()))}
 (BASE/'exception_29_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
