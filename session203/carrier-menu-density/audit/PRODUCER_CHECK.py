#!/usr/bin/env python3
from fractions import Fraction
import itertools,json

def need(q,why):
 if not q:raise RuntimeError(why)

def step(p,op):
 X,e,V,d=p
 return (X,0,2*V,d)if op=='O'else(X+V-e,d,2*V-d,d)

seeds=[(677,27,1061,49),(759,29,927,29),(1547,101,1691,109),(821,33,1079,33)]
sumrecip=sum((Fraction(1,V-d)for X,e,V,d in seeds),Fraction())
need(sumrecip==Fraction(693609731,187977134884),'reciprocal sum')
need(4*sumrecip<Fraction(3,200),'relative density bound')
trials=0
for seed in seeds:
 X0,e0,V0,d=seed;A=V0-d
 for n in range(1,13):
  for w in itertools.product('OZ',repeat=n-1):
   word=w+('Z',);p=seed
   for k,op in enumerate(word):
    p=step(p,op)
    need(p[2]>=2**(k+1)*A+d,'V lower bound')
    need(0<=p[1]<=d and p[0]>=X0,'physical row bounds')
   need(p[0]>=X0+2**(n-1)*A,'last-Z growth')
   trials+=1
print(json.dumps({'status':'PASS','last_Z_words_checked':trials,'max_word_length':12,'upper_all_integer_density':str(2*sumrecip),'upper_odd_relative_density':str(4*sumrecip),'scope':'Only the four selected paired O/Z carrier orbits; F(8,137) remains open.'},sort_keys=True,indent=2))
