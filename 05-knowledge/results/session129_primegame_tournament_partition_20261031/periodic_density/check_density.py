#!/usr/bin/env python3
"""Exact density consumer and independent integer-level cycle controls."""
from fractions import Fraction as Q
from math import gcd,prod
import json
F=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
D2=(3,7,13,17,19,23,29);D6=(3,5,7,17,29)
B2=(11,13);B6=(4301,17765,16302,9867,8671,23023)
m2=prod(D2);m6=prod(D6)
def require(c,s):
 if not c:raise RuntimeError(s)
def step(n):
 for a,b in F:
  if n%b==0:return n//b*a
 raise RuntimeError('unexpected halt')
def permitted(n,period):return gcd(n,m2)==1 if period==2 else gcd(n,m6)==1 and n%38!=0
def codes(n):return [(p,b,n//b) for p,bs in ((2,B2),(6,B6)) for b in bs if n%b==0 and permitted(n//b,p)]
def scale_count(y,ps,forbid38=False):
 total=0
 for mask in range(1<<len(ps)):
  d=prod(ps[i] for i in range(len(ps)) if mask>>i&1)
  total+=(-1)**mask.bit_count()*(y//d-(y//(38*d) if forbid38 else 0))
 return total
rho2=prod((Q(p-1,p) for p in D2));rho6=prod((Q(p-1,p) for p in D6))*Q(37,38)
d2=rho2*sum((Q(1,b) for b in B2));d6=rho6*sum((Q(1,b) for b in B6))
counts={2:0,6:0};limit=100000
for n in range(1,limit+1):
 enc=codes(n);require(len(enc)<=1,'overlapping phase families')
 z=n
 for j in range(6):z=step(z)
 require(bool(enc)==(z==n),'periodic classification mismatch')
 if enc:
  p,b,x=enc[0];z=n
  for j in range(1,p+1):
   z=step(z);require((z==n)==(j==p),'least period mismatch')
  counts[p]+=1
formula2=sum(scale_count(limit//b,D2) for b in B2)
formula6=sum(scale_count(limit//b,D6,True) for b in B6)
require(counts=={2:formula2,6:formula6},'inclusion exclusion count')
print(json.dumps({'status':'PASS','checked_positive_integers':limit,'periodic_point_counts':counts,'scale_moduli':{'period2':m2,'period6_coprime':m6,'period6_forbidden_divisor':38},'scale_densities':{'period2':str(rho2),'period6':str(rho6)},'point_densities':{'period2':str(d2),'period6':str(d6),'total':str(d2+d6)},'base_cycles':{'period2':B2,'period6':B6}},indent=2,sort_keys=True))
