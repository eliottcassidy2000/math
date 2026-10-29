#!/usr/bin/env python3
"""Exact independent controls; no producer imports; assertions survive -O."""
from fractions import Fraction as Q
from itertools import product, permutations
from math import gcd
import json
checks=0

def check(p):
 global checks
 checks+=1
 if not p: raise RuntimeError('check failed')

def compose(A,B): # chronological A then B
 u,v=A;s,t=B
 return s*u,s*v+t

def chart(b,w):
 A=(Q(1),Q(0))
 for a in w: A=compose(A,(Q(b,2**a),Q(1,2**a)))
 return A

def independent(b,w):
 H=sum(w);k=len(w); h=0; C=0
 for i,a in enumerate(w): C+=b**(k-1-i)*2**h;h+=a
 return Q(b**k,2**H),Q(C,2**H)

def rank(A):u,v=A;return (u-1)/v

def cycle(A):u,v=A;return v/(1-u)
def chain(seq):
 out=(Q(1),Q(0))
 for A in seq:out=compose(out,A)
 return out

words=[w for k in range(1,4) for w in product(range(1,5),repeat=k)]
pairs=0; triples=0
for b in (3,5,7,9):
 bank=[chart(b,w) for w in words]
 for w,A in zip(words,bank):check(A==independent(b,w))
 for A,B in product(bank,repeat=2):
  AB=compose(A,B);BA=compose(B,A); r=rank(A);s=rank(B)
  check(AB[0]==BA[0]);check(AB[1]-BA[1]==A[1]*B[1]*(s-r))
  check((AB==BA)==(r==s))
  check(min(r,s)<=rank(AB)<=max(r,s));pairs+=1
 # all triples from charts of lengths 1 and 2 on digits 1,2,3
 small=[chart(b,w) for k in (1,2) for w in product(range(1,4),repeat=k)]
 for seq in product(small,repeat=3):
  values=[chain(p)[1] for p in permutations(seq)]
  check(chain(sorted(seq,key=rank))[1]==max(values))
  check(chain(sorted(seq,key=rank,reverse=True))[1]==min(values));triples+=1

focal=[]
for w in ((1,1,5),(1,3,3)):
 A=chart(3,w);B=chart(3,w+(1,));x=cycle(A);y=cycle(B)
 check(x.denominator==101)
 focal.append({'word':w,'rank':str(rank(A)),'cycle':str(x),'appended_cycle':str(y),'appended_denominator':y.denominator})
check([r['appended_denominator'] for r in focal]==[35,175])
A=chart(3,(2,));AA=compose(A,A);B=chart(3,(1,))
check(rank(A)==rank(AA)==-1);check(compose(A,AA)==compose(AA,A))
left=chain([A,B,AA]);right=chain([AA,B,A])
check(left==(Q(81,128),Q(119,128)))
check(right==(Q(81,128),Q(143,128)))
check(left[1]-right[1]==A[0]*A[1]*B[1]*(rank(A)-rank(B)))
append4=chart(3,(4,))
check(cycle(compose(left,append4))==Q(97,361));check(cycle(compose(right,append4))==Q(557,1805))
check(cycle(compose(A,B))==Q(-7));check(cycle(compose(AA,B))==Q(37,5))
print(json.dumps({'checks':checks,'pair_cases':pairs,'triple_permutation_cases':triples,'multipliers':[3,5,7,9],'word_digits':[1,2,3,4],'word_lengths':[1,2,3],'focal':focal,'nonadjacent_tie':{'left':[str(x) for x in left],'right':[str(x) for x in right],'cycle_left':str(cycle(left)),'cycle_right':str(cycle(right)),'append4_left':str(cycle(compose(left,append4))),'append4_right':str(cycle(compose(right,append4)))},'status':'PASS'},indent=2))
