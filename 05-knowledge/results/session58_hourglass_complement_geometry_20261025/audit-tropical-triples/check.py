#!/usr/bin/env python3
"""Independent row-event-state DP. No candidate implementation imports."""
from itertools import product,combinations
from collections import defaultdict
import json

def need(x,m):
 if not x:raise RuntimeError(m)
def formula(k):return (21**k-6*7**k-3*5**k+14*3**k-6)//48
pairs=list(combinations(range(3),2))
def mask(row):
 flags=[x==0 for x in row]
 for i,j in pairs:flags.extend((row[i]==row[j],row[i]==-row[j]))
 return sum((1<<i) for i,b in enumerate(flags) if b)
rows=[r for r in product((-1,0,1),repeat=3) if sum(x!=0 for x in r)!=1]
need(len(rows)==21,'row alphabet')
# Each bit records whether an unwanted event persists across ALL rows so far.
state={(1<<9)-1:1};answers=[]
for k in range(41):
 signed=state.get(0,0)
 need(signed%48==0,'orbit divisibility')
 need(signed==21**k-6*7**k-3*5**k+14*3**k-6,'DP/formula disagreement')
 answers.append(signed//48)
 need(sum(state.values())==21**k,'all row matrices accounted for')
 new=defaultdict(int)
 for old,n in state.items():
  for row in rows:new[old&mask(row)]+=n
 state=new
# Direct matrix enumeration through k=3 is a second independent exact route.
for k in range(4):
 count=0
 for matrix in product(rows,repeat=k):
  cols=[tuple(row[i] for row in matrix) for i in range(3)]
  if all(any(c) for c in cols) and all(cols[i]!=cols[j] and cols[i]!=tuple(-x for x in cols[j]) for i,j in pairs):count+=1
 need(count==48*answers[k],'direct signed enumeration')
# Polynomial denominator and rational generating identity.
q=[1]
for r in (1,3,5,7,21):
 v=[0]*(len(q)+1)
 for i,a in enumerate(q):v[i]+=a;v[i+1]-=r*a
 q=v
need(q==[1,-37,422,-1982,3801,-2205],'denominator')
conv=[sum(q[j]*answers[k-j] for j in range(min(k,5)+1)) for k in range(41)]
need(conv==[0,0,4,2,-126]+[0]*36,'OGF numerator and recurrence')
# Explicit boundaries and hostile controls.
need(answers[:4]==[0,0,4,150],'boundary counts')
need(21-3*5+2==8 and answers[1]==0,'zero-column-only filter insufficient')
need(all((3**k-1)*(3**k-3)*(3**k-5)//48>=answers[k] for k in range(2,41)),'candidate subset')
try:need(False,'negative guard')
except RuntimeError:pass
else:raise RuntimeError('guard lost under optimization')
# General exponent boxes: quotient only by global sign, NOT projectivization.
def general_formula(exponents):
 from math import prod
 a=exponents
 return (prod(12*x*x+8*x+1 for x in a)-3*prod(4*x+1 for x in a)-6*prod(2*x*x+4*x+1 for x in a)+14*prod(2*x+1 for x in a)-6)//48
def box_count(a):
 vectors=[]
 for v in product(*(range(-x,x+1) for x in a)):
  if any(v) and next(x for x in v if x)>0:vectors.append(v)
 return sum(all(sum(abs(v[i])==max(abs(w[i]) for w in triple) for v in triple)>=2 for i in range(len(a))) for triple in combinations(vectors,3))
box_results={}
for a in ((),(1,),(2,),(3,),(1,2),(2,2),(1,1,2),(2,2,2)):
 state={511:1}
 for exponent in a:
  rows_a=[r for r in product(range(-exponent,exponent+1),repeat=3) if sum(abs(x)==max(map(abs,r)) for x in r)>=2]
  need(len(rows_a)==12*exponent**2+8*exponent+1,'general row count')
  new=defaultdict(int)
  for old,n in state.items():
   for row in rows_a:new[old&mask(row)]+=n
  state=new
 count=state.get(0,0)//48
 need(count==general_formula(a)==box_count(a),'general formula / DP / brute agreement')
 box_results[str(a)]=count
need(box_results['(1, 2)']==15 and box_results['(2, 2)']==54,'general box controls')
# Primitive refinement: each retained row attains its coordinate bound twice.
from math import prod,comb
def primitive_formula(a):
 if not a:return 0
 return (prod(24*x-4 for x in a)-3*4**len(a)-6*prod(4*x+2 for x in a)+14*2**len(a))//48
primitive_boxes={}
for a in ((),(1,),(2,),(3,),(1,2),(2,2),(1,1,1),(1,1,2),(2,2,2)):
 vectors=[v for v in product(*(range(-x,x+1) for x in a)) if any(v) and next(x for x in v if x)>0]
 direct=sum(all(sum(abs(v[i])==a[i] for v in triple)>=2 for i in range(len(a))) for triple in combinations(vectors,3))
 need(direct==primitive_formula(a),'primitive formula/brute agreement')
 primitive_boxes[str(a)]=direct
 need(sum(primitive_formula(tuple(x for x in b if x)) for b in product(*(range(x+1) for x in a)))==general_formula(a),'max allocation decomposition')
P=[primitive_formula((1,)*k) for k in range(41)]
need(P[:4]==[0,0,4,138],'primitive boundaries')
need(all(sum(comb(k,j)*P[j] for j in range(k+1))==answers[k] for k in range(41)),'binomial transform')
pq=[1,-32,284,-928,960]
pc=[sum(pq[j]*P[k-j] for j in range(min(k,4)+1)) for k in range(41)]
need(pc==[0,0,4,10,-120]+[0]*36,'primitive OGF and recurrence starting k5')
need(P[4]!=32*P[3]-284*P[2]+928*P[1]-960*P[0],'primitive recurrence cannot start k4')
print(json.dumps({'status':'PASS','DP_k_range':[0,40],'direct_signed_enumeration_k_range':[0,3],'sequence_prefix':answers[:9],'denominator_coefficients':q,'numerator_coefficients':[0,0,4,2,-126],'all_41_exact_DP_formula_comparisons':True,'negative_controls':'PASS','general_box_checks':box_results,'primitive_box_checks':primitive_boxes,'primitive_prefix':P[:9],'binomial_transform_k_range':[0,40]},indent=2))
